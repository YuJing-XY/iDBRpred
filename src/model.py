import os
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, models

AA_VOCAB = 'ACDEFGHIKLMNPQRSTVWY'
AA_TO_INT = {aa: idx + 1 for idx, aa in enumerate(AA_VOCAB)}
AA_TO_INT['X'] = 21

def apply_rope(x, seq_len, head_dim, dtype):
    pos = tf.range(seq_len, dtype=tf.float32)
    inv_freq = 1.0 / (10000.0 ** (tf.range(0, head_dim, 2, dtype=tf.float32) / tf.cast(head_dim, tf.float32)))
    freqs = tf.einsum('i,j->ij', pos, inv_freq)
    freqs = tf.concat([freqs, freqs], axis=-1)
    sin = tf.cast(tf.sin(freqs), dtype)
    cos = tf.cast(tf.cos(freqs), dtype)
    sin = tf.reshape(sin, [1, seq_len, 1, head_dim])
    cos = tf.reshape(cos, [1, seq_len, 1, head_dim])
    x1, x2 = tf.split(x, 2, axis=-1)
    return x * cos + tf.concat([-x2, x1], axis=-1) * sin

class RoPEWindowedAttention(layers.Layer):
    def __init__(self, embed_dim, num_heads, window_size=31, **kwargs):
        super().__init__(**kwargs)
        self.num_heads = num_heads
        self.head_dim = embed_dim // num_heads
        self.window_size = window_size
        self.q_proj = layers.Dense(embed_dim)
        self.k_proj = layers.Dense(embed_dim)
        self.v_proj = layers.Dense(embed_dim)
        self.o_proj = layers.Dense(embed_dim)

    def call(self, inputs, mask=None):
        B = tf.shape(inputs)[0]
        L = tf.shape(inputs)[1]
        
        q = tf.reshape(self.q_proj(inputs), [B, L, self.num_heads, self.head_dim])
        k = tf.reshape(self.k_proj(inputs), [B, L, self.num_heads, self.head_dim])
        v = tf.reshape(self.v_proj(inputs), [B, L, self.num_heads, self.head_dim])
        
        q_rope = apply_rope(q, L, self.head_dim, inputs.dtype)
        k_rope = apply_rope(k, L, self.head_dim, inputs.dtype)
        
        pad_size = self.window_size // 2
        frames_k = tf.transpose(tf.signal.frame(tf.pad(k_rope, [[0,0],[pad_size,pad_size],[0,0],[0,0]]), self.window_size, 1, axis=1), [0, 1, 3, 2, 4])
        frames_v = tf.transpose(tf.signal.frame(tf.pad(v, [[0,0],[pad_size,pad_size],[0,0],[0,0]]), self.window_size, 1, axis=1), [0, 1, 3, 2, 4])
        
        scores = tf.matmul(tf.expand_dims(q_rope, 3), frames_k, transpose_b=True) / tf.math.sqrt(tf.cast(self.head_dim, inputs.dtype))
        
        if mask is not None:
            attn_mask = tf.reshape(tf.signal.frame(tf.pad(tf.cast(mask, inputs.dtype), [[0,0],[pad_size,pad_size]]), self.window_size, 1, axis=1), [B, L, 1, 1, self.window_size])
            scores += (1.0 - attn_mask) * tf.cast(-1e4, inputs.dtype)
            
        attn_output = tf.matmul(tf.nn.softmax(scores, axis=-1), frames_v)
        attn_output = tf.reshape(attn_output, [B, L, self.num_heads * self.head_dim])
        return self.o_proj(attn_output)

class SwiGLU(layers.Layer):
    def __init__(self, hidden_dim, out_dim, dropout_rate=0.25, **kwargs):
        super().__init__(**kwargs)
        self.w1 = layers.Dense(hidden_dim)
        self.w2 = layers.Dense(hidden_dim)
        self.w3 = layers.Dense(out_dim)
        self.dropout = layers.Dropout(dropout_rate)

    def call(self, x, training=False):
        return self.w3(self.dropout(tf.nn.silu(self.w1(x)) * self.w2(x), training=training))

class TransformerBlock(layers.Layer):
    def __init__(self, embed_dim, num_heads, ff_dim, window_size=31, dropout_rate=0.25, **kwargs):
        super().__init__(**kwargs)
        self.ln1 = layers.LayerNormalization(epsilon=1e-6)
        self.ln2 = layers.LayerNormalization(epsilon=1e-6)
        self.attn = RoPEWindowedAttention(embed_dim, num_heads, window_size)
        self.ffn = SwiGLU(int(ff_dim * 2 / 3), embed_dim, dropout_rate)
        self.drop1 = layers.Dropout(dropout_rate)
        self.drop2 = layers.Dropout(dropout_rate)

    def call(self, inputs, mask=None, training=False):
        x = inputs + self.drop1(self.attn(self.ln1(inputs), mask=mask), training=training)
        return x + self.drop2(self.ffn(self.ln2(x), training=training), training=training)

def build_idbr_model():
    embed_dim = 64
    num_heads = 2
    window_size = 31
    dropout_rate = 0.25
    total_features = 13

    aa_input = layers.Input(shape=(None,), dtype=tf.int32, name='aa_seq')
    feat_input = layers.Input(shape=(None, total_features), dtype=tf.float32, name='features')
    pad_mask = layers.Lambda(lambda x: tf.not_equal(x, 0))(aa_input)

    residue_feats = layers.Lambda(lambda x: tf.gather(x, [0, 1, 2, 3, 4], axis=-1))(feat_input)
    global_feats = layers.Lambda(lambda x: tf.gather(x, [5, 6, 7, 8, 9, 10, 11, 12], axis=-1))(feat_input)

    x_local = layers.Dense(embed_dim, activation='gelu')(residue_feats)
    x_local = TransformerBlock(embed_dim, num_heads, embed_dim * 4, window_size, dropout_rate)(x_local, mask=pad_mask)

    x_global = layers.Dense(16, activation='gelu')(global_feats)
    
    fused_x = layers.Concatenate(axis=-1)([x_local, x_global])
    fused_x = layers.Dense(64, activation='gelu')(fused_x)
    fused_x = layers.Dropout(dropout_rate)(fused_x)
    fused_x = layers.Dense(32, activation='gelu')(fused_x)
    fused_x = layers.Dropout(dropout_rate)(fused_x)

    output = layers.Dense(1, activation='sigmoid', dtype='float32')(fused_x)
    
    return models.Model(inputs=[aa_input, feat_input], outputs=output)

class IDRPredictor:
    def __init__(self, weights_path):
        self.model = build_idbr_model()
        if os.path.exists(weights_path):
            self.model.load_weights(weights_path)
        else:
            raise FileNotFoundError(f"Model weights not found at: {weights_path}")

    @tf.function(reduce_retracing=True)
    def _fast_predict(self, aa, feat):
        return self.model([aa, feat], training=False)

    def predict(self, sequence, features):
        aa_encoded = [AA_TO_INT.get(aa, 21) for aa in sequence]
        aa_tensor = np.array([aa_encoded], dtype=np.int32)
        feat_tensor = np.expand_dims(features, axis=0)
        
        preds = self._fast_predict(aa_tensor, feat_tensor).numpy()
        return preds[0, :, 0]