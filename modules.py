import tensorflow as tf
from tensorflow.keras import layers, models
from configs import FEATURE_NAMES, FEATURE_BOUNDS, RESIDUE_FEAT_NAMES, GLOBAL_FEAT_NAMES

def apply_rope(x, seq_len, head_dim, dtype):
    pos = tf.range(seq_len, dtype=tf.float32)
    inv_freq = 1.0 / (10000.0 ** (tf.range(0, head_dim, 2, dtype=tf.float32) / tf.cast(head_dim, tf.float32)))
    freqs = tf.einsum('i,j->ij', pos, inv_freq)
    freqs = tf.concat([freqs, freqs], axis=-1)
    sin, cos = tf.cast(tf.sin(freqs), dtype), tf.cast(tf.cos(freqs), dtype)
    sin, cos = tf.reshape(sin, [1, seq_len, 1, head_dim]), tf.reshape(cos, [1, seq_len, 1, head_dim])
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
        B, L = tf.shape(inputs)[0], tf.shape(inputs)[1]
        q = tf.reshape(self.q_proj(inputs), [B, L, self.num_heads, self.head_dim])
        k = tf.reshape(self.k_proj(inputs), [B, L, self.num_heads, self.head_dim])
        v = tf.reshape(self.v_proj(inputs), [B, L, self.num_heads, self.head_dim])

        q_rope = apply_rope(q, L, self.head_dim, inputs.dtype)
        k_rope = apply_rope(k, L, self.head_dim, inputs.dtype)
        
        pad_size = self.window_size // 2
        padded_k = tf.pad(k_rope, [[0,0], [pad_size,pad_size], [0,0], [0,0]])
        padded_v = tf.pad(v, [[0,0], [pad_size,pad_size], [0,0], [0,0]])
        
        frames_k = tf.transpose(tf.signal.frame(padded_k, self.window_size, 1, axis=1), [0, 1, 3, 2, 4])
        frames_v = tf.transpose(tf.signal.frame(padded_v, self.window_size, 1, axis=1), [0, 1, 3, 2, 4])
        
        scores = tf.matmul(tf.expand_dims(q_rope, 3), frames_k, transpose_b=True) / tf.math.sqrt(tf.cast(self.head_dim, inputs.dtype))
        
        if mask is not None:
            padded_mask = tf.pad(tf.cast(mask, inputs.dtype), [[0,0], [pad_size,pad_size]])
            attn_mask = tf.reshape(tf.signal.frame(padded_mask, self.window_size, 1, axis=1), [B, L, 1, 1, self.window_size])
            scores += (1.0 - attn_mask) * tf.cast(-1e4, inputs.dtype)
            
        attn_out = tf.matmul(tf.nn.softmax(scores, axis=-1), frames_v)
        return self.o_proj(tf.reshape(attn_out, [B, L, self.num_heads * self.head_dim]))

class SwiGLU(layers.Layer):
    def __init__(self, hidden_dim, out_dim, dropout_rate=0.1, **kwargs):
        super().__init__(**kwargs)
        self.w1 = layers.Dense(hidden_dim)
        self.w2 = layers.Dense(hidden_dim)
        self.w3 = layers.Dense(out_dim)
        self.dropout = layers.Dropout(dropout_rate)

    def call(self, x, training=False):
        return self.w3(self.dropout(tf.nn.silu(self.w1(x)) * self.w2(x), training=training))

class TransformerBlock(layers.Layer):
    def __init__(self, embed_dim, num_heads, ff_dim, window_size=31, dropout_rate=0.1, **kwargs):
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

def build_predictor_model():
    embed_dim = 64
    num_heads = 2
    dropout_rate = 0.25

    aa_input = layers.Input(shape=(None,), dtype=tf.int32, name='aa_seq')
    feat_input = layers.Input(shape=(None, len(FEATURE_NAMES)), dtype=tf.float32, name='num_features')
    pad_mask = layers.Lambda(lambda x: tf.not_equal(x, 0))(aa_input)

    def extract_and_norm(target_feats):
        norm_indices, norm_min, norm_max, pass_indices = [], [], [], []
        for f in target_feats:
            if f in FEATURE_NAMES:
                bounds = FEATURE_BOUNDS.get(f)
                idx = FEATURE_NAMES.index(f)
                if bounds is not None:
                    norm_indices.append(idx)
                    norm_min.append(bounds[0])
                    norm_max.append(bounds[1])
                else:
                    pass_indices.append(idx)
                    
        streams = []
        if norm_indices:
            feat_to_norm = layers.Lambda(lambda x: tf.gather(x, norm_indices, axis=-1))(feat_input)
            min_t = tf.constant(norm_min, tf.float32)
            max_t = tf.constant(norm_max, tf.float32)
            normalized = layers.Lambda(lambda x: tf.clip_by_value((x - tf.cast(min_t, x.dtype)) / (tf.cast(max_t, x.dtype) - tf.cast(min_t, x.dtype) + 1e-7), 0.0, 1.0))(feat_to_norm)
            streams.append(normalized)
            
        if pass_indices:
            streams.append(layers.Lambda(lambda x: tf.gather(x, pass_indices, axis=-1))(feat_input))
            
        return layers.Concatenate(axis=-1)(streams) if len(streams) > 1 else streams[0]

    sliced_residue = extract_and_norm(RESIDUE_FEAT_NAMES)
    x_local = layers.Dense(embed_dim, activation='gelu', name='residue_proj')(sliced_residue)
    x_local = TransformerBlock(x_local.shape[-1], num_heads, x_local.shape[-1]*4, 31, dropout_rate)(x_local, mask=pad_mask)

    sliced_global = extract_and_norm(GLOBAL_FEAT_NAMES)
    x_global = layers.Dense(16, activation='gelu', name='global_proj')(sliced_global)

    fused_x = layers.Concatenate(axis=-1, name='local_global_fusion')([x_local, x_global])

    fused_x = layers.Dense(64, activation='gelu', name='head_dense_64')(fused_x)
    fused_x = layers.Dropout(dropout_rate)(fused_x)
    fused_x = layers.Dense(32, activation='gelu', name='head_dense_32')(fused_x)
    fused_x = layers.Dropout(dropout_rate)(fused_x)

    output = layers.Dense(1, activation='sigmoid', name='b_res_prob', dtype='float32')(fused_x)
    return models.Model(inputs=[aa_input, feat_input], outputs=output)