import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

class FeatureProcessor:
    def __init__(self):
        self.window_size = 31
        self.threshold_a = 1.057653
        self.threshold_b = 3.780319
        self.weights = self._get_linear_weights(self.window_size)
        
        self.raap_map = {
            'A': -0.057, 'C': -0.089, 'D': -0.018, 'E': 0.023, 'F': 0.132,
            'G': -0.154, 'H': -0.046, 'I': 0.103, 'K': 0.162, 'L': 0.065,
            'M': 0.115, 'N': 0.016, 'P': -0.179, 'Q': 0.033, 'R': 0.300,
            'S': -0.218, 'T': -0.119, 'V': -0.052, 'W': 0.249, 'Y': 0.195
        }
        
        self.bounds = {
            'aiupred': (0.26, 1.0),
            'asaquick': (0.0, 0.78),
            'RAAP': (-0.22, 0.3),
            'esm2': (0.0, 1.0),
            'morfchibi': (0.0, 1.0)
        }

    def _get_linear_weights(self, window_size):
        if window_size == 1: 
            return np.array([1.0], dtype=np.float32)
        half_w = window_size // 2
        left = np.linspace(0.1, 1.0, num=half_w + 1)
        return np.concatenate([left, left[:-1][::-1]]).astype(np.float32)

    def _compute_weighted_distance(self, seq, weights):
        half_w = len(weights) // 2
        padded = np.pad(seq, (half_w, half_w), mode='edge')
        windows = sliding_window_view(padded, window_shape=len(weights))
        diff_sq = (windows - 1.0) ** 2
        return np.sqrt(np.sum(diff_sq * weights, axis=1))

    def _normalize(self, seq, feat_name):
        seq = np.array(seq, dtype=np.float64)
        min_v, max_v = self.bounds[feat_name]
        denom = max_v - min_v if (max_v - min_v) != 0 else 1.0
        return np.clip((seq - min_v) / denom, 0.0, 1.0)

    def build_feature_matrix(self, sequence, raw_features):
        seq_len = len(sequence)
        
        raap_seq = np.array([self.raap_map.get(aa, 0.0) for aa in sequence], dtype=np.float32)
        
        aiu_n = self._normalize(raw_features['aiupred'], 'aiupred')
        asa_n = self._normalize(raw_features['asaquick'], 'asaquick')
        raap_n = self._normalize(raap_seq, 'RAAP')
        esm_n = self._normalize(raw_features['esm2'], 'esm2')
        morf_n = self._normalize(raw_features['morfchibi'], 'morfchibi')
        
        feat_a = morf_n
        feat_b = asa_n * raap_n * esm_n
        
        dist_a = self._compute_weighted_distance(feat_a, self.weights)
        dist_b = self._compute_weighted_distance(feat_b, self.weights)
        
        def calc_global(dist_seq, threshold):
            val_max = np.min(dist_seq)
            core_idx = np.argmin(dist_seq)
            val_center = np.abs(np.arange(seq_len) - core_idx) / max(1, seq_len - 1)
            val_content = np.sum(dist_seq <= threshold) / seq_len
            return np.full(seq_len, val_max), val_center, np.full(seq_len, val_content)

        max_a, center_a, content_a = calc_global(dist_a, self.threshold_a)
        max_b, center_b, content_b = calc_global(dist_b, self.threshold_b)
        
        morf_nn_input = np.clip(np.array(raw_features['morfchibi'], dtype=np.float64) / 0.9908, 0.0, 1.0)
        
        features = np.column_stack([
            raw_features['aiupred'],     
            raw_features['asaquick'],    
            raap_seq,                    
            raw_features['esm2'],        
            morf_nn_input,               
            dist_a,
            dist_b,
            max_a,
            center_a,
            content_a,
            max_b,
            center_b,
            content_b
        ])
        
        return features.astype(np.float32)