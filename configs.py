AA_VOCAB = 'ACDEFGHIKLMNPQRSTVWY'
AA_TO_INT = {aa: idx + 1 for idx, aa in enumerate(AA_VOCAB)}
REVERSE_AA_DICT = {idx + 1: aa for idx, aa in enumerate(AA_VOCAB)}
REVERSE_AA_DICT[21] = 'A'  
REVERSE_AA_DICT[0] = ''   

FEATURE_NAMES = [
    'aiupred', 'asaquick', 'RAAP', 'esm2', 'morfchibi',
    'Local_Morf_WeightedEuclid', 'Local_Fusion_WeightedEuclid',
    'Global_Morf_MinEuclid', 'Global_Morf_DistToCore',
    'Global_Morf_BindingRatio', 'Global_Fusion_MinEuclid',
    'Global_Fusion_DistToCore', 'Global_Fusion_BindingRatio'
]

FEATURE_BOUNDS = {
           'aiupred': [0.0026, 1.0371],
            'asaquick': [0.0030, 0.7767],
            'RAAP': [-0.2180, 0.3000],
            'esm2': [0.0000, 1.0000],
            'morfchibi': [0.0000, 0.9908],
            'Local_Morf_WeightedEuclid': None, # [0.0600, 4.0743],
            'Local_Fusion_WeightedEuclid': None, # [2.1702, 4.0729],
            'Global_Morf_MinEuclid': None, #[0.0600, 4.0743],
            'Global_Morf_DistToCore': [0.0000, 1.0000],
            'Global_Morf_BindingRatio': [0.0000, 1.0000],
            'Global_Fusion_MinEuclid': None, # [2.1702, 4.0557],
            'Global_Fusion_DistToCore': [0.0000, 1.0000],
            'Global_Fusion_BindingRatio': [0.0000, 1.0000]
}

RESIDUE_FEAT_NAMES = ['aiupred', 'asaquick', 'RAAP', 'esm2', 'morfchibi']
GLOBAL_FEAT_NAMES = [
    'Local_Morf_WeightedEuclid', 'Local_Fusion_WeightedEuclid',
    'Global_Morf_MinEuclid', 'Global_Morf_DistToCore',
    'Global_Morf_BindingRatio', 'Global_Fusion_MinEuclid',
    'Global_Fusion_DistToCore', 'Global_Fusion_BindingRatio'
]