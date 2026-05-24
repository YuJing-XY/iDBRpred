import os
import argparse
import yaml
import sys
import time
import numpy as np

import tensorflow as tf

gpus = tf.config.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError as e:
        print(f"Error setting memory growth: {e}")

from src.utils import check_and_download_esm2, parse_fasta
from src.features import FeatureProcessor
from src.model import IDRPredictor
from feature_extractors.extractors import (
    AIUPredExtractor, 
    ASAquickExtractor, 
    ESM2Extractor, 
    MoRFchibiExtractor
)

def load_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def standardize_sequence(seq):
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    seq = seq.upper()
    return ''.join([aa if aa in valid_aa else 'A' for aa in seq])

def main():
    parser = argparse.ArgumentParser(description="iDBRpred: An accurate predictor of intrinsically disordered binding residues in protein sequences")
    parser.add_argument('-i', '--input', required=True, help="Input FASTA file path")
    parser.add_argument('-o', '--outdir', required=True, help="Output directory path")
    parser.add_argument('-c', '--config', default='config.yaml', help="Config file path")
    parser.add_argument('-d', '--device', type=str, choices=['cpu', 'gpu'], help="Target computing device")
    parser.add_argument('-t', '--threads', type=int, help="Number of CPU threads to use")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        sys.exit(f"Error: Config file not found at {args.config}")
    
    config = load_config(args.config)

    raw_device = args.device if args.device else config['system']['default_device']
    device = 'cuda' if raw_device in ['gpu', 'cuda'] else 'cpu'
    
    threads = args.threads if args.threads else config['system']['default_threads']

    base_dir = os.path.dirname(os.path.abspath(__file__))
    
    os.makedirs(args.outdir, exist_ok=True)
    out_pred_file = os.path.join(args.outdir, 'binding_predictions.txt')
    out_time_file = os.path.join(args.outdir, 'times.csv')

    check_and_download_esm2(config['esm2'], base_dir)

    print("Initializing feature extractors...")
    aiupred_ext = AIUPredExtractor(force_cpu=(device=='cpu'))
    asaquick_ext = ASAquickExtractor(base_dir)
    
    esm2_ext = ESM2Extractor(
        base_dir=base_dir, 
        model_name=config['esm2']['model_name'], 
        model_dir=config['esm2'].get('model_dir', ''), 
        device=device, 
        threads=threads
    )
    
    morf_ext = MoRFchibiExtractor(base_dir, config['morfchibi']['mode'], config['databases'], threads)

    print("Loading deep learning model...")
    weights_path = os.path.join(base_dir, 'model_weights', 'Best_B_res_Transformer.keras')
    model = IDRPredictor(weights_path)
    processor = FeatureProcessor()

    sequences = parse_fasta(args.input)
    skipped_proteins = [] 
    
    print(f"Starting prediction. Results will be saved to {args.outdir}")
    
    with open(out_pred_file, 'w') as fout, open(out_time_file, 'w') as ftime:
        ftime.write("Protein_ID,Time_ms\n")
        
        for pid, original_seq in sequences.items():
            print(f"Processing {pid}...")
            
            seq = standardize_sequence(original_seq)
            if seq != original_seq.upper():
                print(f"  -> [Info] Non-standard amino acids detected in {pid}, converted to 'A'.")
            
            if len(seq) < 21:
                reason = f"Length ({len(seq)}) < 21"
                print(f"  -> [Warning] Skipped {pid}: {reason}.")
                skipped_proteins.append((pid, reason))
                continue
            
            start_time = time.perf_counter()
            print("  -> Running AIUPred...")
            f_aiu = aiupred_ext.extract(seq)
            print("  -> Running ASAquick...")
            f_asa = asaquick_ext.extract(seq)
            print("  -> Running ESM2 conservation...")
            f_esm = esm2_ext.extract(seq)
            print("  -> Running MoRFchibi...")
            f_morf = morf_ext.extract(seq)
            
            raw_features = {
                'aiupred': f_aiu,
                'asaquick': f_asa,
                'esm2': f_esm,
                'morfchibi': f_morf
            }
            
            input_features = processor.build_feature_matrix(seq, raw_features)
            predictions = model.predict(seq, input_features)
            
            end_time = time.perf_counter()
            elapsed_ms = (end_time - start_time) * 1000
            
            fout.write(f">{pid}\n")
            fout.write(f"{seq}\n")
            fout.write(" ".join([f"{p:.4f}" for p in predictions]) + "\n")
            
            ftime.write(f"{pid},{elapsed_ms:.2f}\n")

    print("\n" + "="*40)
    print("Prediction Run Summary")
    print("="*40)
    print(f"Total processed: {len(sequences) - len(skipped_proteins)}")
    print(f"Total skipped:   {len(skipped_proteins)}")
    
    if skipped_proteins:
        print("\nSkipped Proteins Details:")
        for pid, reason in skipped_proteins:
            print(f" - {pid}: {reason}")
            
    print(f"\nPrediction complete. Results and profiling saved to {args.outdir}")

if __name__ == "__main__":
    main()
