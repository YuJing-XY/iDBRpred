import os
import argparse
import time
import numpy as np
import pandas as pd
import tensorflow as tf

from configs import AA_TO_INT, REVERSE_AA_DICT, FEATURE_NAMES
from modules import build_predictor_model

def parse_dataset_table(txt_path):
    raw_data = []
    if not os.path.exists(txt_path):
        print(f"[Error] Dataset file not found: {txt_path}")
        return raw_data

    with open(txt_path, 'r') as f:
        curr_pid = None
        curr_aa, curr_feats, curr_bres, curr_dreg = [], [], [], []
        col_map = {}
        feat_indices = []

        for line in f:
            line = line.strip()
            if not line:
                if curr_pid and len(curr_aa) > 0:
                    raw_data.append({
                        'pid': curr_pid,
                        'aa_seq': np.array(curr_aa, dtype=np.int32),
                        'features': np.array(curr_feats, dtype=np.float32),
                        'b_res': np.array(curr_bres, dtype=np.float32),
                        'loss_mask': np.ones(len(curr_aa), dtype=np.float32),
                        'raw_d_reg': curr_dreg
                    })
                    curr_aa, curr_feats, curr_bres, curr_dreg = [], [], [], []
                continue

            if line.startswith('>'):
                curr_pid = line.replace('>', '').strip()
            elif line.startswith('Idx'):
                col_map = {col: i for i, col in enumerate(line.split())}
                feat_indices = [col_map[feat] for feat in FEATURE_NAMES]
            else:
                parts = line.split()
                if len(parts) < len(col_map):
                    continue
                curr_aa.append(AA_TO_INT.get(parts[col_map['AA']], 21))
                curr_bres.append(float(parts[col_map['B_res']]))
                curr_dreg.append(parts[col_map['D_reg']])
                curr_feats.append([float(parts[idx]) for idx in feat_indices])

        # Catch the last sequence
        if curr_pid and len(curr_aa) > 0:
            raw_data.append({
                'pid': curr_pid,
                'aa_seq': np.array(curr_aa, dtype=np.int32),
                'features': np.array(curr_feats, dtype=np.float32),
                'b_res': np.array(curr_bres, dtype=np.float32),
                'loss_mask': np.ones(len(curr_aa), dtype=np.float32),
                'raw_d_reg': curr_dreg
            })
            
    return raw_data


def run_predictor(input_file, weights_file, output_dir, batch_size=64):
    os.makedirs(output_dir, exist_ok=True)
    
    print("[Info] Building model structure and injecting pre-trained weights...")
    model = build_predictor_model()
    
    try:
        model.load_weights(weights_file)
        print("[Info] Model weights loaded successfully.")
    except Exception as e:
        print(f"[Error] Failed to load weights from {weights_file}. Details: {e}")
        return

    print(f"[Info] Parsing input sequence sheet: {os.path.basename(input_file)}")
    parsed_data = parse_dataset_table(input_file)
    if not parsed_data:
        print("[Error] Input file is empty or format is invalid. Terminating inference.")
        return
        
    print(f"[Info] Total proteins loaded: {len(parsed_data)}")

    @tf.function(reduce_retracing=True)
    def graph_predict(aa, feat):
        return model([aa, feat], training=False)

    predictions_pool = []
    timing_records = []

    print("[Info] Commencing batched inference...")
    for i in range(0, len(parsed_data), batch_size):
        batch = parsed_data[i : i+batch_size]
        max_len = max(len(item['aa_seq']) for item in batch)
        
        batch_aa = [np.pad(item['aa_seq'], (0, max_len - len(item['aa_seq']))) for item in batch]
        batch_feat = [np.pad(item['features'], ((0, max_len - len(item['aa_seq'])), (0, 0))) for item in batch]

        start_t = time.perf_counter()
        
        preds = graph_predict(np.array(batch_aa, dtype=np.int32), np.array(batch_feat, dtype=np.float32)).numpy()

        end_t = time.perf_counter()
        
        batch_time_ms = (end_t - start_t) * 1000.0
        avg_time_per_protein_ms = batch_time_ms / len(batch)

        for j, item in enumerate(batch):
            v_idx = np.where(item['loss_mask'] == 1.0)[0]
            predictions_pool.append({
                'pid': item['pid'],
                'aa_seq': item['aa_seq'],
                'y_pred': preds[j, :len(item['aa_seq']), 0][v_idx]
            })
            
            timing_records.append({
                'Protein_ID': item['pid'],
                'Inference_Time_ms': round(avg_time_per_protein_ms, 4)
            })

    predictions_txt_path = os.path.join(output_dir, "binding_predictions.txt")
    print(f"[Info] Exporting sequence predictions to: {predictions_txt_path}")
    with open(predictions_txt_path, 'w') as f:
        for prot in predictions_pool:
            f.write(f">{prot['pid']}\n")
            f.write(f"".join([REVERSE_AA_DICT.get(x, 'X') for x in prot['aa_seq']]) + "\n")
            f.write(f" ".join([f"{p:.4f}" for p in prot['y_pred']]) + "\n")

    # Export millisecond inference timing log
    timing_csv_path = os.path.join(output_dir, "inference_times.csv")
    print(f"[Info] Exporting performance profiling to: {timing_csv_path}")
    pd.DataFrame(timing_records).to_csv(timing_csv_path, index=False)
    
    print("\n[Success] Inference pipeline completed flawlessly.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DisoBindModel General Prediction Engine')
    parser.add_argument('--input_path', type=str, required=True, help='Path to input sequence table file (.txt)')
    parser.add_argument('--weights_path', type=str, required=True, help='Path to pre-trained model weights (.keras)')
    parser.add_argument('--output_dir', type=str, default='./output', help='Directory to save outputs')
    parser.add_argument('--device', type=str, choices=['cpu', 'gpu'], default='gpu', help='Target device for inference execution')
    args = parser.parse_args()
    
    # Configure hardware target
    if args.device.lower() == 'cpu':
        tf.config.set_visible_devices([], 'GPU')
        print("[Config] Hardware target set to CPU exclusively.")
    else:
        print("[Config] Hardware target set to GPU (if available).")
        
    run_predictor(args.input_path, args.weights_path, args.output_dir)