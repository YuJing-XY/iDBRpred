import os
import sys
import urllib.request

def parse_fasta(file_path):
    sequences = {}
    if not os.path.exists(file_path):
        sys.exit(f"Error: Input file not found: {file_path}")
        
    with open(file_path, 'r') as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    sequences[header] = "".join(seq)
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if header:
            sequences[header] = "".join(seq)
    return sequences

def check_and_download_esm2(config_esm2, base_dir):
    model_name = config_esm2.get('model_name', 'esm2_t33_650M_UR50D')
    custom_dir = config_esm2.get('model_dir', '')
    
    files_to_download = {
        f"{model_name}.pt": f"https://dl.fbaipublicfiles.com/fair-esm/models/{model_name}.pt",
        f"{model_name}-contact-regression.pt": f"https://dl.fbaipublicfiles.com/fair-esm/regression/{model_name}-contact-regression.pt"
    }
    
    use_custom = False
    if custom_dir and os.path.exists(custom_dir):
        all_files_exist = all(
            os.path.exists(os.path.join(custom_dir, fname)) 
            for fname in files_to_download.keys()
        )
        if all_files_exist:
            use_custom = True

    if use_custom:
        print(f"Found existing ESM-2 weights in custom directory: {custom_dir}")
        return
        
    if custom_dir:
        print(f"Weights missing or incomplete in custom directory. Falling back to default directory.")

    default_hub_dir = os.path.join(base_dir, 'feature_extractors', 'Alignment_free_conservation_ESM2', 'ESM2_weights')
    default_checkpoints_dir = os.path.join(default_hub_dir, 'checkpoints')
    
    os.makedirs(default_checkpoints_dir, exist_ok=True)
    
    for file_name, url in files_to_download.items():
        file_path = os.path.join(default_checkpoints_dir, file_name)
        if not os.path.exists(file_path):
            print(f"Downloading {file_name} to {default_checkpoints_dir}. This may take a while...")
            try:
                urllib.request.urlretrieve(url, file_path)
            except Exception as e:
                sys.exit(f"Error downloading {file_name}: {str(e)}")