# iDBRpred: An accurate predictor of intrinsically disordered binding residues in protein sequences

iDBRpred is an advanced deep learning framework based on a dual-branch Transformer architecture (featuring RoPE and SwiGLU), designed for high-accuracy prediction of binding residues within intrinsically disordered regions (IDRs) of protein sequences. 

---

## Directory Structure

```text
iDBRpred/
├── config.yaml              # Global configuration
├── environment.yml          # The required dependencies
├── predict.py               # Main prediction script
├── src/                     # Model modules
│   ├── __init__.py
│   ├── features.py          
│   ├── model.py            
│   └── utils.py             
├── feature_extractors/      # Feature generation
│   ├── extractors.py        
│   ├── AUIPred_v3/          
│   ├── GENN+ASAquick/       
│   ├── MoRFchibi_v1/        
│   ├── ESpritz/            
│   └── ncbi-blast-2.17.0+/  
├── model_weights/           # Model weights
│   └── Best_B_res_Transformer.keras
├── examples/                # Example inputs and outputs
│   ├── example_input.fasta
│   └── expected_output.txt
├── LICENSE   
└── README.md                

```

---

## Quick Start & Usage

### 1. Installation & Environment

Clone the repository and ensure you have the required dependencies installed:

```bash
git clone https://github.com/yujing-xy/DisoBindModel.git
cd iDBRpred
conda env create -f environment.yml -n iDBRpred_env
conda activate iDBRpred_env

```

### 2. Initialization (Grant Permissions)

You must grant execution permissions before your first run the predictor:

```bash
chmod +x feature_extractors/GENN+ASAquick/bin/*
chmod +x feature_extractors/GENN+ASAquick/asaquick
chmod -R +x feature_extractors/MoRFchibi_v1/MCS1.03/mc*
chmod +x feature_extractors/ncbi-blast-2.17.0+/bin/*
chmod -R +x feature_extractors/ESpritz/

```

### 3. Configuration Setup

Before running predictions, open `config.yaml` and update the absolute paths for your local environment:

* **Databases:** Set swissprot_path and uniref90_path under the 'databases' section to your local SwissProt and UniRef90 database installations.

* **ESM-2 Weights:** Set `model_dir` under the `esm2` section to your local esm2_t33_650M_UR50D weights directory. The weights must include 'esm2_t33_650M_UR50D.pt' and 'esm2_t33_650M_UR50D-contact-regression.pt'. Leave it empty (`""`) if you want the pipeline to auto-download and cache the weights in the default directory: `feature_extractors/Alignment_free_conservation_ESM2/ESM2_weights/checkpoints/`.

### 4. Run Predictions

Execute the main predictor script to perform batch inference directly from a FASTA file containing single or multiple protein sequences:

```bash
python predictor.py \
    -i ./examples/example_input.fasta \
    -o ./output_results \
    -c config.yaml \
    -d gpu \
    -t 2

```

**Command-line Arguments:**

* `-i`, `--input`: Input FASTA file path (Required).
* `-o`, `--outdir`: Output directory path where results will be saved (Required).
* `-c`, `--config`: Configuration file path (Default: `config.yaml`).
* `-d`, `--device`: Target computing device (`cpu` or `gpu`).
* `-t`, `--threads`: Number of CPU threads to use.

---

## Expected Outputs

The script will automatically generate the specified output directory (`-o`) containing:

* **`binding_predictions.txt`**: The predicted probability of binding for each residue formatted under FASTA headers.

```text
>Protein_ID_001
MKVIFLALLVSTISSVFAAA...
0.0142 0.0215 0.0956 0.7845 ... 0.0031

```

* **`times.csv`**: A performance profiling file recording the execution time (in milliseconds) required to process each protein sequence.

---
