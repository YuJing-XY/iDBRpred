# DisoBindModel: A Dual-Branch Transformer for Protein Disordered Binding Residue Prediction

DisoBindModel is an advanced deep learning framework based on a dual-branch Transformer architecture, designed for high-accuracy prediction of binding residues within intrinsically disordered regions (IDRs) of protein sequences. 
---

## Directory Structure

```text
DisoBindModel/
├── configs.py               # Configuration
├── feature_extractors/      # Feature generation
│   ├── aiupred/                
│   ├── asaquick/                
│   ├── morfchibi/               
│   ├── esm2_conservation/      
│   └── raap/                    
├── modules.py                # Model framework
├── model_weights/
│   └── Best_B_res_Transformer.keras # Model weights
├── predict.py                # Main predictor script
├── LICENSE   
├── environment.yml           # The required dependencies      
└── README.md               

```

---

## Quick Start & Usage

### 1. Installation & Environment

Clone the repository and ensure you have the required dependencies installed (TensorFlow 2.x, NumPy, Pandas):

```bash
git clone [https://github.com/yujing-xy/DisoBindModel.git](https://github.com/yujing-xy/DisoBindModel.git)
cd DisoBindModel
conda env create -f environment.yml -n DisoBindModel

```

### 2. Run Predictions

Execute the standalone prediction script to run batch inference over your plain-text protein feature sheets. You can specify whether to run the model on a CPU or GPU using the `--device` argument:

```bash
python predict.py \
    --input_path ./data/sample_input.txt \
    --weights_path ./weights/Best_B_res_Transformer.keras \
    --output_dir ./output \
    --device gpu

```

*(Note: Use `--device cpu` if you want to force computation on the CPU, which is useful for environments without compatible CUDA drivers).*

### 3. Expected Outputs

After a successful run, the script will automatically save two result files in your specified output directory:

* binding_predictions.txt: A plain-text prediction file structured as follows:
```text
>Protein_ID_001
MKVIFLALLVSTISSVFAAA...
0.0142 0.0215 0.0956 0.7845 ... 0.0031

```


* times.csv: A performance profiling file recording the inference time (in milliseconds) required to process each individual protein sequence.
