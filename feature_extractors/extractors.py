import os
import sys
import subprocess
import tempfile
import shutil
import numpy as np
import torch

class AIUPredExtractor:
    def __init__(self, force_cpu=False, gpu_num=0):
        self.aiupred_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'feature_extractors', 'AIUPred_v3')
        if self.aiupred_dir not in sys.path:
            sys.path.insert(0, self.aiupred_dir)
        from aiupred import AIUPred
        self.predictor = AIUPred(force_cpu=force_cpu, gpu_num=gpu_num)

    def extract(self, sequence: str) -> np.ndarray:
        return self.predictor.predict_disorder(sequence)

class ASAquickExtractor:
    def __init__(self, base_dir: str):
        self.asaquick_dir = os.path.join(base_dir, 'feature_extractors', 'GENN+ASAquick')
        self.bin_dir = os.path.join(self.asaquick_dir, 'bin')
        self.executable = os.path.join(self.asaquick_dir, 'asaquick')

    def extract(self, sequence: str) -> np.ndarray:
        seq_len = len(sequence)
        fd, temp_fasta_path = tempfile.mkstemp(suffix='.fasta', text=True)
        with os.fdopen(fd, 'w') as f:
            f.write(f">temp_seq\n{sequence}\n")

        fasta_basename = os.path.basename(temp_fasta_path)
        out_dir_path = os.path.join(self.asaquick_dir, f"asaq.{fasta_basename}")

        env = os.environ.copy()
        env["PATH"] = f"{self.bin_dir}:{env.get('PATH', '')}"

        try:
            subprocess.run(
                [self.executable, temp_fasta_path],
                cwd=self.asaquick_dir,
                env=env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True
            )

            rasa_values = []
            rasa_pred_file = os.path.join(out_dir_path, 'rasaq.pred')
            
            if os.path.exists(rasa_pred_file):
                with open(rasa_pred_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith('#'): continue
                        parts = line.split()
                        if len(parts) >= 3:
                            rasa_values.append(float(parts[2]))
            
            if len(rasa_values) != seq_len:
                return np.zeros(seq_len, dtype=np.float32)

            return np.array(rasa_values, dtype=np.float32)

        except Exception:
            return np.zeros(seq_len, dtype=np.float32)

        finally:
            if os.path.exists(temp_fasta_path):
                os.remove(temp_fasta_path)
            if os.path.exists(out_dir_path):
                shutil.rmtree(out_dir_path)

class ESM2Extractor:
    def __init__(self, base_dir: str, model_name: str = 'esm2_t33_650M_UR50D', model_dir: str = '', device: str = 'cuda', threads: int = 8):
        self.kibby_dir = os.path.join(base_dir, 'feature_extractors', 'Alignment_free_conservation_ESM2', 'kibby')
        if self.kibby_dir not in sys.path:
            sys.path.insert(0, self.kibby_dir)
            
        default_target_dir = os.path.join(base_dir, 'feature_extractors', 'Alignment_free_conservation_ESM2', 'ESM2_weights')
        target_dir = default_target_dir 

        if model_dir and os.path.exists(os.path.join(model_dir, f"{model_name}.pt")):
            if model_dir.rstrip('/').endswith("checkpoints"):
                target_dir = os.path.dirname(model_dir.rstrip('/'))
            else:
                print(f"  -> [Warning] ESM2 custom weights found in {model_dir}, but the folder is not named 'checkpoints'.")
                print("  -> [Warning] PyTorch Hub requires this specific folder name. Falling back to default directory to enforce consistency.")
  

        torch.hub.set_dir(target_dir)
            
        from my_library import ESM_Model, RegressionModel
        self.estimate_full_length = __import__('my_library').estimate_full_length
        
        self.device = device if torch.cuda.is_available() and device == 'cuda' else 'cpu'
        self.threads = threads
        
        self.esm = ESM_Model()
        self.esm.load(model_name)
        
        npz_file = os.path.join(self.kibby_dir, 'linear_models', f'{model_name}.npz')
        self.regression_model = RegressionModel.from_npz(npz_file)

    def extract(self, sequence: str) -> np.ndarray:
        try:
            def estimate_chunk(x):
                embeddings = self.esm.encode(x, device=self.device, threads=self.threads)
                return self.regression_model.predict(embeddings[1:-1])
                
            conservation = self.estimate_full_length(
                sequence, 
                estimate_chunk, 
                chunk_size=1022, 
                min_overlap=300
            )
            return np.array(conservation, dtype=np.float32)
        except Exception as e:
            print(f"Warning: ESM2 extraction failed. Padding with zeros. Error: {e}")
            return np.zeros(len(sequence), dtype=np.float32)


class MoRFchibiExtractor:
    def __init__(self, base_dir: str, mode: str, databases_config: dict, threads: int = 2):
        self.mcs_dir = os.path.join(base_dir, 'feature_extractors', 'MoRFchibi_v1', 'MCS1.03')
        self.blast_bin = os.path.join(base_dir, 'feature_extractors', 'ncbi-blast-2.17.0+', 'bin')
        
        self.espritz_dir = os.path.join(base_dir, 'feature_extractors', 'ESpritz')
        if not self.espritz_dir.endswith('/'):
            self.espritz_dir += '/'
            
        self.mode = mode

        mode_map = {
            'web': 'mcw',
            'light': 'mcl',
            'normal': 'mc'
        }
        self.executable = os.path.join(self.mcs_dir, mode_map.get(mode, 'mcw'))
        
        self.input_file = os.path.join(self.mcs_dir, 'input.fasta')
        self.output_file = os.path.join(self.mcs_dir, 'output.txt')
        
        self._generate_properties_file(databases_config, threads)

    def _generate_properties_file(self, db_config, threads):
        prop_file_path = os.path.join(self.mcs_dir, 'MoRFchibi.properties')
        
        swissprot_path = db_config.get('swissprot_path', '')
        uniref90_path = db_config.get('uniref90_path', '')
        
        properties_content = f"""#Debug mode
debug\t0
#ESpritz path
ESpritz\t{self.espritz_dir}
#BLAST path 
PSIBLAST\t{self.blast_bin}
#input fasta file
Input\tinput.fasta
#output file
Output\toutput.txt   
# UniProt Swiss-Prot database
SwissProt\t{swissprot_path}
# UniProt UniRef90 database
UniRef90\t{uniref90_path}
# number of threads
threads\t{threads}
"""
        with open(prop_file_path, 'w') as f:
            f.write(properties_content)

    def extract(self, sequence: str) -> np.ndarray:
        with open(self.input_file, 'w') as f:
            f.write(f">temp_seq\n{sequence}\n")

        if os.path.exists(self.output_file):
            try:
                os.remove(self.output_file)
            except OSError:
                pass

        env = os.environ.copy()
        if os.path.exists(self.blast_bin):
            env["PATH"] = f"{self.blast_bin}:{env.get('PATH', '')}"

        try:
            subprocess.run(
                [self.executable],
                cwd=self.mcs_dir,
                env=env,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

            scores = []
            
            # Normal (MC), Light (MCL), Web (MCW) 
            target_idx = 2 

            if os.path.exists(self.output_file):
                with open(self.output_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith('#') or line.startswith('>'):
                            continue
                        
                        parts = line.split()
                        if len(parts) > target_idx:
                            scores.append(float(parts[target_idx]))

            if len(scores) != len(sequence):
                return np.zeros(len(sequence), dtype=np.float32)

            return np.array(scores, dtype=np.float32)

        except Exception as e:
            print(f"Warning: MoRFchibi extraction failed. Padding with zeros. Error: {e}")
            return np.zeros(len(sequence), dtype=np.float32)

        finally:
            if os.path.exists(self.input_file):
                os.remove(self.input_file)
            if os.path.exists(self.output_file):
                try:
                    os.remove(self.output_file)
                except OSError:
                    pass