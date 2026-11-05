#!/bin/bash
#SBATCH --job-name=fes-gpu-test
#SBATCH --partition=mig
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -e

module purge
module load anaconda3/2026.7
conda activate fes-esm

python - <<'PY'
import torch

print("PyTorch:", torch.__version__)
print("CUDA available:", torch.cuda.is_available())

if not torch.cuda.is_available():
    raise RuntimeError("CUDA is unavailable")

print("GPU:", torch.cuda.get_device_name(0))

x = torch.randn(1000, 1000, device="cuda")
print("Test output:", (x @ x).shape)
PY
