#!/bin/bash
#SBATCH --output=annotate_amen_var_%j.out
#SBATCH --error=annotate_amen_var_%j.err
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -n 12

set -euo pipefail

# sbatch slurm_annotate_vcf.sh <input.vcf>            single VCF
# sbatch slurm_annotate_vcf.sh -I <vcf_list.txt>      one VCF per line

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
batch_size=${BATCH_SIZE:-500}

if [ -n "${CONDA_ENV:-}" ]; then
    source "${CONDA_BASE:-$HOME/anaconda3}"/etc/profile.d/conda.sh
    conda activate "${CONDA_ENV}"
fi

export CUDA_VISIBLE_DEVICES=${SLURM_JOB_GPUS:-${CUDA_VISIBLE_DEVICES:-0}}
export TF_FORCE_GPU_ALLOW_GROWTH=true

time python "${script_dir}"/annotate_vcf.py "$@" -b "${batch_size}"
