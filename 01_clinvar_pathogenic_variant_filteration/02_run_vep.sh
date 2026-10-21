#!/bin/bash
#SBATCH --output=VEP_%j.out
#SBATCH --error=VEP_%j.err
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=4G

set -euo pipefail

# Required environment:
#   VEP_SIF        VEP singularity image
#   VEP_CACHE_DIR  VEP offline cache directory
#   VEP_PLUGIN_DIR VEP plugin directory (contains UTRannotator/ and the EVE file)
#   VEP_FASTA_DIR  directory holding the bgzipped GRCh38 primary assembly fasta
#   REVEL_DIR      directory holding the tabixed REVEL table
# Optional: VEP_FASTA, REVEL_TSV, EVE_VCF, FORK, CONDA_BASE, CONDA_ENV

input=${1:?usage: $0 <input.vcf>}

VEP_SIF=${VEP_SIF:?}
VEP_CACHE_DIR=${VEP_CACHE_DIR:?}
VEP_PLUGIN_DIR=${VEP_PLUGIN_DIR:?}
VEP_FASTA_DIR=${VEP_FASTA_DIR:?}
REVEL_DIR=${REVEL_DIR:?}
VEP_FASTA=${VEP_FASTA:-Homo_sapiens.GRCh38.dna.primary_assembly.fa.bgz}
REVEL_TSV=${REVEL_TSV:-new_tabbed_revel_grch38.tsv.gz}
EVE_VCF=${EVE_VCF:-eve_merged.vcf.bgz}
FORK=${FORK:-8}

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
data_dir=$(cd "$(dirname "${input}")" && pwd)
input_base=$(basename "${input}")
output_base=${input_base//.vcf/.VEP.noskip.vcf}
output=${data_dir}/${output_base}

source /etc/profile.d/modules.sh 2>/dev/null && module load singularity || true

time singularity exec \
  -B "${data_dir}":/opt/vep/.vep \
  -B "${VEP_CACHE_DIR}":/opt/vep/.vep_cache \
  -B "${VEP_PLUGIN_DIR}":/opt/vep/.vep_plugins \
  -B "${REVEL_DIR}":/opt/vep/REVEL \
  -B "${VEP_FASTA_DIR}":/opt/vep/.vep_fasta \
  "${VEP_SIF}" env PERL5LIB="${PERL5LIB:-}:/opt/vep/.vep_plugins/UTRannotator" \
  /opt/vep/src/ensembl-vep/vep \
  -i /opt/vep/.vep/"${input_base}" \
  -o /opt/vep/.vep/"${output_base}" \
  --offline --fork "${FORK}" \
  --dir_cache /opt/vep/.vep_cache \
  --dir_plugins /opt/vep/.vep_plugins \
  --assembly GRCh38 \
  --fasta /opt/vep/.vep_fasta/"${VEP_FASTA}" \
  --vcf --use_given_ref --dont_skip --numbers --hgvs --symbol --transcript_version --canonical --mane --check_existing --keep_csq \
  --plugin REVEL,/opt/vep/REVEL/"${REVEL_TSV}" \
  --plugin NMD

if [ -n "${CONDA_ENV:-}" ]; then
    source "${CONDA_BASE:-$HOME/anaconda3}"/etc/profile.d/conda.sh
    conda activate "${CONDA_ENV}"
fi

time python "${script_dir}"/03_curate_vep.py -i "${output}"
