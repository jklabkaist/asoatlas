#!/bin/bash
# Link the reference files required by the vendored jkspliceai / jklib packages.
#
#   setup_data.sh <source_dir> [--copy]
#
# <source_dir> holds the files below; they are not distributed with this repository:
#   jkspliceai/models/spliceai{1..5}.h5                  SpliceAI models (Illumina, GPL-3.0)
#   jkspliceai/annotations/MANE.GRCh38.v1.0.summary.txt  ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.0/
#   jkspliceai/annotations/gencode.v41.annotation.gtf.refFlat.MANE.v1.arranged.txt
#                                                        refFlat of the MANE transcripts of GENCODE v41
#   jkspliceai/annotations/labranchor.gencode_v40.top.sorted.bed.gz(.tbi)
#                                                        LaBranchoR branchpoint predictions (GENCODE v40)
#   jklib/datafiles/hg38.2bit                            hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/

set -euo pipefail

src=${1:?usage: $0 <source_dir> [--copy]}
mode=${2:-}

src=$(cd "${src}" && pwd)
vendor=$(cd "$(dirname "${BASH_SOURCE[0]}")/vendor" && pwd)

# accept either <src>/<pkg>/... or the source checkout layout <src>/<pkg>/<pkg>/...
sai_src=${src}/jkspliceai
[ -d "${sai_src}/jkspliceai" ] && sai_src=${sai_src}/jkspliceai
jklib_src=${src}/jklib
[ -d "${jklib_src}/jklib" ] && jklib_src=${jklib_src}/jklib

link() {
    local from=$1 to=$2
    if [ ! -e "${from}" ]; then
        echo "missing: ${from}" >&2
        return
    fi
    mkdir -p "$(dirname "${to}")"
    if [ "${mode}" == "--copy" ]; then
        cp -f "${from}" "${to}"
    else
        ln -sf "${from}" "${to}"
    fi
    echo "  ${to} -> ${from}"
}

for model in 1 2 3 4 5; do
    link "${sai_src}/models/spliceai${model}.h5" "${vendor}/jkspliceai/models/spliceai${model}.h5"
done
link "${sai_src}/models/LICENSE" "${vendor}/jkspliceai/models/LICENSE"

for annot in MANE.GRCh38.v1.0.summary.txt \
             gencode.v41.annotation.gtf.refFlat.MANE.v1.arranged.txt \
             labranchor.gencode_v40.top.sorted.bed.gz \
             labranchor.gencode_v40.top.sorted.bed.gz.tbi; do
    link "${sai_src}/annotations/${annot}" "${vendor}/jkspliceai/annotations/${annot}"
done

link "${jklib_src}/datafiles/hg38.2bit" "${vendor}/jklib/datafiles/hg38.2bit"

echo "done."
