#!/bin/bash
# ASO amenability of the 6 VEP-annotated variants in example.vcf
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

python ../annotate_vcf.py -i example.vcf -o example.amenability.vcf -b 100

if diff -q example.amenability.vcf expected/example.amenability.vcf > /dev/null; then
    echo "example.amenability.vcf matches expected/example.amenability.vcf"
else
    echo "example.amenability.vcf differs from expected/example.amenability.vcf"
    diff <(grep -v '^#' expected/example.amenability.vcf) <(grep -v '^#' example.amenability.vcf) | head
fi
