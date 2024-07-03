#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 50
#SBATCH -p msilarge,msismall
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=4-120

set -ue
set -o pipefail

sample_file=bfpgfp_samples.txt
vcf=BFPGFP_annotated_filtered.vcf.gz
line=${SLURM_ARRAY_TASK_ID}

# load modules, use local for newer bcftools
module use /home/selmecki/shared/software/modulefiles.local

module load bcftools/1.17
module load gatk/4.1.2

# remove int files
function finish {
  arr=("${sample}".table "${sample}".ann.table "${sample}".first_ann.table)
  for f in "${arr[@]}"; do
    if [ -e "$f" ]; then
        rm "$f"
    fi
  done
}

trap finish EXIT

sample=$(awk -v val=$line 'NR == val { print $0}' $sample_file)

bcftools view -s "${sample}"  "${vcf}" \
    | bcftools view -i 'FORMAT/AD[0:1] > 4' \
    | bcftools view -i 'FORMAT/AF[0:0] >= 0.07' > "${sample}"_final_filtered.vcf

[[ $(grep -vc "^#"  "$sample"_final_filtered.vcf ) -ne 0 ]] && {

# get info in tabular form
  gatk VariantsToTable \
    -V "${sample}"_final_filtered.vcf \
    -F CHROM -F POS -F GENE -F REF -F ALT -GF GT -GF AF \
    -O "${sample}".table

# grab annotation
  bcftools query -f '%INFO/ANN\n' \
    "${sample}"_final_filtered.vcf > "${sample}".ann.table

# keep only the first
  awk 'BEGIN{FS="|"; OFS="\t"} {print $2, $3, $4, $10, $11}' \
    "${sample}".ann.table > "${sample}".first_ann.table

# fix the header
  sed -i '1i Annotation\tImpact\tORF\tCoding_Change\tAA_change' \
    "${sample}".first_ann.table

# do row numbers still match?
  [[ $(wc -l < "$sample.table") -ne $(wc -l < "$sample.first_ann.table") ]] && { >&2 echo "table lengths don't match"; exit 1; }

# final csv output
  paste "${sample}".table "${sample}".first_ann.table \
    | awk 'BEGIN{OFS=","} {print $1, $2, $3, $10, $4, $5, $6, $7, $8, $9, $11,$12}' > "${sample}".csv
  sed -i "s/^/$sample,/" "${sample}".csv
}
