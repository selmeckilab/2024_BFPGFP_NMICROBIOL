#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=
#SBATCH -t 50
#SBATCH -p msilarge,msismall
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

list=ams5178_to_merge.txt
merged=ams5178_tmp.vcf
regions=SC5314_A21_region_file.bed
vcf_file=bfpgfp_ams5178.vcf
snpeff=/home/selmecki/shared/software/snpEff/snpEff.jar
snpeff_config=/home/selmecki/shared/software/snpEff/snpEff.config
snpeff_db=SC5314_s02m09r08
out_file=BFPGFP_AMS5178_annotated.vcf

# load modules, use local for newer bcftools
module use /home/selmecki/shared/software/modulefiles.local

module load bcftools/1.17
module load htslib/1.9

# remove int files
function finish {
  rm "${vcf_file}"
  rm "${merged}"*
}

trap finish EXIT

# merge individual VCFs
bcftools merge -l "${list}" -o "${merged}"
bgzip "${merged}"
tabix "${merged}".gz

bcftools view "${merged}".gz -R "${regions}" \
  | bcftools view -i "FORMAT/AD[*:1] >4" \
  | bcftools view -i "FORMAT/F1R2[*] > 0 & FORMAT/F2R1[*] > 0" \
  | bcftools view -i "F_PASS(GT='alt') < 0.5" \
  | bcftools view -i "FORMAT/AF[*] >=0.2" \
  -o "${vcf_file}"

# Annotate using snpeff with manual database (CTG)
java -Xmx4g -jar "${snpeff}" \
    -ud 100 \
    -c "${snpeff_config}" \
    "${snpeff_db}" \
    "${vcf_file}" \
    > "${out_file}"

bgzip "${out_file}"
tabix "${out_file}".gz
