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

list=manual_vcfs_to_merge.txt
merged=final_tmp.vcf
gene_annotation=Calbicans_SC5314_A21_sorted_genes.bed.gz
final_vcf=BFPGFP_annotated_filtered.vcf

# load modules, use local for newer bcftools
module use /home/selmecki/shared/software/modulefiles.local

module load bcftools/1.17
module load htslib/1.9

# merge individual VCFs
bcftools merge -l "${list}" -o "${merged}"
bgzip "${merged}"
tabix "${merged}".gz

bcftools annotate -a "${gene_annotation}" \
  -c CHROM,FROM,TO,GENE \
  -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') \
  "${merged}".gz > "${final_vcf}"

bgzip "${final_vcf}"
tabix "${final_vcf}".gz

samples=$(bcftools query -l "${final_vcf}".gz | sort)

for s in $samples; do
    bcftools view -s $s "${final_vcf}".gz \
    | bcftools view -i "FORMAT/AD[0:1] >4" \
    | bcftools view -i "FORMAT/AF[0:0] >=0.2" \
    -o "${s}"_filtered.vcf
done
