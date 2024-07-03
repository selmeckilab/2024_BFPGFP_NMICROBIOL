#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=3gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=
#SBATCH -t 40
#SBATCH -p msismall,msilarge
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=1-120

set -ue
set -o pipefail

ref_fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta
sample_file=vcfs_to_pass.txt
line=${SLURM_ARRAY_TASK_ID}

# Load modules
module use /home/selmecki/shared/software/modulefiles.local

module load bcftools/1.17
module load htslib/1.9

#Get strain ID froms sample file line equal to array task ID
sample=$(awk -v val=$line 'NR == val { print $0}' $sample_file)
strain=$(basename $sample | cut -d "_" -f 1)

bcftools view -s "${strain}" -f "PASS" "${sample}" \
    | bcftools norm -m -any -f "${ref_fasta}" \
    > pass_vcfs/"${strain}"_mutect2_filtered_pass_norm.vcf

bgzip pass_vcfs/"${strain}"_mutect2_filtered_pass_norm.vcf
tabix pass_vcfs/"${strain}"_mutect2_filtered_pass_norm.vcf.gz
