#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=3gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=
#SBATCH --time=8:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=1-65

set -ue
set -o pipefail

ref_fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta
sample_file=ams5178_bam_files.txt
normal_bam=bam/AMS5178_trimmed_bwa_sorted_markdup.bam
normal_name=AMS5178
line=${SLURM_ARRAY_TASK_ID}

# Load modules
module load htslib/1.9
module load gatk/4.1.2

#Get strain ID froms sample file line equal to array task ID
file=$(awk -v val=$line 'NR == val { print $1}' $sample_file)
strain=$(basename "$file" | cut -d "_" -f 1)

function finish {
  rm -f "${strain}"*unfiltered.vcf*
  rm -f "${strain}"*.idx
}
trap finish EXIT

gatk Mutect2 \
    -R "${ref_fasta}" \
    -I "${file}" \
    -I "${normal_bam}" \
    -normal "${normal_name}" \
    -O "${strain}"_unfiltered.vcf

gatk FilterMutectCalls \
    -R "${ref_fasta}" \
    -V "${strain}"_unfiltered.vcf \
    -O "${strain}"_mutect2_filtered.vcf

bgzip "${strain}"_mutect2_filtered.vcf
tabix "${strain}"_mutect2_filtered.vcf.gz

mv "${strain}"_mutect2_filtered.vcf.gz* raw_vcfs/
