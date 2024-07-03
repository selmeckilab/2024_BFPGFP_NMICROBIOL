#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=6gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --time=9:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

bam_dir=bam
logs_dir=logs
temp_dir=/scratch.global/scot0854/bfpgfp  #scratch directory
file_name="$(date +%F)_Calbicans_bfpgfp_qc"

# Use local modules
module use /home/selmecki/shared/software/modulefiles.local

# Load modules
module load fastqc/0.11.9
module load qualimap/20231012
module load multiqc/20230928

# raw fastq file qc
find "$temp_dir" -name "*trimmed_*P.fq" -exec fastqc  -t 8 "$temp_dir" {} \;

# bam qc
unset DISPLAY
find "$bam_dir" -name "*.bam" \
-exec qualimap bamqc -bam {} -outdir $temp_dir/{} --java-mem-size=4G  \;

# multiqc
multiqc "${temp_dir}" "${logs_dir}" "${bam_dir}" -o logs -n "$file_name"
