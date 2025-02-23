## Single-cell detection of copy number changes reveals dynamic mechanisms of adaptation to antifungals in Candida albicans

<br>

Parallel evolution experiments of fluorescently-tagged *C. albicans* were performed in the presence and absence of antifungal drugs.

This repository summarizes variant detection from WGS data and visualization of population dynamics from flow cytometry data.

### Illumina read trimming, alignment and QC
Adapter trimming, contaminant checking and quality trimming were performed for all fastq files using BBDuk (part of BBTools, v38.94). Reads were aligned to *C. albicans* SC5314 reference genome, version A21-s02-m09-r08 using BWA-MEM (v0.7.17). Aligned reads were sorted by coordinate, duplicates were marked and bam files indexed by samtools (v1.10) (see script "array_bbduk_align_sort.sh").

Raw read and alignment quality were assessed with FastQC (v0.11.9), qualimap (v.2.2.2-dev) and samtools (v1.10) flagstat, and summarized with MultiQC (v1.16) (see script "basic_qc.sh").

### Variant calling

Variant calling and preliminary filtering was performed with Mutect2 and FilterMutectCalls (GATK v4.1.2), both with default parameters.

Variant calling was run in 3 batches corresponding to progenitor strains

1.  SC5314 (AMS2401) progenitor treated as "normal"
2.  Chr3 BFP/GFP Progenitor (AMS5178) progenitor treated as "normal"
3.  Chr5 BFP/GFP progenitor (AMS5192) progenitor treated as "normal"

See scripts "array_mutect2_sc5314.sh", "array_mutect2_ams5178.sh" and "array_mutect2_ams5192.sh" for details.

### Variant filtering
Additional VCF filtering was performed with bcftools (v1.17). Individual mutect2 vcf files were then subset to remove the progenitor strain, and filtered on quality = "PASS" (see "array_filter_individual_vcf.sh"). Groups of individual VCFs were merged into 3 different files based on their progenitor (SC5314, AMS5178 and AMS5192). Merged files were subset to exclude known repetitive regions as annotated in the SC5314 A21-s02-m09-r08 GFF (rRNA, repeat_region, retrotransposon) and telomere-proximal regions, defined here as extending from each chromosome end to the first non-repetitive genome feature, and then filtered on the following parameters:

-   At least 5 supporting reads for alternate alleles
-   Supporting reads in both directions
-   Alternate alleles found in less than half of analyzed progeny (half or more suggesting the variant originated in the progenitor)
-   alternate allele frequency of at least 20% in at least 1 isolate (diploid single colonies) or at least 5% (polyploid samples and populations)

Variants were annotated with snpEff (v5.0e, database built from SC5314 version A21-s02-m09-r08, with alternate yeast nuclear codon table specified) (see scripts "merge_evolved_sc5314.sh", "merge_evolved_ams5178.sh", "merge_lowfreq_ams5178.sh", "merge_evolved_ams5192.sh" and "merge_lowfreq_ams5192.sh").

### Manual curation of variants
Variants were visually inspected in IGV (v2.16.1) for each of the three groups of evolved isolates. Working copies of each VCF file were edited to exclude false positives identified during manual review. The manually-curated VCFs were merged to a single, final variant file including all 120 strains, and the merged VCF was annotated with gene names (see "merge_filtered_vcfs.sh"). A CSV file of manually curated, annotated variants was produced from the final VCF (see "variant_subset_table.sh").

### Visualization of flow cytometry data
Flow cytometry data of Chr3 BFP/GFP and Chr5 BFP/GFP evolved lineages were collected at passages 0, 1, 5, 10, 12 and 15. Fractions of five different types of fluorescent cells were calculated and plotted in R (see "Stacked_area_script.R").
