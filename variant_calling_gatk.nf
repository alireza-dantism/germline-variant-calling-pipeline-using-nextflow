#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process variant_calling_gatk {
    input:
    path bam_file

    output:
    path "${bam_file.baseName}.vcf" // Output VCF file

    script:
    """
	/Users/msg/Documents/bioinformatics-project/gatk/gatk  --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' HaplotypeCaller \
    	-I ${bam_file} \
    	-O ${bam_file.baseName}.vcf \
    	-R /Users/msg/Documents/bioinformatics-project/germline_variant_calling/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    	-L /Users/msg/Documents/bioinformatics-project/germline_variant_calling/data/ROI.bed
    """
}

// Define the process to convert VCF files to table format
process convert_vcf_to_table {
    // Define input and output
    input:
    path vcf_file

    output:
    path "${vcf_file.baseName}.tsv" // Output TSV file

    script:
    """
    # Convert VCF to TSV using bcftools or any other tool
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' $vcf_file > ${vcf_file.baseName}.tsv
    """
}


workflow {
    // Define your samples (BAM files)
    samples = Channel.fromPath('data/*.bam')

    // Variant calling with GATK
    vcf_files = variant_calling_gatk(samples)

    //convert_vcf_to_table(samples);
    table_output = vcf_files | convert_vcf_to_table
}
