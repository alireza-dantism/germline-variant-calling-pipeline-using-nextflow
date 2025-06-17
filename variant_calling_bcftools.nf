#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process variant_calling {
    input:
    path bam_file
    path roi_bed
    path reference

    output:
    path "output.bcf"
    path "output.vcf"
    path "output.tsv"

    script:
    """
    # Check if the BAM file exists
    if [ ! -f "${bam_file}" ]; then
        echo "Error: BAM file not found: ${bam_file}"
        exit 1
    fi

    # Create index for the BAM file if it doesn't exist
    samtools index ${bam_file}

    # Generate a BCF file from the pileup with increased max depth
    bcftools mpileup -f ${reference} -R ${roi_bed} ${bam_file} -Ou -o temp.bcf

    # Call variants and save to a new output BCF file
    bcftools call -mv -Ob -o output.bcf temp.bcf

    # Convert BCF to VCF
    bcftools view output.bcf > output.vcf

    # Convert VCF to TSV
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' output.vcf > output.tsv

    # Clean up temporary files
    rm temp.bcf
    """
}

workflow {
    // Define your samples (BAM files, ROI file and refrencec_genome file
    samples = Channel.fromPath('/Users/msg/Documents/bioinformatics-project/germline_variant_calling/data/*.bam')
    roi_file = file('/Users/msg/Documents/bioinformatics-project/germline_variant_calling/data/ROI.bed')
    ref = file('/Users/msg/Documents/bioinformatics-project/germline_variant_calling/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa')

    variant_calling(samples, roi_file, ref)
}
