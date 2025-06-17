# Germline variant calling pipeline using nextflow

Understanding germline variants is crucial for unraveling the complexities of human genetics and their implications for health and disease. Germline variants are heritable genetic changes that can influence an individual’s susceptibility to various conditions, including cancer, cardiovascular diseases, and other hereditary disorders. By accurately identifying these variants, researchers and clinicians can make informed decisions about diagnosis, treatment, and prevention strategies.

This repository presents a streamlined bioinformatics workflow designed specifically for germline variant calling. Leveraging aligned genomic data for each sample, our pipeline utilizes Nextflow alongside the Genome Analysis Toolkit (GATK) and BCFtools to ensure comprehensive variant identification.

## Key Features

**Nextflow Workflow Management:** Nextflow is a powerful workflow management system that enables reproducible and scalable data analysis. It allows users to define complex workflows in a simple, concise manner while ensuring compatibility across different computing environments, from local machines to cloud infrastructures.

**GATK for Variant Calling:** The Genome Analysis Toolkit (GATK) is a widely used toolkit for variant discovery in high-throughput sequencing data. It provides state-of-the-art algorithms for variant calling and is particularly known for its best practices in the analysis of germline variants, ensuring high sensitivity and specificity.

**BCFtools for Variant Manipulation:** BCFtools is a set of utilities that facilitate the manipulation of variant call format (VCF) and binary call format (BCF) files. It allows users to filter, annotate, and manipulate variant data efficiently, making it an essential tool for post-variant calling analysis.

### BCFtools Code 

``
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
    samples = Channel.fromPath('data/*.bam')
    roi_file = file('data/ROI.bed')
    ref = file('data/Homo_sapiens.GRCh38.dna.primary_assembly.fa')

    variant_calling(samples, roi_file, ref)
}
``
