# Germline variant calling pipeline using nextflow

Understanding germline variants is crucial for unraveling the complexities of human genetics and their implications for health and disease. Germline variants are heritable genetic changes that can influence an individual’s susceptibility to various conditions, including cancer, cardiovascular diseases, and other hereditary disorders. By accurately identifying these variants, researchers and clinicians can make informed decisions about diagnosis, treatment, and prevention strategies.

This repository presents a streamlined bioinformatics workflow designed specifically for germline variant calling. Leveraging aligned genomic data for each sample, our pipeline utilizes Nextflow alongside **the Genome Analysis Toolkit (GATK)** and **BCFtools** to ensure comprehensive variant identification.

## Key Features

**Nextflow Workflow Management:** Nextflow is a powerful workflow management system that enables reproducible and scalable data analysis. It allows users to define complex workflows in a simple, concise manner while ensuring compatibility across different computing environments, from local machines to cloud infrastructures.

**GATK for Variant Calling:** The Genome Analysis Toolkit (GATK) is a widely used toolkit for variant discovery in high-throughput sequencing data. It provides state-of-the-art algorithms for variant calling and is particularly known for its best practices in the analysis of germline variants, ensuring high sensitivity and specificity.

**BCFtools for Variant Manipulation:** BCFtools is a set of utilities that facilitate the manipulation of variant call format (VCF) and binary call format (BCF) files. It allows users to filter, annotate, and manipulate variant data efficiently, making it an essential tool for post-variant calling analysis.

## BCFtools Code 

```

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

```

## Running Code

```
nextflow run variant_calling_bcftools.nf
```

![Running code](/command-bcftools.png)


## BCFtools Result

### Sample 1 - Output

|   1   |   67395837   |   .  |   C     |   A     |   225.417  |   .  |
|-------|--------------|------|---------|---------|------------|------|
|   1   |   179551371  |   .  |   G     |   A     |   225.417  |   .  |
|   10  |   98459557   |   .  |   G     |   A     |   225.417  |   .  |
|   11  |   16111862   |   .  |   G     |   A     |   14.3183  |   .  |
|   11  |   16111867   |   .  |   A     |   G     |   124.156  |   .  |
|   13  |   33129519   |   .  |   T     |   C     |   160.894  |   .  |
|   13  |   38859469   |   .  |   A     |   G     |   222.229  |   .  |
|   14  |   67575857   |   .  |   A     |   G     |   220.803  |   .  |
|   17  |   59886176   |   .  |   A     |   G     |   225.417  |   .  |
|   18  |   23833905   |   .  |   T     |   C     |   225.417  |   .  |
|   18  |   62360008   |   .  |   C     |   T     |   4.97023  |   .  |
|   18  |   63987229   |   .  |   A     |   G     |   225.417  |   .  |
|   20  |   6119441    |   .  |   A     |   G     |   225.417  |   .  |
|   21  |   42903480   |   .  |   T     |   G     |   225.417  |   .  |
|   22  |   38782696   |   .  |   A     |   G     |   225.417  |   .  |
|   3   |   4362083    |   .  |   A     |   G     |   225.417  |   .  |
|   5   |   83538811   |   .  |   T     |   C     |   222.334  |   .  |
|   8   |   120216440  |   .  |   A     |   C     |   222.406  |   .  |
|   8   |   123975238  |   .  |   T     |   C     |   222.141  |   .  |
|   9   |   76706955   |   .  |   C     |   T     |   225.417  |   .  |
|   9   |   76706991   |   .  |   TTGT  |   TT    |   9.79536  |   .  |
|   9   |   97428498   |   .  |   A     |   G     |   225.417  |   .  |
|   9   |   101371346  |   .  |   C     |   T     |   222.317  |   .  |
|   9   |   122629130  |   .  |   A     |   G     |   225.417  |   .  |
|   X   |   110451457  |   .  |   T     |   A     |   225.417  |   .  |
|   X   |   112454808  |   .  |   T     |   C     |   222.284  |   .  |
|   X   |   112454829  |   .  |   CGG   |   CGGG  |   22.3479  |   .  |

### Sample 2 - Output

|   1   |   67395837   |   .  |   C     |   A   |   222.262  |   .  |
|-------|--------------|------|---------|-------|------------|------|
|   1   |   179551371  |   .  |   G     |   A   |   184.361  |   .  |
|   10  |   98459557   |   .  |   G     |   A   |   225.417  |   .  |
|   11  |   16111867   |   .  |   A     |   G   |   104.228  |   .  |
|   12  |   884764     |   .  |   C     |   T   |   222.024  |   .  |
|   13  |   33129519   |   .  |   T     |   C   |   225.417  |   .  |
|   14  |   67575857   |   .  |   A     |   G   |   225.417  |   .  |
|   15  |   34236747   |   .  |   G     |   A   |   225.417  |   .  |
|   17  |   14102122   |   .  |   G     |   A   |   216.895  |   .  |
|   17  |   14102147   |   .  |   TTTG  |   T   |   17.3983  |   .  |
|   17  |   59886176   |   .  |   A     |   G   |   225.417  |   .  |
|   17  |   59886219   |   .  |   GGGC  |   G   |   18.7088  |   .  |
|   17  |   73200670   |   .  |   A     |   G   |   225.417  |   .  |
|   17  |   73201609   |   .  |   G     |   A   |   225.417  |   .  |
|   18  |   23833905   |   .  |   T     |   C   |   225.417  |   .  |
|   18  |   63987229   |   .  |   A     |   G   |   225.417  |   .  |
|   2   |   227032260  |   .  |   C     |   T   |   119.405  |   .  |
|   20  |   6119441    |   .  |   A     |   G   |   190.352  |   .  |
|   20  |   6119481    |   .  |   TGC   |   T   |   9.51099  |   .  |
|   20  |   54169680   |   .  |   G     |   A   |   220.803  |   .  |
|   21  |   42903480   |   .  |   T     |   G   |   107.648  |   .  |
|   22  |   38782696   |   .  |   A     |   G   |   152.405  |   .  |
|   3   |   4362083    |   .  |   A     |   G   |   215.229  |   .  |
|   3   |   4362118    |   .  |   GGAC  |   G   |   6.76067  |   .  |
|   4   |   5748177    |   .  |   T     |   C   |   225.417  |   .  |
|   4   |   5748201    |   .  |   GTTT  |   GT  |   23.1835  |   .  |
|   4   |   5748202    |   .  |   TTTA  |   T   |   12.0832  |   .  |
|   5   |   83538811   |   .  |   T     |   C   |   204.352  |   .  |
|   8   |   93923709   |   .  |   T     |   C   |   167.316  |   .  |
|   8   |   120216440  |   .  |   A     |   C   |   222.406  |   .  |
|   8   |   123975238  |   .  |   T     |   C   |   222.263  |   .  |
|   9   |   97428498   |   .  |   A     |   G   |   222.229  |   .  |
|   9   |   101371346  |   .  |   C     |   T   |   225.417  |   .  |
|   X   |   110451457  |   .  |   T     |   A   |   222.398  |   .  |
|   X   |   112454805  |   .  |   G     |   A   |   7.66682  |   .  |
|   X   |   112454808  |   .  |   T     |   C   |   59.8929  |   .  |
|   X   |   136349199  |   .  |   C     |   T   |   203.352  |   .  |

***

## GATK Code 

```

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process variant_calling_gatk {
    input:
    path bam_file

    output:
    path "${bam_file.baseName}.vcf" // Output VCF file

    script:
    """
	/Users/alireza/bioinformatics-project/gatk/gatk  --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' HaplotypeCaller \
    	-I ${bam_file} \
    	-O ${bam_file.baseName}.vcf \
    	-R /data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    	-L /data/ROI.bed
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


```

## Running Code

```
nextflow run variant_calling_gatk.nf
```

![Running code](/command-gatk.png)


## GATK Result

### Sample 1 - Output

|   1   |   67395837   |   .  |   C     |   A     |   530.06  |   .  |
|-------|--------------|------|---------|---------|-----------|------|
|   1   |   179551371  |   .  |   G     |   A     |   406.06  |   .  |
|   10  |   98459557   |   .  |   G     |   A     |   604.06  |   .  |
|   11  |   16111867   |   .  |   A     |   G     |   154.64  |   .  |
|   13  |   38859469   |   .  |   A     |   G     |   132.64  |   .  |
|   17  |   59886176   |   .  |   A     |   G     |   481.06  |   .  |
|   18  |   23833905   |   .  |   T     |   C     |   258.06  |   .  |
|   18  |   23833944   |   .  |   GTCT  |   G,GT  |   163.02  |   .  |
|   18  |   63987229   |   .  |   A     |   G     |   332.06  |   .  |
|   20  |   6119441    |   .  |   A     |   G     |   109.96  |   .  |
|   20  |   6119476    |   .  |   TCTG  |   T     |   40.71   |   .  |
|   21  |   42903480   |   .  |   T     |   G     |   515.06  |   .  |
|   22  |   38782696   |   .  |   A     |   G     |   382.06  |   .  |
|   3   |   4362083    |   .  |   A     |   G     |   361.06  |   .  |
|   5   |   83538811   |   .  |   T     |   C     |   223.64  |   .  |
|   8   |   120216440  |   .  |   A     |   C     |   251.64  |   .  |
|   8   |   123975238  |   .  |   T     |   C     |   367.64  |   .  |
|   9   |   76706955   |   .  |   C     |   T     |   332.06  |   .  |
|   9   |   97428498   |   .  |   A     |   G     |   564.06  |   .  |
|   9   |   101371346  |   .  |   C     |   T     |   80.64   |   .  |
|   9   |   122629130  |   .  |   A     |   G     |   627.06  |   .  |
|   X   |   112454808  |   .  |   T     |   C     |   98.73   |   .  |


### Sample 2 - Output

|   1   |   67395837   |   .  |   C           |   A     |   117.64  |   .  |
|-------|--------------|------|---------------|---------|-----------|------|
|   1   |   179551371  |   .  |   G           |   A     |   230.64  |   .  |
|   10  |   98459557   |   .  |   G           |   A     |   728.06  |   .  |
|   11  |   16111867   |   .  |   A           |   G     |   104.64  |   .  |
|   12  |   884764     |   .  |   C           |   T     |   392.64  |   .  |
|   14  |   67575857   |   .  |   A           |   G     |   530.06  |   .  |
|   15  |   34236747   |   .  |   G           |   A     |   307.06  |   .  |
|   17  |   14102122   |   .  |   G           |   A     |   173.64  |   .  |
|   17  |   59886176   |   .  |   A           |   G     |   307.06  |   .  |
|   17  |   73200670   |   .  |   A           |   G     |   431.06  |   .  |
|   17  |   73201609   |   .  |   G           |   A     |   357.06  |   .  |
|   17  |   73201641   |   .  |   G           |   GCTT  |   81.6    |   .  |
|   18  |   23833905   |   .  |   T           |   C     |   307.06  |   .  |
|   18  |   23833944   |   .  |   G           |   GA    |   51.6    |   .  |
|   18  |   63987229   |   .  |   A           |   G     |   382.06  |   .  |
|   2   |   227032260  |   .  |   C           |   T     |   104.64  |   .  |
|   20  |   54169680   |   .  |   G           |   A     |   160.64  |   .  |
|   20  |   54169706   |   .  |   GCAAAAGACA  |   G     |   167.92  |   .  |
|   3   |   4362083    |   .  |   A           |   G     |   151.64  |   .  |
|   4   |   5748177    |   .  |   T           |   C     |   703.06  |   .  |
|   5   |   83538811   |   .  |   T           |   C     |   142.64  |   .  |
|   8   |   93923709   |   .  |   T           |   C     |   136.64  |   .  |
|   8   |   120216440  |   .  |   A           |   C     |   173.64  |   .  |
|   8   |   123975238  |   .  |   T           |   C     |   477.64  |   .  |
|   9   |   97428498   |   .  |   A           |   G     |   229.64  |   .  |
|   9   |   101371346  |   .  |   C           |   T     |   589.06  |   .  |
|   X   |   112454808  |   .  |   T           |   C     |   136.64  |   .  |
|   X   |   136349199  |   .  |   C           |   T     |   185.64  |   .  |


## Why bcftools and GATK showing different results

The discrepancy between bcftools and GATK showing different results regarding variant calls can arise from several factors. Here are some common reasons and considerations to help you troubleshoot the situation:

**1. Different Calling Algorithms**
bcftools and GATK use different algorithms and methodologies for variant calling. GATK is often considered more sophisticated due to its focus on germline variant discovery and can include advanced features like joint calling across samples.
bcftools is generally simpler and faster for quick calls but might miss variants that GATK identifies, especially in low-coverage regions or in the presence of complex variants.

**2. Default Parameters and Settings**
The default parameters used in bcftools and GATK may differ. For example, the minimum quality score or the filtering criteria applied by default can affect what is considered a variant.
Check the specific parameters you’re using in both tools. For example, ensure you’re not applying stricter filters in bcftools that might exclude variants.

**3. Coverage and Read Depth**
If your BAM file has low coverage or poor quality reads, bcftools may not call variants that GATK finds, especially if GATK uses advanced techniques like local realignment or base quality score recalibration.
