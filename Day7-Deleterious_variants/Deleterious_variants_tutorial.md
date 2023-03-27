# Deleterious variants analyses

### Introduction

SnpEff and SnpSift are powerful tools designed for annotating and analyzing genetic variants, with a focus on deleterious variants. This tutorial will guide you through the process of using these tools to work with non-model organisms. We'll provide insights on how to map the genome against a reference that has a SnpEff database available, perform the analysis, and finish with enrichment analysis options using g:Profiler.

1. Preparing the reference genome and annotation

Before using SnpEff and SnpSift, you need to have a reference genome and annotation files. If you're working with non-model organisms, you might have to create a custom SnpEff database using the organism's reference genome and annotation files. The annotation files usually include gene models in GFF3, GTF, or Gencode format. Another option is to map your reads to a closely related species that already has its genome in the SnpEff database.

1. Creating a custom SnpEff database

To create a custom SnpEff database, follow these steps:

a. Download the SnpEff software and set up the environment:

```bash
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
```

b. Create a new directory for your non-model organism in the `data` folder:

```bash
mkdir data/your_organism
```

c. Copy the reference genome (FASTA) and annotation (GFF3, GTF, or Gencode) files into the new directory:

```bash
cp /path/to/reference_genome.fasta data/your_organism/
cp /path/to/annotation_file.gff3 data/your_organism/
```

d. Modify the `snpEff.config`  file to include the new genome:

```bash
echo "your_organism.genome : Your_Organism" >> snpEff.config
```

e. Build the SnpEff database:

```bash
java -jar snpEff.jar build -gff3 -v your_organism
```

1. Mapping the genome against a reference with SnpEff database

To map the genome of your non-model organism against a reference with a SnpEff database, you will first need to align the reads to the reference genome using an aligner like BWA or Bowtie2. Once you have the alignments in BAM or SAM format, you can call variants using a tool like GATK, FreeBayes, or Samtools.

1. Annotating variants with SnpEff

With the custom SnpEff database created and the variants called, you can now annotate the variants using SnpEff:

```bash
java -jar snpEff.jar your_organism input.vcf > annotated_output.vcf
```

1. Analyzing deleterious variants with SnpSift

SnpSift is a collection of tools that can be used to filter, sort, and analyze the annotated VCF files. You can filter deleterious variants using SnpSift Filter:

```bash
java -jar SnpSift.jar filter "ANN[0].EFFECT has 'missense_variant' | ANN[0].EFFECT has 'frameshift_variant'" annotated_output.vcf > deleterious_variants.vcf
```

Some other examples of filters you can use with SnpSift

a) Filtering by quality (QUAL) and depth (DP):

```bash
java -jar SnpSift.jar filter "(QUAL > 30) & (DP > 10)" input.vcf > filtered_output.vcf
```

b) Filtering by minor allele frequency (AF) for a specific population in a 1000 Genomes Project VCF file:

```bash
java -jar SnpSift.jar filter "AF < 0.05" input.vcf > filtered_output.vcf
```

c) Filtering by impact, retaining only HIGH impact variants:

```bash
java -jar SnpSift.jar filter "ANN[0].IMPACT has 'HIGH'" input.vcf > high_impact_output.vcf
```

d) Filtering by specific gene:

```bash
java -jar SnpSift.jar filter "ANN[0].GENE has 'BRCA1'" input.vcf > brca1_output.vcf
```

e) Filtering variants that are either missense, frameshift, or stop gained, and have a SIFT score below 0.05:

```bash
java -jar SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant' | ANN[0].EFFECT has 'frameshift_variant' | ANN[0].EFFECT has 'stop_gained') & (ANN[0].SIFT_SCORE < 0.05)" input.vcf > filtered_output.vcf
```

### Enrichment analysis with g:Profiler

Performing gene enrichment analysis is a crucial step after identifying genes containing deleterious variants, as it provides essential insights into the biological context and functional implications of these genes. Gene enrichment analysis allows researchers to understand the roles of genes within a broader biological framework, facilitating the elucidation of underlying molecular mechanisms and pathways associated with a particular phenotype or disease.

g:Profiler is a web-based tool that allows you to perform functional enrichment analysis. It also has a python package, you can see an script example here

```python
import argparse
import pandas as pd
from gprofiler import GProfiler

def read_gene_list(file_path):
    with open(file_path, "r") as file:
        genes = [line.strip() for line in file.readlines()]
    return genes

def run_gprofiler(genes, output_prefix):
    gp = GProfiler(return_dataframe=True)

    enrichment_results = gp.profile(organism='hsapiens', query=genes)

    # Save the raw results to a CSV file
    enrichment_results.to_csv(f"{output_prefix}_raw_enrichment_results.csv", index=False)

    return enrichment_results

def parse_and_save_results(enrichment_results, output_prefix):
    filtered_results = enrichment_results[
        (enrichment_results["p_value"] <= 0.05) &
        (enrichment_results["source"].isin(["GO:BP", "GO:MF", "GO:CC", "KEGG", "HP"]))
    ]
    
    sorted_results = filtered_results.sort_values("p_value")

    # Save all significant categories to a new CSV file
    sorted_results.to_csv(f"{output_prefix}_significant_enrichment_categories.csv", index=False)

    # Save all significant categories to a new Excel file
    sorted_results.to_excel(f"{output_prefix}_significant_enrichment_categories.xlsx", index=False)

    return sorted_results

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Run g:Profiler and parse enrichment analysis results")
parser.add_argument("gene_list_file", help="Input file containing the list of genes")
parser.add_argument("output_prefix", help="Output file prefix for significant enrichment categories")

args = parser.parse_args()

# Read the gene list from the input file
genes = read_gene_list(args.gene_list_file)

# Run g:Profiler and save the raw results
enrichment_results = run_gprofiler(genes, args.output_prefix)

# Parse the results and save significant categories
significant_categories = parse_and_save_results(enrichment_results, args.output_prefix)

# Display significant categories
print(significant_categories)
```
