# Heterozygosity tutorial

Heterozygosity is a simple and informative statistic that can be obtained by analyzing whole-genome data. You can calculate the average heterozygosity of an individual or assess the local heterozygosity throughout the genome to answer more specific questions.

There are many methods available to achieve this. In this tutorial, we will use two strategies. To obtain the average heterozygosity, we will manipulate the VCF file using bcftools and vcftools. To explore heterozygosity variation across the genome, we will use ANGSD and other associated tools.

### Average heterozygosity

```bash
#!/bin/bash

#Replace these with the actual file names
VCF_FILE="/pool/genomics/figueiroh/SMSC_2023_class/vcf/NN_6samples_HD_PASS_DP5.vcf.gz"
OUTPUT_FILE="NN_heterozygosity_v5.tsv"
GENOME_LENGTH=2468345093 #you can use the fasta index (.fai) to sum the total length of the genome

#Get a list of sample names from the VCF file
SAMPLES=$(bcftools query -l $VCF_FILE)

#Write a header line to the output file
echo -e "Sample\tHeterozygous_sites\tHeterozygosity" > $OUTPUT_FILE

#Loop through each sample and calculate the heterozygosity
for SAMPLE in $SAMPLES; do
  HETEROZYGOUS=$(bcftools view -s $SAMPLE $VCF_FILE | grep -v "#" | grep -o "0/1" | wc -l)
  HETEROZYGOSITY=$(echo "scale=7; $HETEROZYGOUS / $GENOME_LENGTH" | bc)
  echo -e "$SAMPLE\t$HETEROZYGOUS\t$HETEROZYGOSITY" >> $OUTPUT_FILE
done

```

This is a bash script that calculates the heterozygosity of each sample in a VCF file. It uses the bcftools command line tool to extract the list of sample names and calculate the number of heterozygous variants for each sample. The results are written to a tab-separated file with the name "heterozygosity.tsv" and two columns: "Sample" and "Heterozygosity". The script assumes that the VCF file is named "your_vcf_file.vcf" and is located in the same directory as the script. Copy this script to a file, give a name that ends with `.sh`, and run on interactive mode on Hydra using `bash command.sh`.

| Individual | Heterozygous sites | Heterozygosity |
| --- | --- | --- |
| NN114296 | 1690684 | 0.0006762 |
| NN114297 | 1690576 | 0.0006762 |
| NN114393 | 1599970 | 0.0006399 |
| NN115950 | 1645258 | 0.0006581 |
| NN190240 | 1413972 | 0.0005655 |

```r
library(ggplot2)

# Replace this with the actual file name
INPUT_FILE <- "heterozygosity.tsv"

# Read in the data from the input file
data <- read.table(INPUT_FILE, header=TRUE, sep="\t")

# Create a scatter plot of heterozygosity values
plot <- ggplot(data, aes(x=Sample, y=Heterozygosity)) +
    geom_point(size=2) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
    labs(x="Sample", y="Heterozygosity")

# Save the plot to a file
ggsave("heterozygosity_plot.png", plot, width=10, height=5, dpi=300)
```

This is an R code block that creates a scatter plot of heterozygosity values from a tab-separated value (tsv) file. The ggplot2 library is loaded, and the name of the input file is specified as "heterozygosity.tsv". The data is read in using the read.table function, and a scatter plot is created using ggplot with "Sample" on the x-axis and "Heterozygosity" on the y-axis. The resulting plot is saved as a png file named "heterozygosity_plot.png".

### Genome-wide heterozygosity

Another way to display heterozygosity is by using a window-based approach. This approach allows for the visualization of variation throughout the genome and the ability to focus on regions of particularly high or low diversity. It also permits comparison of results with other methods, such as Runs of Homozygosity (RoH). There are several software programs that can perform this task in a similar way, including ANGSD.

ANGSD (Analysis of Next Generation Sequencing Data) is an open-source software designed to analyze low-depth Next Generation Sequencing (NGS) data. Developed by researchers at the University of Copenhagen, ANGSD was specifically designed to handle large-scale sequencing data, particularly from non-model organisms. ANGSD provides a flexible and efficient framework for analyzing complex genomic data, enabling researchers to perform various types of population genetics and genomics analyses.

Key capabilities of ANGSD include:

1. Genotype calling: ANGSD can estimate genotype likelihoods from sequencing data without actually calling genotypes, which helps reduce biases and errors introduced by hard genotype calls.
2. SNP discovery: The software can identify Single Nucleotide Polymorphisms (SNPs) and estimate their allele frequencies while accounting for sequencing errors and varying levels of sequencing depth.
3. Genetic association studies: ANGSD can perform Genome-Wide Association Studies (GWAS) and estimate genotype-phenotype associations using mixed linear models or generalized linear models.
4. Population structure: ANGSD can estimate population structure and admixture proportions using principal component analysis (PCA) or model-based approaches.
5. Selection scans: The software can detect signals of positive or balancing selection, using various statistics like FST, Tajima's D, and nucleotide diversity.
6. Handling of various file formats: ANGSD can work with various input file formats, including BAM, CRAM, and VCF files, and output results in multiple formats suitable for downstream analyses.

Overall, ANGSD is a versatile and powerful tool for analyzing NGS data, particularly for non-model organisms or low-coverage sequencing projects. It provides a comprehensive suite of analysis tools that cater to various research objectives in population genetics and genomics.

### How to run heterozygosity using ANGSD

```bash
module load bioinformatics/angsd/0.921

angsd -P <threads> -i <input_bam_file> -anc <ancestral_fasta_file> -dosaf <dosaf_value> -gl <genotype_likelihood_method> -C <base_quality_adjustment> -minQ <min_base_quality> -minmapq <min_mapping_quality> -fold <fold_value> -out <output_file> -ref <reference_fasta_file> -r <region_of_interest>
```

- `P <threads>` - Sets the number of threads to be used in parallel.
- `i <input_bam_file>` - Specifies the input file as a BAM file.
- `anc <ancestral_fasta_file>` - Specifies the ancestral fasta reference file.
- `dosaf <dosaf_value>` - Computes the Site Frequency Spectrum (SFS) based on the genotype likelihoods.
- `gl <genotype_likelihood_method>` - Specifies the method used for calculating genotype likelihoods.
- `C <base_quality_adjustment>` - Adjusts the base quality score by a specified value before using it.
- `minQ <min_base_quality>` - Sets the minimum base quality score required.
- `minmapq <min_mapping_quality>` - Sets the minimum mapping quality score required.
- `fold <fold_value>` - Indicates whether you are analyzing folded SFS or unfolded SFS.
- `out <output_file>` - Specifies the output file path and name.
- `ref <reference_fasta_file>` - Specifies the reference fasta file.
- `r <region_of_interest>` - Specifies the region of interest.

```bash
realSFS -nsites <number_of_sites> <input_saf_idx_file> > <output_est_ml_file>
```

- `nsites <number_of_sites>` - Specifies the number of sites to be considered for the estimation. Replace `<number_of_sites>` with the desired number.
- `<input_saf_idx_file>` - Specifies the input .saf.idx file, which is the output from the ANGSD program that contains information about the SFS.
- `> <output_est_ml_file>` - Specifies the output file path and name for the maximum likelihood estimate of the SFS.

ANGSD does not differentiate the chromosomes if you run the whole genome at once, that is why we need to use the ‘region’ variable when running ANGSD to specify the chromosome/scaffold. This will be useful to look the heterozygosity throughout the genome. 

```bash
#!/bin/bash

#Set variables

input_bam_file="/pool/genomics/figueiroh/SMSC_2023/mapping/NN114296_cloud_leopard_sorted.bam"
ancestral_fasta_file="/pool/genomics/figueiroh/SMSC_2023/reference/mNeoNeb1.pri.cur.20220520.fasta"
reference_fasta_file="/pool/genomics/figueiroh/SMSC_2023/reference/mNeoNeb1.pri.cur.20220520.fasta"
output_directory="/pool/genomics/figueiroh/SMSC_2023/heterozygosity"
SAMPLE="NN114296"

#Loop through scaffolds 1 to 18

for i in {1..18}; do
# Run ANGSD command
angsd -P 10 -i ${input_bam_file} -anc ${ancestral_fasta_file} -dosaf 1 -gl 1 -C 50 -minQ 20 -minmapq 30 -fold 1 -out ${output_directory}/$SAMPLE.SUPER_${i} -ref ${reference_fasta_file} -r SUPER_${i}

# Run realSFS command
realSFS -nsites 200000 ${output_directory}/$SAMPLE.SUPER_${i}.saf.idx > ${output_directory}/$SAMPLE.SUPER_${i}.est.ml

done
```

Replace the `<placeholders>` with the desired values for your specific analysis, and update the paths for input and output files as needed. The `$SAMPLE` variable should also be set to the appropriate sample name.

This script will loop through scaffolds 1 to 19, running the ANGSD command and then the realSFS command for each scaffold. The results will be saved in the specified output directory.

Now we can add the sample name and scaffold number for each line in our output. This will make our work easier when we want to plot our results. 

```bash
#!/bin/bash

output_directory="path/to/your/output_directory"
SAMPLE="your_sample_name"

# Loop through scaffolds 1 to 18
for i in {1..18}; do
  # Calculate the number of lines in the output file
  num_lines=$(wc -l < ${output_directory}/$SAMPLE.SUPER_${i}.est.ml)

  # Add the number of lines, sample name, and scaffold number to each line of the output file
  awk -v lines="$num_lines" -v sample="$SAMPLE" -v scaffold="$i" '{print lines, sample, "SUPER_" scaffold, $0}' ${output_directory}/$SAMPLE.SUPER_${i}.est.ml > ${output_directory}/$SAMPLE.SUPER_${i}.est.ml.annotated

  # Optional: Move the annotated file to the original file
  mv ${output_directory}/$SAMPLE.SUPER_${i}.est.ml.annotated ${output_directory}/$SAMPLE.SUPER_${i}.est.ml

done
```

Concatenate files

```bash
#!/bin/bash

Set the input and output directory

input_directory="path/to/your/input_directory"
output_directory="path/to/your/output_directory"

Create a new output file

output_file="${output_directory}/all_est_ml_concatenated.txt"

Remove the output file if it already exists

if [ -f "${output_file}" ]; then
rm "${output_file}"
fi

Loop through all "[est.ml](http://est.ml/)" files and concatenate them
for file in ${input_directory}/*.est.ml; do
cat "${file}" >> "${output_file}"
done
```

Plot the results for one chromosome using R

```bash
# Load required libraries
library(tidyverse) # Collection of packages for data manipulation and visualization
library(viridis)   # Package for generating color palettes
library(scales)    # Package for scaling and formatting axis labels

# Read data file and store it in the variable 'het_master'
het_master <- read.table("/path/to/file/NN_SUPER.est.ml")

# Manipulate the data
het_master %>% 
  rename(sample=V2) %>%          # Rename V2 as 'sample'
  rename(chromosome = V3) %>%    # Rename V3 as 'chromosome'
  mutate(heterozygosity = V5/(V4 + V5)) %>% # Calculate heterozygosity as V5 / (V4 + V5)
  mutate(position = ((V1*200000)-200000))   %>% # Calculate position as (V1 * 200000) - 200000
  filter(chromosome == "SUPER_2") %>% # Filter data to include only rows where 'chromosome' is 'SUPER_2'

  # Create a ggplot2 plot
  ggplot(aes(x=position, y=heterozygosity)) + # Set x-axis as 'position' and y-axis as 'heterozygosity'
  geom_line(colour="grey",alpha=0.5) + # Add a line plot with grey color and 0.5 alpha (transparency)
  geom_point(aes(colour=factor(chromosome))) + # Add points, color them based on 'chromosome' factor variable
  scale_color_viridis(discrete = TRUE) + # Use viridis color palette for the points
  facet_grid(sample ~ chromosome,scales = "free_x") + # Create a facet grid with 'sample' on the y-axis and 'chromosome' on the x-axis, set x-axis scales to be free
  labs(x = NULL, y = "Heterozygosity\n") + # Remove x-axis labels and set y-axis label to "Heterozygosity\n"
  scale_y_continuous(labels = comma) + # Format y-axis labels with commas
  scale_x_continuous(labels = comma) + # Format x-axis labels with commas
  theme_minimal() + # Apply a minimal theme to the plot
  theme(legend.position = "none", # Remove legend
        strip.text.x = element_text(face = "bold"), # Set strip text for x-axis to bold
        strip.text.y = element_text(face = "bold"), # Set strip text for y-axis to bold
        panel.grid.major.x = element_blank(), # Remove major x-axis gridlines
        panel.grid.minor.x = element_blank(), # Remove minor x-axis gridlines
        panel.spacing.x = unit(0, "line"), # Set panel spacing to zero
        panel.border = element_rect(color = "black", fill = NA, size = 0.25)) # Add a black border around the panels
```

### How to run on an HPC

```bash
#!/bin/sh
#---------------Parameters----------------------#
#$ -S /bin/sh
#$ -pe mthread 4
#$ -q mThC.q
#$ -l mres=4G,h_data=1G,h_vmem=1G
#$ -cwd
#$ -j y
#$ -N ANGSD_realSFS_job
#$ -o ANGSD_realSFS_job.log
#---------------Modules-------------------------#
module load bioinformatics/angsd/0.921
module load other_required_modules
#---------------Your Commands-------------------#
echo + date job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME

#Set variables

input_bam_file="path/to/your/input_bam_file.bam"
ancestral_fasta_file="path/to/your/ancestral_fasta_file.fasta"
reference_fasta_file="path/to/your/reference_fasta_file.fasta"
output_directory="path/to/your/output_directory"
SAMPLE="your_sample_name"

Loop through scaffolds 1 to 19

for i in {1..18}; do
# Run ANGSD command
angsd -P <threads> -i ${input_bam_file} -anc ${ancestral_fasta_file} -dosaf <dosaf_value> -gl <genotype_likelihood_method> -C <base_quality_adjustment> -minQ <min_base_quality> -minmapq <min_mapping_quality> -fold <fold_value> -out ${output_directory}/$SAMPLE.scaffold${i} -ref ${reference_fasta_file} -r HiC_scaffold_${i}

# Run realSFS command
realSFS -nsites <number_of_sites> ${output_directory}/$SAMPLE.scaffold${i}.saf.idx > ${output_directory}/$SAMPLE.scaffold${i}.est.ml

done

echo = date job $JOB_NAME done
```

```bash
#!/bin/sh
#---------------Parameters----------------------#
#$ -S /bin/sh
#$ -pe mthread 4
#$ -q mThC.q
#$ -l mres=4G,h_data=1G,h_vmem=1G
#$ -cwd
#$ -j y
#$ -N Annotate_EST_ML_job
#$ -o Annotate_EST_ML_job.log
#---------------Modules-------------------------#
#Load any required modules if necessary (uncomment the following line and replace 'module_name')
#module load module_name
#---------------Your Commands-------------------#
echo + date job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#Set variables

output_directory="path/to/your/output_directory"
SAMPLE="your_sample_name"

#Loop through scaffolds 1 to 19

for i in {1..19}; do
# Add sample name and scaffold number to each line of the output file
awk -v sample="$SAMPLE" -v scaffold="$i" '{print sample, "scaffold" scaffold, $0}' ${output_directory}/$SAMPLE.scaffold${i}.est.ml > ${output_directory}/$SAMPLE.scaffold${i}.est.ml.annotated

# Optional: Move the annotated file to the original file
mv ${output_directory}/$SAMPLE.scaffold${i}.est.ml.annotated ${output_directory}/$SAMPLE.scaffold${i}.est.ml

done

echo = date job $JOB_NAME done
```
