# SNP density tutorial

SNP density is a valuable tool for visualizing genetic diversity across the genome. This technique involves representing the number of SNPs in a heatmap, which is particularly useful for comparing genetic diversity between different populations or closely related species.

To calculate SNP density, count the number of SNPs over a certain length of the genome, typically ranging from 100 kilobases to 1 megabase depending on the diversity and analyzed region. One way to perform this analysis is by using vcftools with the snpden function (Danecek et al., 2011). This function takes a VCF file as input and outputs a table with the SNP count for each region interval.

Finally, the results can be visualized by generating a plot using your preferred programming language. In this case, we will use the R programming language with the tidyverse package.

The VCF file used for this tutorial can be found in the folder below

```bash
/pool/genomics/figueiroh/SMSC2023_class/vcf
```

The first step is to create the folder where we are going to work with the vcf files in order to generate the SNP density plots.

```bash
mkdir snpden

cd snpden

mkdir plots
```

Inside the `snpden` folder, the first step is to create new VCF files that only contain heterozygous positions. There are several ways to accomplish this, but for simplicity, we will use vcftools.

```bash
vcftools --gzvcf NN114296_HD_PASS_DP5.vcf.gz --recode --out NN114296_HD_PASS_DP5_hetsites --maf 0.1

vcftools --gzvcf NN114297_HD_PASS_DP5.vcf.gz --recode --out NN114297_HD_PASS_DP5_hetsites --maf 0.1

vcftools --gzvcf NN114393_HD_PASS_DP5.vcf.gz --recode --out NN114393_HD_PASS_DP5_hetsites --maf 0.1

vcftools --gzvcf NN115950_HD_PASS_DP5.vcf.gz --recode --out NN115950_HD_PASS_DP5_hetsites --maf 0.1

vcftools --gzvcf NN190240_HD_PASS_DP5.vcf.gz --recode --out NN190240_HD_PASS_DP5_hetsites --maf 0.1
```

1. `-gzvcf input.vcf.gz`: Specifies the input VCF file in compressed gzip format, located in the parent directory's "filtered" folder. The file is named "NN114296_HD_PASS_DP5.vcf.gz".
2. `-recode`: Indicates that you want to create a new VCF file with the filtered results.
3. `-out outuput`: Specifies the output file prefix, which will be "NN114296_HD_PASS_DP5_hetsites". The output file will be named "NN114296_HD_PASS_DP5_hetsites.recode.vcf" in the current directory.
4. `-maf 0.1`: Sets a minimum allele frequency (MAF) filter of 0.1, meaning that only variants with a MAF of at least 10% will be included in the output file.

Now we can run the snp density function on vcftools. This will count the number of heteroyzgous snps breaking the genome into windows. 

```bash
vcftools --vcf NN114296_HD_PASS_DP5_hetsites.recode.vcf --SNPdensity 1000000 --out NN114296_HD_PASS_DP5_hetsites
```

1. `-vcf NN114296_HD_PASS_DP5_hetsites.recode.vcf`: Specifies the input VCF file named "NN114296_HD_PASS_DP5_hetsites.recode.vcf" in the current directory.
2. `-SNPdensity 1000000`: Calculates the SNP density in non-overlapping windows of 1,000,000 base pairs (1 Mb) across the genome. SNP density is the number of SNPs per window.
3. `-out NN114296_HD_PASS_DP5_hetsites`: Specifies the output file prefix, which will be "NN114296_HD_PASS_DP5_hetsites". The output file will be named "NN114296_HD_PASS_DP5_hetsites.snpden" in the current directory.

You can do the same thing for the other available samples

```bash
vcftools --vcf NN114297_HD_PASS_DP5_hetsites.recode.vcf --SNPdensity 1000000 --out NN114297_HD_PASS_DP5_hetsites

vcftools --vcf NN114393_HD_PASS_DP5_hetsites.recode.vcf --SNPdensity 1000000 --out NN114393_HD_PASS_DP5_hetsites

vcftools --vcf NN115950_HD_PASS_DP5_hetsites.recode.vcf --SNPdensity 1000000 --out NN115950_HD_PASS_DP5_hetsites

vcftools --vcf NN190240_HD_PASS_DP5_hetsites.recode.vcf --SNPdensity 1000000 --out NN190240_HD_PASS_DP5_hetsites
```

Here is an example bash loop that can be used to run the `vcftools` commands on multiple input files:

```
for file in *.vcf.gz
do
    name=$(basename $file)
    name=${name%.vcf.gz}
    vcftools --vcf $file --SNPdensity 1000000 --out ${name}_density
done

```

This loop will iterate through all files in the directory that end in ".vcf.gz", extract the base name and use it to name the output files. It will then run the same `vcftools` command on each input file, breaking the genome into 1-megabase windows, and outputting the results to a file with the format "input_file_density".

Then we add the sample names to each table

```bash
awk -v sample="NN114296" 'NR==1{print $0"\tIndiv"} NR>1{print $0"\t"sample}' NN114296_HD_PASS_DP5_hetsites.snpden > NN114296_HD_PASS_DP5_hetsites_id.snpden
```

1. `v sample="NN114296"`: Sets a variable called `sample` with the value "NN114296".
2. `'NR==1{print $0"\tIndiv"} NR>1{print $0"\t"sample}'`: This is an `awk` script that processes the input file line by line. It checks if the current line number (NR) is equal to 1, and if it is, it prints the entire line (`$0`) followed by a tab and the string "Indiv". For all other lines (NR > 1), it prints the entire line followed by a tab and the value of the `sample` variable.
3. `NN114296_HD_PASS_DP5_hetsites.snpden`: Specifies the input file, which is the SNP density output file generated by VCFtools in your previous command.
4. `> NN114296_HD_PASS_DP5_hetsites_id.snpden`: Redirects the output to a new file named "NN114296_HD_PASS_DP5_hetsites_id.snpden" in the current directory.

```bash
awk -v sample="NN114297" 'NR==1{print $0"\tIndiv"} NR>1{print $0"\t"sample}' NN114297_HD_PASS_DP5_hetsites.snpden > NN114297_HD_PASS_DP5_hetsites_id.snpden 

awk -v sample="NN114393" 'NR==1{print $0"\tIndiv"} NR>1{print $0"\t"sample}' NN114393_HD_PASS_DP5_hetsites.snpden > NN114393_HD_PASS_DP5_hetsites_id.snpden

awk -v sample="NN115950" 'NR==1{print $0"\tIndiv"} NR>1{print $0"\t"sample}' NN115950_HD_PASS_DP5_hetsites.snpden > NN115950_HD_PASS_DP5_hetsites_id.snpden

awk -v sample="NN190240" 'NR==1{print $0"\tIndiv"} NR>1{print $0"\t"sample}' NN190240_HD_PASS_DP5_hetsites.snpden > NN190240_HD_PASS_DP5_hetsites_id.snpden
```

Concatenate all samples into a single file

```bash
tail -q -n +2 *_id.snpden > NN_hetsites.snpden
```

Here is a possible bash script that consolidates the code and makes it more efficient:

```
#!/bin/bash

# create directories
mkdir -p snpden/plots

# set parameters
WINDOW_SIZE=1000000
MAF_THRESHOLD=0.1
SAMPLES=(NN114296 NN114297 NN114393 NN115950 NN190240)

# process input files
for sample in "${SAMPLES[@]}"; do
  # extract heterozygous sites
  vcftools --gzvcf filtered/${sample}_HD_PASS_DP5.vcf.gz --recode --out snpden/${sample}_hetsites --maf ${MAF_THRESHOLD}

  # calculate SNP density
  vcftools --vcf snpden/${sample}_hetsites.recode.vcf --SNPdensity ${WINDOW_SIZE} --out snpden/${sample}_hetsites

  # add sample name to output
  awk -v sample="${sample}" 'NR==1{print $0"\\tIndiv"} NR>1{print $0"\\t"sample}' snpden/${sample}_hetsites.snpden > snpden/${sample}_hetsites_id.snpden
done

# combine output files
tail -q -n +2 snpden/*_id.snpden > snpden/NN_hetsites.snpden

```

This script automates the process of extracting heterozygous sites, calculating SNP density, and adding sample names to the output files. The script also uses variables to set the window size and minor allele frequency threshold, so these values can easily be changed if needed. Finally, the script combines the output files into a single file, which can be used to generate the plot.

You can now download the output to your machine and use RStudio to plot the results.

There are two scripts available for this task. You can find the scripts to build the plots for this tutorial in the Github folder.

`Neofelis_snpdenPlot_chr.R`

This script is written in R and does the following tasks:

1. Sets the working directory to the folder where the data files are located.
2. Loads the required packages (tidyverse).
3. Reads the SNP density data file (NN_hetsites.snpden.txt) and stores it in a variable called 'snpden'.
4. Defines the order of the scaffolds and the chromosomes to be used in the visualization.
5. Reorders the chromosome column of the data frame according to the target order.
6. Subsets the data from chromosomes that are not "NA".
7. Categorizes the SNP density data based on the number of variants per kilobase.
8. Renames the levels of the CHROM column to match the defined chromosome order.
9. Iterates through each chromosome in the dataset and performs the following steps:
a. Subsets the data for the current chromosome.
b. Creates a ggplot object with the SNP density data and applies appropriate formatting, colors, and facet grid.
c. Saves the generated plot as a .png file in the "plots/" folder, with the file name including the chromosome name (e.g., 'Neofelis_chr1.1Mb.snpden.png').

This script creates a series of SNP density plots for each chromosome in the dataset, showing the distribution of heterozygous SNPs along the chromosome length. The plots are saved as .png files in the "plots/" folder.

`Neofelis_snpdenPlot_sample.R`

1. Set the working directory to the location where the data files are stored.
2. Load the required packages, specifically the tidyverse package.
3. Read in the SNP density data file (NN_hetsites.snpden.txt) using the read.table function.
4. Define the order of the scaffolds to be used in the visualization and the order of the chromosomes to be used in the visualization.
5. Copy the data frame to another variable snpden.master.
6. Reorder the chromosome column of the data frame according to the target order.
7. Subset the data from chromosomes that are not "NA."
8. Create a new column called groups that cuts the numeric values in the VARIANTS.KB column into different groups based on their values.
9. Rename the CHROM levels to have chromosome names that start with "chr."
10. Convert the BIN_START column from a factor to a numeric column.
11. Loop over each individual in the dataset.
12. Subset the data for the current individual.
13. Define a title for the plot.
14. Create a ggplot object using the ggplot function, specifying the data frame and the aesthetic mapping.
15. Add a tile layer to the plot using geom_tile to represent the SNP density.
16. Facet the plot by chromosome and scaffolds.
17. Set the axis labels and plot titles using the labs function.
18. Customize the plot using theme functions, such as theme_minimal, element_text, and element_blank.
19. Scale the fill color using the scale_fill_manual function.
20. Scale the x-axis using the scale_x_continuous function.
21. Save the plot using the ggsave function.