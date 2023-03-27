# Runs of homozygosity tutorial

Runs of Homozygosity (ROH) are defined as the uninterrupted stretches of homozygous genotypes within an individual's genome. These regions can provide valuable insights into the demographic history, inbreeding levels, and disease susceptibility of a population. By analyzing the length and distribution of ROH, we can infer the population structure, migration patterns, and effective population size of a group.

Moreover, ROH can also be used to identify deleterious mutations and genomic regions under positive selection. Inbreeding depression, which is caused by the accumulation of deleterious alleles, can be estimated by measuring the frequency and length of ROH. Longer ROH segments are associated with increased homozygosity and reduced genetic diversity, which can lead to reduced fitness and increased risk of disease.

In this tutorial, we will demonstrate how to estimate ROH using the 'bcftools roh' plugin, which is a widely used tool for detecting ROH from VCF files. We will also discuss the interpretation and application of ROH results in different research contexts, such as conservation genetics, human population genetics, and animal breeding. By the end of this tutorial, you will have a better understanding of the biological significance and practical utility of ROH analysis.

Please see below the location of the VCF file that we will use for the tutorial.

```r
/pool/genomics/figueiroh/SMSC_2023_class/vcf/NN_6samples_HD_PASS_DP5.vcf.gz
```

For the RoHs plots, we are only interested in the “RG” portion of the files, where it contains the homozygous blocks in the genome. These blocks are important because they are indicative of long stretches of DNA that are identical in the two chromosomes, which can occur when the parents are related.

There are several ways of showing the results, and it will mostly depend on your main question. For example, if you are interested in the frequency of RoHs in different populations, you can create histograms that show the distribution of the length of these blocks. On the other hand, if you want to study the relationship between RoHs and disease, you may want to compare the number and length of RoHs between cases and controls, and perform statistical tests to determine if there is an association. In either case, it is important to consider the study design and the underlying biological mechanisms that could affect the results.

```bash
module load bioinformatics/bcftools

bcftools roh --AF-dflt 0.4 -I -G30 --rec-rate 1.4e-9 /pool/genomics/figueiroh/SMSC_2023/vcf/NN_6samples_HD_PASS_DP5.vcf.gz > /path/to/your/folder/NN_6samples_HD_PASS_DP5.roh.txt
```

- `bcftools roh`: This command runs the `roh` plugin from `bcftools` to detect runs of homozygosity.
- `-AF-dflt 0.4`: This option sets the default allele frequency to 0.4 when the frequency is missing in the input file.
- `I`: This option tells the plugin to perform the imputation of missing genotypes.
- `G30`: This option sets the phred-scaled genotype quality threshold to 30. Genotypes below this quality threshold will be treated as missing.
- `-rec-rate 1.4e-9`: This option sets the recombination rate to 1.4e-9.

After running the 'bcftools roh' plugin, you may want to filter and process the output file to retain specific information. For example, you can extract the chromosome, start, and end positions of ROH using the following command:

```bash
grep "RG" NN_6samples_HD_PASS_DP5.roh.txt | cut -f 2,3,6 > NN_6samples_HD_PASS_DP5.roh.edited.txt
```

You can run this command on interactive mode. 

Two of the most basic statistics we can obtain from this analysis are the number of runs of homozygosity blocks (NROH) and the total length of genome that showed runs of homozygosity (SROH). These two values provide highly informative data, and the ratio between them is known as the inbreeding coefficient (FROH).

You can use the following R script to estimate these values and plot the results.

```r
# Set the working directory
setwd("/Users/henrique/Dropbox/Documentos/Pós-Doc/Smithsonian/SMSC_workshop/roh")

# Load libraries and read data
library(tidyverse)
library(ggrepel)

# Read data with read_delim() for better control over input file parsing
clouded_roh <- read_delim("NN_6samples_HD_PASS_DP5.roh.edited.txt", delim = "\t", skip = 1, col_names = c("Sample", "Chromosome", "RoH_length"))

# Compute NROH and SROH
clouded_nroh <- clouded_roh %>% 
  group_by(Sample) %>% 
  summarize(NROH = n())

clouded_sroh <- clouded_roh %>% 
  group_by(Sample) %>% 
  summarize(SROH = sum(RoH_length))

# Compute FROH
clouded_froh <- inner_join(clouded_nroh, clouded_sroh, by = "Sample") %>% 
  mutate(FROH = NROH / SROH)

# Create a table with NROH, SROH, and FROH for each sample
summary_table <- clouded_froh

# Display the table
print(summary_table)

# Save the table to a CSV file
write_csv(summary_table, "summary_table.csv")

# Plot NROH vs. SROH and save the plot to a file
froh_plot <- inner_join(clouded_nroh, clouded_sroh, by = "Sample") %>% 
  ggplot(aes(x = SROH, y = NROH)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "NROH vs. SROH")

ggsave("froh_plot.png", froh_plot, width = 8, height = 6, dpi = 300)

# Create a boxplot of RoH lengths for each sample and save the plot to a file
roh_boxplot <- clouded_roh %>% 
  ggplot(aes(x = Sample, y = RoH_length)) +
  geom_boxplot() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "RoH Lengths per Sample")

ggsave("roh_boxplot.png", roh_boxplot, width = 8, height = 6, dpi = 300)
```
