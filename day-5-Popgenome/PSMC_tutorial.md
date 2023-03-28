# PSMC

Developed by Henrique V. Figueiró

Contact: henriquevf@gmail.com

The PSMC (Pairwise Sequential Markovian Coalescent) analysis is a widely used tool in genomics for inferring the demographic history of a population from genomic data. It is based on the coalescent theory, a framework for understanding how genetic diversity is shaped by the history of a population. By modeling the coalescent process in a pairwise fashion, PSMC can provide insights into changes in population size over time, allowing researchers to reconstruct the demographic history of a population.

Mapping files are located in the directory below. You can also use your own mapping files.

```bash
/pool/genomics/figueiroh/SMSC_2023_class/mapping
/pool/genomics/figueiroh/SMSC_2023_class/reference
```

Intermediate data are located in the following directory:

```bash
/pool/genomics/figueiroh/SMSC_2023_class/psmc
```

On your user create the following folder: 

```bash
mkdir psmc
```

### Consensus building

The initial step in the process is to generate a consensus sequence from the bam files. Although the command for this step is relatively old, it still functions effectively. However, it is important to ensure that the command is still operational and functioning correctly before proceeding. This step should take between 9 and 12 hours to run on Hydra with 8 Gb of RAM. 

```bash
module load bioinformatics/bcftools
bcftools mpileup -Ou -f <reference_genome> <bam_file> | bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > <output.fq.gz>
```

1. `bcftools mpileup -C50 -uf <reference_genome> <bam_file>`: This command generates a textual pileup format of the input BAM file (`<bam_file>`) using the given reference genome (`<reference_genome>`). The `C50` option applies a coefficient to adjust the base alignment quality, and the `u` flag outputs the results in the uncompressed BCF format, which is required for piping to `bcftools`. The `f` flag specifies the reference genome file.
2. `bcftools call -c`: This command performs variant calling on the input data received from the `bcftools mpileup` command (indicated by `` as input). The `c` option uses the consensus caller, which is suitable for calling a diploid consensus sequence.
3. `vcfutils.pl vcf2fq -d 10 -D 100`: This command is part of the `bcftools` package and converts the output from `bcftools call` (in VCF format) to a FastQ format. The `d 10` and `D 100` options set the minimum and maximum depth thresholds for filtering variants, respectively.
4. `gzip > <output.fq.gz>`: This part of the command compresses the final output using `gzip` and saves it as a `.fq.gz` file (`<output.fq.gz>`).

With the consensus sequence we generate the input file in order to run PSMC. For the sake of time you can copy the data from the original folder.

```bash
module load bioinformatics/psmc
fq2psmcfa -q20 <input.fq.gz> > <output.psmcfa>
```

This command will take a gzipped FastQ file `input.fq.gz` as input and create a PSMC input file `output.psmcfa` with bases that have quality scores below 20 masked (`-q20` ).The higher the quality threshold, the more stringent the filtering process, leading to more masked bases in the output.

### Run PSMC

```bash
module load bioinformatics/psmc
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o <output.psmc> <input.psmcfa>
```

1. `psmc`: This is the main command to run the PSMC tool.
2. `N25`: This flag sets the effective population size (`N`) to 25. The effective population size is a measure of the genetic diversity in a population and is used to calculate the time to the most recent common ancestor of the population.
3. `t15`: This flag sets the scaled mutation rate per generation (`t`) to 15. The scaled mutation rate is the product of the mutation rate per base pair per generation and the effective population size.
4. `r5`: This flag sets the scaled recombination rate per generation (`r`) to 5. The scaled recombination rate is the product of the recombination rate per base pair per generation and the effective population size.
5. `p "4+25*2+4+6"`: This flag sets the time intervals (`p`) for the PSMC model. The specified pattern, "4+25*2+4+6", means that there are 4 intervals of equal size at the start, followed by 25 intervals with twice the size of the previous intervals, and then 4 more intervals of equal size, and finally 6 more intervals of increasing size. This allows the model to have higher time resolution near the present and lower resolution in the more distant past.
6. `o <output.psmc>`: This flag specifies the output file name for the PSMC results. Replace `<output.psmc>` with the desired output file name.
7. `<input.psmcfa>`: This is the input file in PSMCFA format, which contains the sequence data to be analyzed. Replace `<input.psmcfa>` with the name of the input file.

### Simulation

You can also create the command to run simulations based on your results. We are not using this for the tutorial, but you have the option if you need.

```bash
module load bioinformatics/psmc
psmc2history.pl <input.psmc> | history2ms.pl > <output_ms-cmd.sh>
```

This command consists of two parts that involve the use of Perl scripts from the PSMC utilities package. The scripts convert the PSMC output file into an ms command, which can then be used to simulate genetic data using the ms software.

1. `utils/psmc2history.pl <input.psmc>`: This script takes the PSMC output file (`<input.psmc>`) and converts it into a simple history format. This format represents the inferred demographic history of the population.
2. `utils/history2ms.pl > <output_ms-cmd.sh>`: This script takes the output from the previous script (the simple history format) and generates an ms command. The ms command is then saved in an output shell script file (`<output_ms-cmd.sh>`).

The resulting shell script file (`<output_ms-cmd.sh>`) contains the ms command that you can use to simulate genetic data based on the demographic history inferred by the PSMC analysis. The ms software is a widely used tool for simulating coalescent-based gene genealogies and genetic data under various demographic scenarios.

### Plotting the results

Finally you can plot the results using the provided script. After that you should download the output to your machine. 

```bash
module load bioinformatics/psmc
psmc_plot.pl -g 7 -u 1e-8 -p <input_name> </path/to/file/psmc/input.psmc>
```

1. `-g 7`: This option specifies the maximum generation time in years (tmax) that will be used for the analysis. In this case, the value is 7 years.
2. `-u 1e-8`: This option sets the mutation rate per site per generation (u) used in the analysis. The value is set to 1x10^-8.
3. `-p`: This option is used to save the generated plot in the pdf format. 
4. `input_name`: This option specifies a prefix to be used for the output files generated by the script. In this case, the prefix is "NN114296".
5. `input.psmc`: This argument specifies the input file for the script. In this case, the input file is "NN114296.psmc".

The process outlined above is relatively straightforward, with the main point of concern being the selection of appropriate values for the generation time and mutation rate of the species being studied. However, obtaining accurate information on these parameters can be challenging, particularly for taxonomic groups where such data is scarce or non-existent. In such cases, researchers may need to rely on estimates or approximations to proceed with their analysis.

Now, try to run with a generation time of four years. What is the difference you can observe in the plot?

```bash
module load bioinformatics/psmc
psmc_plot.pl -g 4 -u 1e-8 -p NN114296_g4 /path/to/file/psmc/NN114296.psmc
```
### Bootstrapping

In this tutorial, we will walk you through the process of running the Pairwise Sequentially Markovian Coalescent (PSMC) model with bootstrap replicates to infer the population size history from a single diploid individual's genome sequence. We assume that you have already installed PSMC and its dependencies.

**Step 1: Prepare your input data**

Before you start, make sure you have a FASTQ file containing the diploid genome sequence of the individual you want to analyze. This file should be compressed using gzip and have a .fq.gz extension.

**Step 2: Convert FASTQ to PSMCFA format**

Convert the FASTQ file to PSMC input format (.psmcfa) using the following command:

```bash
fq2psmcfa -q20 /path/to/input/FASTQ_file.fq.gz > /path/to/output/psmcfa_file.psmcfa
```

Replace the paths and file names with the appropriate values for your analysis.

**Step 3: Split the PSMCFA file**

Split the .psmcfa file into multiple smaller files for PSMC program input:

```bash
splitfa /path/to/input/psmcfa_file.psmcfa > /path/to/output/split_psmcfa_file.split.psmcfa
```

**Step 4: Run the PSMC program**

Run the PSMC program with specified parameters on the .psmcfa file to infer population size history:

```bash
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /path/to/output/psmc_output_file.psmc /path/to/input/psmcfa_file.psmcfa
```

**Step 5: Run bootstrap replicates**

To perform bootstrap replicates, run the following command:

```bash
seq 100 | xargs -P 4 -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-{}.psmc /path/to/input/split_psmcfa_file.split.psmcfa | sh
```

This command will run 100 bootstrap replicates of the PSMC analysis in parallel, using **4 processes** at a time. The output files will be named round-1.psmc, round-2.psmc, ..., round-100.psmc.

**Step 6: Combine the PSMC output with bootstrap replicates**

Combine the original PSMC output file with the 100 bootstrap replicate output files into a single combined file:

```bash
cat /path/to/input/psmc_output_file.psmc round-*.psmc > /path/to/output/combined_psmc_output_file.combined.psmc
```

**Step 7: Generate a plot of the inferred population size history**

Generate a plot of the inferred population size history using the combined PSMC output file:

```bash
psmc_plot.pl -g 4 -u 1e-8 -X 1000000 /path/to/input/combined_psmc_output_file.combined /path/to/output/psmc_plot_file
```

This command will create a PDF file containing the PSMC plot, which shows the inferred population size history.

**Note on using an HPC:**

Here is the consolidated commands to be used on a job file. Remember to give at least 4 threads, or edit the `xargs -P 4` to use more or less threads. Some steps can take a long time, particularly the bootstrap, so use processing times longer than 24 hours.

```bash
fq2psmcfa -q20 /path/to/input/FASTQ_file.fq.gz > /path/to/output/psmcfa_file.psmcfa

splitfa /path/to/input/psmcfa_file.psmcfa > /path/to/output/split_psmcfa_file.split.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /path/to/output/psmc_output_file.psmc /path/to/input/psmcfa_file.psmcfa

seq 100 | xargs -P 4 -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-{}.psmc /path/to/input/split_psmcfa_file.split.psmcfa | sh

cat /path/to/input/psmc_output_file.psmc round-*.psmc > /path/to/output/combined_psmc_output_file.combined.psmc

psmc_plot.pl -g 4 -u 1e-8 -X 1000000 /path/to/input/combined_psmc_output_file.combined /path/to/output/psmc_plot_file
