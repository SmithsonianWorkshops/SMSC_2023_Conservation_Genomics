# Variant Calling - SMSC (Day 4)

<!-- TOC depthFrom:2 -->

 * [Folder structure](#Folder-structure)
 * [Runing Minimap2](#Runing-Minimap2)
 * [Input files](#Input-file)
 * [Map to the reference genome](#Map-to-the-reference-genome)
 * [The Genome Analysis Toolkit (GATK)](#The-Genome-Analysis-Toolkit-(GATK))
 


<!-- /TOC -->

### Folder structure

Let's create a new folder  called `pop_genomics` in your  `/scratch/genomics/your_username`. 
our folder structure for this workshop. It's easier to troubleshoot any issues if we are all working within the same framework. First, we will create a folder called `genome_qc_annotation` in your `/scratch/genomics/your_username` folder. 

```
cd /scratch/genomics/your_username/
mkdir pop_genomics
```

Next, we will change directories to the pop_genommics folder and we will create a several folders that we will use today and tomorrow. Here's the list of folders:

- variant_calling
- pop_structure
- PSMC

<details><summary>SOLUTION</summary>
<p>

```
mkdir variant_calling pop_structure PSMC 

```
If you type the command `tree`, this is what you should see:

```
|__ PSMC  
|__ pop_genomics
|__ variant_calling

 ```

</p>
</details>


If we follow this folder structure, we will have all the results organized by software or analysis, which will facilitate finding everything later. 

**Submitting jobs**: With this folder structure, we will save all the job files with each program folder and the job files will be submitted from there.   

### Input files
For this session, we will use the reference cloud leopard assembly. You can find this file here: `/data/genomics/workshops/smsc_2023/clouded_leopard_pacbio/mNeoNeb1.pri.cur.20220520.fasta`

We will also need the illumina data (fastq) for our 7 Clouded Leopard individuals this will represent our population. Remember that the data is locted here: `/data/genomics/workshops/smsc_2023/clouded_leopard_illumina`


### Map to the reference genome

The firs step to call for variants is to map all our individuals to our reference genome. Similar to yesterday, we will use minimap2 and samtools to map and convert from SAM file to BAM file. SAM files are tab-delimited text files that are usually very large in size. To save space and make it easy for software to hand large aliments, SAM is usually converted to a binary format BAM. 


#### Job file: minimap_pop_clouded_leopard.job
- Queue: medium
- PE: multi-thread
- Number of CPUs: 10
- Memory: 6G (6G per CPU, 60G total)
- Module: 
```
module load bio/minimap/
```
- Commands:

```
minimap2 -ax sr -t 20 /data/genomics/workshops/smsc_2023/clouded_leopard_pacbio/mNeoNeb1.pri.cur.20220520.fasta /data/genomics/workshops/smsc_2023/clouded_leopard_illumina/NN114296_1.fastq.gz /data/genomics/workshops/smsc_2023/clouded_leopard_illumina/NN114296_2.fastq.gz | samtools sort -@20 -O BAM -o NN114296_cloud_leopard_sorted.bam -
```

##### Explanation:

```
minimap2
-ax: preset configuration to map illumina short reads to genomes.
-t: number of threads to use.
Samtools
view: convert comand 
-b  output format BAM
sort: sort command
-@: number of threads to use.
-O: output format.
-o: name of the outputformat
```

Note! this step is quite time-consuming. I have run this for you before so we do not have to wait 3-5 hours to complete each of seven mappings. I have created a folder from where you can call all the BAM files `/pool/genomics/ariasc/SMSC_2023/mapping/`


### The Genome Analysis Toolkit (GATK)

The GATK (Genome Analysis Toolkit) is one of the most used programs for genotype calling in sequencing data in model and non model organisms. However, the GATK was designed to analyze human genetic data and all its pipelines are optimized for this purpose. Variant calling with GATK currently requieres 4 steps: Mark duplicates, 




####  Mark Duplicates with picard-tools

Potential PCR duplicates need to be marked. Marking duplicates make sense even if you used a PCR-free library preparation procedure because reads identified as duplicates are not removed and can be included in the subsequent analyses if needed.

#### Job file: markduplicates_1.job
- Queue: medium
- PE: multi-thread
- Number of CPUs: 5
- Memory: 10G (10G per CPU, 50G total)
- Module: 
```
  module load bio/gatk/4.1.3.0 
  
```
- Commands:

```
rungatk MarkDuplicates -I=NN114296_cloud_leopard_sorted.bam -O=NN114296_sorted_reads_duplMarked.bam -M=sorted_reads_duplMarked.metrics
```

##### Explanation:

```
-I: Input file (sorted_bam)
-O: Output file (duplmarked bam)
-M: Output for metrics file
```

#### add or replace group arguments on duplicate marked BAM file

The GATK requires read group information in BAM files. It is used to differentiate samples and to detect artifacts associated with sequencing techniques. 

#### Job file: add_group_arguments_1.job
- Queue: medium
- PE: multi-thread
- Number of CPUs: 5
- Memory: 10G (10G per CPU, 50G total)
- Module: 
```
  module load bio/gatk/4.1.3.0 
```
- Commands:

```
rungatk AddOrReplaceReadGroups -I=NN114296_sorted_reads_duplMarked.bam -O=NN114296_sorted_reads_duplMarked_SM.bam -ID=1 -LB=lib1 -PL=ILLUMINA -PU=unit1 -SM=NN114296

```

##### Explanation:

```
-I: Input file (sorted_bam)
-O: Output file (duplmarked bam)
-ID: Read Group ID, lane number
-LB: librari ID
-PL: Plataform (i.e. illumina)
-PU: run Barcode 
-SM: Sample Name
```

#### create index for your reference genome

first we will copy our reference to our variant calling folder and them we will run CreateSequenceDictionary.

```
cp /data/genomics/workshops/smsc_2023/clouded_leopard_pacbio/mNeoNeb1.pri.cur.20220520.fasta .
rungatk CreateSequenceDictionary -R=mNeoNeb1.pri.cur.20220520.fasta -O=mNeoNeb1.pri.cur.20220520.dict
```
##### Explanation:

```
-R: Reference Genome (Fasta)
-O: Output reference index dictionary

```

#### Base quality recalibration

The GATK best practices recommend performing Base Quality Score Recalibration. This procedure detects systematic errors in your data by comparing it to the reference training data set. The problem with non-mode organisms is that there is usually no large high confidence genotype database that is required for training. One of the solutions to this problem is to create a custom training database from the same data you analyze but using very stringent filtering criteria and re-calibrate the data. Others have addressed this question at the GATK forum and although the approach was supported by the GATK team, in practice it did not works well. It simply shifts the distribution to lower scores and in many cases biases result. Thus people suggest to not performed this is step on non-model organism.

#### Variant Calling

The next step in the GATK best practices workflow is to proceed with the variant calling. There are a couple of workflows to call variants using GATK4. Here we will follow the Genomic Variant Call Format (GVCF) workflow which is more suited for scalable variant calling i.e. allows
incremental addition of samples for joint genotyping. 

This also requires several steps:

a) HaplotypeCaller 
HaplotypeCaller is the focal tool within GATK4 to simultaneously call germline SNVs and small Indels using local de novo assembly of haplotype regions. Briefly, the HaplotypeCaller works by: 1. Identify regions with evidence of variations. 2. Re-asssemble the active regions.  3. Determine likelihoods of the haplotypes given the data. 4. assigns samples genotypes. For each potential variant site, the program applies Baye's rule using the likelihood of the allele given the data.

#### Job file: HaplotypeCaller_1.job
- Queue: medium
- PE: multi-thread
- Number of CPUs: 5
- Memory: 10G (10G per CPU, 50G total)
- Module: 
```
  module load bio/gatk/4.1.3.0 
```
- Commands:

```
rungatk HaplotypeCaller -I NN114296_sorted_reads_duplMarked_SM.bam -R /data/genomics/workshops/smsc_2023/clouded_leopard_pacbio/mNeoNeb1.pri.cur.20220520.fasta -ERC GVCF -O NN114296.g.vcf.gz

```

##### Explanation:

```
-I: Input file (sorted_duplicatmak_SM.bam)
-R: Reference Genome (Fasta)
-O: Output file (g.vcf file)
-ERC: preset parameters for likelihood calculations and is more efficient and recommended for large number of samples and therefore more scalable.
```

b) Apply CombineGVCFs
The CombineGVCFs tool is applied to combine multiple single sample GVCF files, merging them into a single multi-sample GVCF file.
We have pre-processed two additional samples (NN and NN) up to the HaplotypeCaller step (above).

#### Job file: CombineGVCFs.job
- Queue: medium
- PE: multi-thread
- Number of CPUs: 5
- Memory: 10G (10G per CPU, 50G total)
- Module: 
```
  module load bio/gatk/4.1.3.0 
```
- Commands:

```
rungatk CombineGVCFs -R /data/genomics/workshops/smsc_2023/clouded_leopard_pacbio/mNeoNeb1.pri.cur.20220520.fasta \
-V NN114296.g.vcf.gz \
-V NN114296.g.vcf.gz \
-V NN114296.g.vcf.gz \
-O clouded_leopard_cohort.g.vcf.gz
```

##### Explanation:

```
-V: each of the individual g.vcf files. One -V for each.
-R: Reference Genome (Fasta)
-O: Output file with all the individuals combine (i.e. X_cohort.g.vcf file)
```

Now that we have a combine GVCF file, we are ready to perform genotyping.

c) Apply GenotypeGVCFs
This tool is designed to perform joint genotyping on a single input, which may contain one or many samples. In any case, the input samples must possess genotype likelihoods produced by HaplotypeCaller.

#### Job file: GenotypeGVCFs.job
- Queue: medium
- PE: multi-thread
- Number of CPUs: 5
- Memory: 10G (10G per CPU, 50G total)
- Module: 
```
  module load bio/gatk/4.1.3.0 
```
- Commands:

```
rungatk GenotypeGVCFs \
-R /data/genomics/workshops/smsc_2023/clouded_leopard_pacbio/mNeoNeb1.pri.cur.20220520.fasta \
-V clouded_leopard_cohort.g.vcf.gz \
-O clouded_leopard_output.vcf.gz

```

##### Explanation:

```
-V: cohort_g.vcf files.
-R: Reference Genome (Fasta)
-O: Output file with the genotype in vcf format (i.e. X_output.vcf file)
```

#### SNP filtering

The VariantFiltration tools is designed for hard-filtering variant calls based on custom quality criteria such as sequencing depth, mapping quality etc.

#### Job file: GenotypeGVCFs.job
- Queue: medium
- PE: multi-thread
- Number of CPUs: 5
- Memory: 10G (10G per CPU, 50G total)
- Module: 
```
  module load bio/gatk/4.1.3.0 
```
- Commands:

```
rungatk VariantFiltration \
-R /data/genomics/workshops/smsc_2023/clouded_leopard_pacbio/mNeoNeb1.pri.cur.20220520.fasta \
-V clouded_leopard_output.vcf.gz \
-O cloduded_leopard_varfilter.vcf \
--filter-name "Low_depth10" \
--filter-expression "DP < 10"
```

##### Explanation:

```
-V: genotyped output vcf files
-R: Reference Genome (Fasta)
-O: Output file with the variant filter results in vcf format (i.e. X_variantfilter.vcf file)
--filter-name: a given filer name that will be printed in the vcf file for your information
-filter-expression: hard filter that you want to perform.
```

Yeah!!!! Now we have a final VCF clean and ready for our population genomic analyses.


