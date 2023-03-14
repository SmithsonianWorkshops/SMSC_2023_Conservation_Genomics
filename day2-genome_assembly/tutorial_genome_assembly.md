### 1. Running FastQC on raw data
* FastQC is a program that can quickly scan your raw data to help figure out if there are adapters or low quality reads present. Create a job file to run FastQC on one of the fastq files here: ```/scratch/genomics/dikowr/cloud_leopard_raw_data```
	+ **module**: ```bioinformatics/fastqc```
	+ **command**: ```fastqc <FILE.fastq>```
	+ after your job finishes, find the results and download some of the images, e.g. ```per_base_quality.png``` to your local machine using ffsend (load ```bioinformatics/ffsend``` module) and then the command ```ffsend upload <FILE>```.


### 2. Trimming adapters with TrimGalore! 
* PacBio data will be error-corrected solely by the assembler, but Illumina data trimming and thinning are common.
* Most assemblers these days don't want you to trim/thin for quality before assembling, but trimming is important for downstream applications. TrimGalore will auto-detect what adapters are present and remove very low quality reads (quality score <20) by default.  
* Create a job file to trim adapters and very low quality reads for the Illumina data here: ```/scratch/genomics/dikowr/cloud_leopard_data/illumina_raw```
	+ **command**: ```trim_galore --paired --retain_unpaired <FILE_1.fastq> <FILE_2.fastq>```  
	+ **module**: ```bio/trim_galore```
	+ You can then run FastQC again to see if anything has changed.

### 3. Run Genomescope

* Genomescope can be used to estimate genome size from short read data: 
	+ [Genomescope](http://qb.cshl.edu/genomescope/) 

* To run Genomescope, first you need to generate a Jellyfish histogram.

* You'll need two job files for Jellyfish, one to count the kmers and the second to generate a histogram to give to Genomescope: 
* Here is a copy of the Cloud Leopard Illumina data: ```/scratch/genomics/dikowr/cloud_leopardd_raw_data/illumina_trimmed```
	+ Hint: don't copy these data to your own space.

* First job file: kmer count:
	+ Module: ```bio/jellyfish```
	+ Commands: ```jellyfish count -C -m 21 -t $NSLOTS -s 800000000 *.fastq -o reads.jf```
	+ ```-m``` = kmer length  
	+ ```-s``` = RAM  
	+ ```-C``` = "canonical kmers" don't change this 
	+ Hint: this job needs to run on the high memory queue. 

* This will take a while, so we can move on and then come back to it when it finishes.

* Second job file: histogram:
	+ Module: ```bio/jellyfish```
	+ Commands: ```jellyfish histo -t $NSLOTS reads.jf > reads.histo```

* Download the histogram to your computer (e.g. using ffsend again), and put it in the Genomescope webservice: [Genomescope](http://qb.cshl.edu/genomescope/)

* let's run the analysis the same analysis with the HiFi PacBio data.

### 4. Hifiasm Assembly

Hifiasm is a fast haplotype-resolved de novo assembler for PacBio HiFi reads. it can produce partially phased assemblies of great quality. Hifiasm can use trio (Parental) short reads data or Hi-C data to produce haplotype-resolved assemblies.

#### Assembly using only HiFi data

* Job file: hifi_only.job

  + **PE:** multi-thread
  + **Queue:** Long, himem 
  + **Number of CPUs:** 32
  + **RAMMemory:** 10G (10G per CPU, 300G total)
  + **Module:** `module load bio/hifiasm`
  + **Command:**
```hifiasm -o cloud_leopard_only.asm -t 32 cloud_leopard_hifi.fq.gz```

##### Comand explanation:
```
-o: name and path of the output file in asm format
-t: sets the number of CPUs in use

cloud_leopard_hifi.fq.gz Input reads. Input sequences should be FASTA 
or FASTQ format, uncompressed or compressed with gzip (.gz). The quality scores of reads 
in FASTQ are ignored by hifiasm. Hifiasm outputs assemblies in `GFA <https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md>`_ format.
```
 
* This job should complete in a few hours

#### Assembly using both HiFi data and Hi-C

* Job file: hifi_Hi_C.job

  + **PE:** multi-thread
  + **Queue:** Long, himem 
  + **Number of CPUs:** 32
  + **RAMMemory:** 10G (10G per CPU, 300G total)
  + **Module:** `module load bio/hifiasm`
  + **Command:**
  ```hifiasm -o cloud_leopard_hic.asm -t32 --h1 hi_c_read1.fq.gz --h2 hi_c_read2.fq.gz cloud_leopard_hifi.fq.gz```

* This job should complete in a few hours
* In this mode, each contig represents a haplotig, which means that comes from one parental haplotype only.
* Hifiasm doen not performed scaffolding. For this run scaffolder such as SALSA.

##### Command explanation:
```
-o: name and path of the output file in asm format
-t: sets the number of CPUs in use
--h1: path to read 1 of Hi_C data
--h2: path to read 2 of Hi_C data
cloud_leopard_hifi.fq.gz Input reads. Input sequences should be FASTA 
or FASTQ format, uncompressed or compressed with gzip (.gz). The quality scores of reads 
in FASTQ are ignored by hifiasm. Hifiasm outputs assemblies in `GFA <https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md>`_ format.
```

### 5. Run the fasta metadata parser to get statistics about both the primary assembly and haplotype-resolved assemblies.
* We use a python script to grab some statistics from assembly files. These are the stats it produces:  

Total number of base pairs:    
Total number of contigs:   
N90:  
N80:  
N70:  
N60:  
N50:  
L90:  
L80:  
L70:  
L60:  
L50:  
GC content:  
Median contig size:  
Mean contig size:  
Longest contig is:  
Shortest contig is: 

* We have put a finished assembly here: ```/data/genomics/workshops/STRI_genomics/cloud_leopard_final_assembly.fasta```
	+ Module: ```bio/assembly_stats```
	+ Commands: ```assembly_stats <ASSEMBLY> > assembly_stats.out```

* How long is the longest contig and scaffold?
* What is the contig and scaffold N50?
* Are there differences between all the assemblies?
