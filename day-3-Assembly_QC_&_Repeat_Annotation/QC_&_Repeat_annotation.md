# Genome assembly QC and Genome annotation - SMSC (Day 3)

<!-- TOC depthFrom:2 -->

 * [Folder structure](#Folder-structure)
 * [Input file](#Input-file)
 * [Run BUSCO](#Run-BUSCO)
 * [Run Blobtoolskit2](#Run-Blobtooskit2)
 * [Masking and annotating repetitive elements with Repeatmodeler and RepeatMasker](#Masking-and-annotating-repetitive-elements-with-Repeatmodeler-and-RepeatMasker)


<!-- /TOC -->

### Folder structure

Let's create our folder structure for this workshop. It's easier to troubleshoot any issues if we are all working within the same framework. First, we will create a folder called `genome_qc_annotation` in your `/scratch/genomics/your_username` folder. 

```
cd /scratch/genomics/your_username/
mkdir genome_qc_annotation
```

Next, we will change directories to the genome_annot folder and we will create a several folders that we will use today and tomorrow. Here's the list of folders:

- Busco
- Bloobtools
- Repeat_annotation
- Gemoma

<details><summary>SOLUTION</summary>
<p>

```
mkdir assembly busco blobtools repeat_annotation gemoma

```
If you type the command `tree`, this is what you should see:

```
|__ assembly
|__ busco  
|__ blobtools
|__ gemoma
|__ repeat_annotation
 ```

</p>
</details>


If we follow this folder structure, we will have all the results organized by software, which will facilitate finding everything later. 

**Submitting jobs**: With this folder structure, we will save all the job files with each program folder and the job files will be submitted from there.   

### Input file
For this session, we will use the cloud leopard assembly. You can find this file here: `/data/genomics/workshops/smsc_2023/mNeoNeb1.pri.cur.20220520.fasta.gz`

Copy this file to your assembly folder and gunzip the file.

```
cp /data/genomics/workshops/smsc_2023/mNeoNeb1.pri.cur.20220520.fasta.gz /scratch/genomics/your_username/assembly/
cd /scratch/genomics/your_username/assembly/
gunzip mNeoNeb1.pri.cur.20220520.fasta.gz
```

**If you want to run thing quckly you can run the programs by Extracting some scaffolds**
<details><summary>SOLUTION</summary>
<p>

To generate this file, we used `bioawk` and `samtools` to extract the sequences from the original assembly:

`module load bio/bioawk`
`module load bio/samtools`

Create a list with the 10 largest sequences. The number of sequences is determined by the number following the `head` command in the end. 

`cat assembly.fasta | bioawk -c fastx '{ print length($seq), $name }' | sort -k1,1rn | head -10 > 10largest.list`

Use `samtools` to extract the list of sequences from the original assembly:

`xargs samtools faidx assembly.fasta < 10largest.list > assembly_10largest.fasta` 

</p>
</details>

### Run BUSCO

BUSCO (Simão et al. 2015; Waterhouse et al. 2017) assesses completeness by searching the genome for a selected set of single copy orthologous genes. There are several databases that can be used with BUSCO and they can be downloaded from here: [https://buscos.ezlab.org](https://buscos.ezlab.org). 


#### Job file: busco_cloud_leopard.job
- Queue: medium
- PE: multi-thread
- Number of CPUs: 10
- Memory: 6G (6G per CPU, 60G total)
- Module: `module load bio/busco/5.4.3`
- Commands:

```
busco  -o clouded_leopard -i ../assembly/mNeoNeb1.pri.cur.20220520.fasta -l mammalia_odb10 -c $NSLOTS -m genome
```

##### Explanation:
```
-o: name of the output folder and files
-i: input file (FASTA)
-l: name of the database of BUSCOs (This will automatically connected and dowloand the database from the BUSCO website).
-c: number of CPUs
-m: mode (options are genome, transcriptome, proteins)

```

**Note about Databases:**

If you do not have internet connection on the node where running the software you can download the database and run the program offline. For instance, to download the Mammalia database you can use the command `wget` and extract it.It is important to dowlod and untar the folder on your busco folder. Let's `cd` to the directory `busco` first.

	wget https://busco-data.ezlab.org/v5/data/lineages/mammalia_odb10.2021-02-19.tar.gz
	tar -zxf aves_odb9.tar.gz

The command to run busco will have to change to: 

```
busco  -o clouded_leopard -i ../assembly/mNeoNeb1.pri.cur.20220520.fasta -l mammalia_odb10 -c $NSLOTS -m genome --offline --download_path /path/to/datasets
```

### Run Bloobtools2

BlobTools2 is a command line tool designed for interactive quality assessment of genome assemblies and contaminant detection and filtering.

#### Preparing files for blobtools2

First, you need to blast your assembly to know nt databases. For this we will used blastn program. 


#### Job file: blast_clouded_leopard.job
- Queue: medium
- PE: multi-thread
- Number of CPUs: 10
- Memory: 6G (6G per CPU, 60G total)
- Module: `module load bio/blast/2.13.0`
- Commands:

```
blastn -db /data/genomics/db/ncbi/db/latest_v4/nt/nt -query ../assembly/mNeoNeb1.pri.cur.20220520.fasta -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 20 -max_hsps 1 -evalue 1e-20 -num_threads $NSLOTS -out clouded_leopard_blast.out
```

##### Explanation:
```
-db: ncbi nucleotide database
-query: input file (FASTA)
-outfmt: format of the output file (important to for blobtools) 
-max_target_seqs: Number of aligned sequences to keep.
-max_hsps: Maximum number of HSPs (alignments) to keep for any single query-subject pair.
-num_threads: number of CPUs
-out: name of the output file
```

Second, you need to map raw reads to the genome assembly. We will use minimap2 for this. Minimap2 is a versatile sequence alignment program that aligns DNA or mRNA sequences against reference database. Typical use cases include: (1) mapping PacBio or Oxford Nanopore genomic reads to a reference genome; or (2) aligning Illumina single- or paired-end reads to a reference genome. After mapping the reads we need to convert the output file SAM into a BAM fil and sort this file. For this we will the program samtools. Samtools is a suite of programs for interacting with high-throughput sequencing data. 

#### Job file: minimap_clouded_leopard.job
- Queue: medium
- PE: multi-thread
- Number of CPUs: 10
- Memory: 6G (6G per CPU, 60G total)
- Module: 
```
  module load bio/samtools/1.9 
  module load bio/tools/conda     
  start-conda  
  conda activate minimap2
```
- Commands:

```
minimap2 -ax map-hifi -t 20 ../assembly/mNeoNeb1.pri.cur.20220520.fasta  /path/to_each/hifi.fastq | samtools sort -@20 -O BAM -o cloud_leopard_sorted.bam -
```

##### Explanation:

```
minimap2
-ax: preset configuration to map hifi reads to genomes.
-t: number of threads to use.
Samtools
sort: sort command
-@: number of threads to use.
-O: output format.
-o: name of the outputformat
```

#### Creating blobtools2 data base

Now that we have the blast and mapping results we can create the BlobTools2 database. The minimum requirement to create a new database with BlobTools2 is an assembly FASTA file. This runs very fast so we do can use an interactive node.

```
qrsh -pe mthread 3
module load bio/blobtools/2.6.3 
blobtools create --fasta ../assembly/mNeoNeb1.pri.cur.20220520.fasta clouded_leopard_blobt_nb
```

##### Explanation:

```
create: command to creat the database with blobtools2.
--fasta: path and name of the assembly fasta file.
clouded_leopard_blobt_nb: name of the database.
```
#### Adding data to a database

Once you have a BlobDir database, other data can be added by parsing analysis output files into one or more fields using the `blobtools add` command. Preset parsers are available for a range of analysis types, such as BLAST or Diamond hits with taxonomic assignments for scaffolds/contigs; read mapping files provide base and read coverage and BUSCO results that show completeness metrics for the genome assembly.
Again, since this is not computationally expensive we can use an interactive node.


```
qrsh -pe mthread 3
module load bio/blobtools/2.6.3 
blobtools add --threads 3 --hits clouded_leopard_blast.out --taxrule bestsumorder --taxdump /pool/genomics/ariasc/SMSC_2023/blobtools/taxdump clouded_leopard_blobt_nb 
blobtools add --threads 3 --cov cloud_leopard_sorted.bam clouded_leopard_blobt_nb # this can take a several minutes 
blobtools add --threads 3 --busco full_table.tsv clouded_leopard_blobt_nb
```

#### Create interactive html page with the blobtools results

After you finish creating and adding data to the database in order to visulize the results you need to install blobtools2 on your personal machine and  download the database folder. First lets tar zip the folder with the command ```tar -czvf name-of-archive.tar.gz /path/to/directory-or-file``` and  download the folder using the ffsend (load ```bio/ffsend``` module). 

Now let's install blobtools on your machine. For this we will use conda...

Please download your folder from the ffsend link, and move the file to a new folder on your machine. After downlodaing the folder you need to untar the folder

<details><summary>SOLUTION</summary>
<p>
Please download your folder from the ffsend link, and move the file to a new folder on your machine. After downlodaing the folder you need to untar the folder:

```
mkdir blobtools_results
cd blobtools_results
#move or ffsend folder to this folder
tar -xzvf archive.tar.gz
```
</p>
</details>

After installing blobtools2 make sure that you have activate your conda environment and from the program folder run the following command on your database folder.

```
conda activate btk_env
./blobtools view --remote path/to/clouded_leopard_nb
```
The cool thing about this is that you can interact with the results and visualize them in several ways. It also give you nice images that can be use on your publication.

* How long is the longest contig and scaffold?
* What is the contig and scaffold N50?
* Are there any contamination?
* if yes, what taxa are contaminant of your assembly?


### Masking and annotating repetitive elements with Repeatmodeler and RepeatMasker

Repeatmodeler is a repeat-identifying software that can provide a list of repeat family sequences to mask repeats in a genome with RepeatMasker. Repeatmasker is a program that screens DNA sequences for interspersed repeats and low complexity DNA sequences. the output of the program is a detailed annotation of the repeats that are present in the query sequence as well as a modified version of the query sequence in which all the annotated repeats have been masked. 
Things to consider with Repeatmodeler software is that it can take a long time with large genomes (>1Gb==>96hrs on a 16 cpu node). You also need to set the correct parameters in repeatmodeler so that you get repeats that are not only grouped by family, but are also annotated.

Repeatmodeler http://www.repeatmasker.org/RepeatModeler/  
RepeatMasker http://www.repeatmasker.org/RMDownload.html

The first step to run Repeatmodeler is that you need to build a Database. The Build Database step is quick (several seconds at most).

#### Job file: repeatmodeler_database.job
- Queue: medium
- PE: multi-thread
- Number of CPUs: 1
- Memory: 10G
- Module: `module load bio/repeatmodeler`
- Commands:

```
BuildDatabase -name cloud_leopard ../assembly/mNeoNeb1.pri.cur.20220520.fasta
# usage:BuildDatabase -name {database_name} {genome_file-in_fasta_format}
```
##### Explanation:
```
-name: name to be given to the database
genome file in fasta format
```

##### Output files:
- Database folder with the structure to populate and run repeatmodeler program.

The second step is two actually run RepeatModeler. Again this can take several days depending on the genome size.

#### Job file: repeatmodeler_cloud_leopard.job
- Queue: high
- PE: multi-thread
- Number of CPUs: 36
- Memory: 10G
- Module: `module load bio/repeatmodeler`
- Commands:

```
# Usage: RepeatModeler -database {database_name} -pa {number of cores} -LTRStruct > out.log
RepeatModeler -database cloud_leopard -pa 36 -LTRStruct -engine ncbi > out.log
```

##### Explanation:
```
-database: The prefix name of a XDF formatted sequence database containing the genomic sequence to use when building repeat models.
-pa: number of cpus
-engine: The name of the search engine we are using. I.e abblast/wublast or ncbi (rmblast version).
-LTRStruct: enables the optional LTR structural finder.
```

##### Output files:
- consensi.fa.classified: complete database to be in RepeatMasker.

The last step to get a repeat annotation is to run ReapeatMasker. 


#### Job file: repeatmasker_cloud_leopard.job
- Queue: high
- PE: multi-thread
- Number of CPUs: 36
- Memory: 10G
- Module: `module load bioi/repeatmodeler`
- Commands:

```
# usage: RepeatMasker -pa 36 -gff -lib {consensi_classified} -dir {dir_name} {genome_in_fasta}

RepeatMasker -pa $NSLOTS -xsmall -gff -lib consensi.fa.classified -dir ../repeatmasker/cloud_leopard_RM ../assembly/cloud_leopard_10largest.fasta 
```
##### Explanation:
```
-lib: repeatmodeler repbase database to be search
-pa: number of cpus
-xsmall: softmasking (instead of hardmasking with N)
-dir ../repmasker: writes the output to the directory repmasker
-gff: output format of the annotated repeats
```

##### Output files:
- cloud_leopard_10largest.fasta.tbl: summary information about the repetitive elements
- cloud_leopard_10largest.fasta.masked: masked assembly (in our case, softmasked)
- cloud_leopard_10largest.fasta.out: detailed information about the repetitive elements, including coordinates, repeat type and size.
