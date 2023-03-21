# Genome assembly QC and Genome annotation - SMSC (Day 3)

<!-- TOC depthFrom:2 -->

 * [Folder structure](#Folder-structure)
 * [Input file](#Input-file)
 * [Run BUSCO](#Run-BUSCO)
 * [Run Blobtoolskit](#Run-Blobtoolskit)
 * [Run Blobtoolskit2](#Run-Blobtoolskit2)
 * [Masking and annotating repetitive elements with Repeatmodeler and RepeatMasker](#Masking-and-annotating-repetitive-elements-with-Repeatmodeler-and-RepeatMasker)
 * [Runing Gemoma](#Runing-Gemoma)


<!-- /TOC -->

### Folder structure

Let's create our folder structure for this workshop. It's easier to troubleshoot any issues if we are all working within the same framework. First, we will create a folder called `genome_qc_annotation` in your `/scratch/genomics/your_username` folder. 

```
cd /scratch/genomics/your_username/
mkdir genome_qc_annotation
```

Next, we will change directories to the genome_annot folder and we will create a several folders that we will use today and tomorrow. Here's the list of folders:

- busco
- bloobtools
- repeat_annotation
- gemoma

<details><summary>SOLUTION</summary>
<p>

```
mkdir busco blobtools repeat_annotation gemoma

```
If you type the command `tree`, this is what you should see:

```
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


**If you want to run thing quickly you can run the programs by Extracting some scaffolds**
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
busco  -o clouded_leopard -i path/to_assembly/mNeoNeb1.pri.cur.20220520.fasta -l mammalia_odb10 -c $NSLOTS -m genome
```

##### Explanation:
```
-o: name of the output folder and files
-i: input file (FASTA)
-l: name of the database of BUSCOs (This will automatically connected and dowloand the database from the BUSCO website).
-c: number of CPUs
-m: mode (options are genome, transcriptome, proteins)

```

##### *** IMPORTANT ***

BUSCO doesn't have an option to redirect the output to a different folder. For that reason, we will submit the BUSCO job from the `busco` folder. Assuming you just created the job file in the `jobs` folder:

```
cd ../busco
qsub busco_cloud_leopard.job
```

**Note about Databases:**

If you do not have internet connection on the node where running the software you can download the database and run the program offline. For instance to download the Mammalia database you can use the command `wget` and extract it.It is important to dowlod and untar the folder on your busco folder. Let's `cd` to the directory `busco` first.

	wget https://busco-data.ezlab.org/v5/data/lineages/mammalia_odb10.2021-02-19.tar.gz
	tar -zxf aves_odb9.tar.gz

The command to run busco will have to change to: 

```
busco  -o clouded_leopard -i /path/to_assembly/mNeoNeb1.pri.cur.20220520.fasta -l mammalia_odb10 -c $NSLOTS -m genome --offline --download_path /path/to/datasets
```


### Run Bloobtools

BlobTools is a command line tool designed for interactive quality assessment of genome assemblies and contaminant detection and filtering.


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
blastn -db /data/genomics/db/ncbi/db/latest_v4/nt/nt -query /path/to_assembly/mNeoNeb1.pri.cur.20220520.fasta -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 20 -max_hsps 1 -evalue 1e-20 -num_threads $NSLOTS -out clouded_leopard_blast.out
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
minimap2 -ax map-hifi -t 20 /path/to_assembly/mNeoNeb1.pri.cur.20220520.fasta  /path/to_each/hifi.fastq | samtools sort -@20 -O BAM -o cloud_leopard_sorted.bam -
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

#### Creating blobtools database

Now that we have the blast and mapping results we can create the BlobTools database. This can take a few minutes depending on how much coverage you have for your genome assembly.

- Queue: medium
- PE: multi-thread
- Number of CPUs: 1
- Memory: 6G (6G per CPU, 6G total)
- Module: `module load bio/blobtools`


- Commands:

```
blobtools create -i /path/to_assembly/mNeoNeb1.pri.cur.20220520.fasta -b clouded_leopard_mapped.bam -t clouded_leopard_blast.out -o my_first_blobplot
```

##### Explanation:

```
-i: genome assembly (fasta)
-b: mapped reads to genome assembly (bam)
-t: hits output file from a search algorith (i.e blastn). hit file is a TSV file which links sequence IDs in a assembly to NCBI TaxIDs, with a given score.
-o: path and/or name of the blobtools database.

```


#### blob and cov plots

Once you have a BlobDir database, we can plot the blobplot and the covplot. Since this is not computationally speaking intense, we can use an interactive node to run blobtools plot command.

```
qrsh
module load bio/blobtools
mkdir plots
blobtools plot -i my_first_blobplot.blobDB.json -o plots/

```

This comand generates three files:

* my_first_blobplot.blobDB.json.bestsum.phylum.p7.span.100.blobplot.bam0.png
* my_first_blobplot.blobDB.json.bestsum.phylum.p7.span.100.blobplot.read_cov.bam0.png
* my_first_blobplot.blobDB.json.bestsum.phylum.p7.span.100.blobplot.stats.txt

Please download these files to your machine. Remember that you can used the ffsend module.

* Is your genome assembly contaminated?
* if yes, What taxa are contaminant of your assembly?
* All your reads are mapping to your genome?


### Run Bloobtools2

Similar to blobtools 1.1, blobtools2 requires an assembly (Fasta), blast hit file (blast.out) and a mapping reads file (bam). You can see above on blobtools section how to create those files.

#### Creating blobtools2 data base

Now that we have the blast and mapping results we can create the BlobTools2 database. The minimum requirement to create a new database with BlobTools2 is an assembly FASTA file. This runs very fast so we do can use an interactive node.

```
qrsh -pe mthread 3
module load bio/blobtools/2.6.3 
blobtools create --fasta /path/to_assembly/mNeoNeb1.pri.cur.20220520.fasta clouded_leopard_blobt_nb
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

The first step to run Repeatmodeler is that you need to build a Database. The Build Database step is quick (several seconds at most). you can ither use and interactive node or a job file.

#### Job file: repeatmodeler_database.job
- Queue: medium
- PE: multi-thread
- Number of CPUs: 1
- Memory: 10G
- Module: `module load bio/repeatmodeler`
- Commands:

```
BuildDatabase -name cloud_leopard /path/to_assembly/mNeoNeb1.pri.cur.20220520.fasta
# usage:BuildDatabase -name {database_name} {genome_file-in_fasta_format}
```
##### Explanation:
```
-name: name to be given to the database
genome file in fasta format
```

##### Output files:
- several files with a specific structure to be populate and used by repeatmodeler command.

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
RepeatModeler -database clouded_leopard -pa 36 -engine ncbi > repeatmodeler_cl_out.log

```

##### Explanation:
```
-database: The prefix name of a XDF formatted sequence database containing the genomic sequence to use when building repeat models.
-pa: number of cpus
-engine: The name of the search engine we are using. I.e abblast/wublast or ncbi (rmblast version).
#Note the new version of repeatmodeler calls another program called LTRStruct.
-LTRStruct: enables the optional LTR structural finder.
```

##### Output files:
- consensi.fa.classified: complete database to be used in RepeatMasker.

The last step to get a repeat annotation is to run ReapeatMasker. 


#### Job file: repeatmasker_cloud_leopard.job
- Queue: high
- PE: multi-thread
- Number of CPUs: 30
- Memory: 10G
- Module: `module load bio/repeatmodeler`
- Commands:

```
# usage: RepeatMasker -pa 30 -gff -lib {consensi_classified} -dir {dir_name} {genome_in_fasta}

RepeatMasker -pa $NSLOTS -xsmall -gff -lib consensi.fa.classified -dir ../repeatmasker /path/to_assembly/mNeoNeb1.pri.cur.20220520.fasta
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
- mNeoNeb1.pri.cur.20220520.fasta.tbl: summary information about the repetitive elements
- mNeoNeb1.pri.cur.20220520.masked.fasta: masked assembly (in our case, softmasked)
- mNeoNeb1.pri.cur.20220520.fasta.out: detailed information about the repetitive elements, including coordinates, repeat type and size.


If you do not want to run Repeatmodeler for discovery of new possible repeat clases you can use Repeatmasker with available repeatdatabases.
lets create a new folder to run this job and to hold the results.



#### Job file: repeatmasker_alone.job
- Queue: medium
- PE: multi-thread
- Number of CPUs: 10
- Memory: 4G
- Module: `module load bio/repeatmasker`
- Commands:

```
RepeatMasker -species mammalia -pa $NSLOTS -xsmall -dir ../repeatmasker /path/to_assembly/mNeoNeb1.pri.cur.20220520.fasta

```
##### Explanation:
```
-species: species/taxonomic group repbase database (browse available species here: 
 https://www.girinst.org/repbase/update/browse.php)
-pa: number of cpus
-xsmall: softmasking (instead of hardmasking with N)
-dir ../repeatmasker: writes the output to the directory repeatmasker

```

##### Output files:
- mNeoNeb1.pri.cur.20220520.fasta.tbl: summary information about the repetitive elements
- mNeoNeb1.pri.cur.20220520.fasta.masked: masked assembly (in our case, softmasked)
- mNeoNeb1.pri.cur.20220520.fasta.out: detailed information about the repetitive elements, including coordinates, repeat type and size.

##### About the species:
- You can use the script `queryTaxonomyDatabase.pl` from the RepeatMasker module to search for your species of interest. 

	`queryTaxonomyDatabase.pl -species cat`
	
	**Output:**
	
	```
	RepeatMasker Taxonomy Database Utility
	======================================
	Species = cat
	Lineage = Felis catus
	          Felis
	          Felinae
	          Felidae
	          Feliformia
	          Carnivora
	          Laurasiatheria
	          Boreoeutheria
	          Eutheria
	          Theria Mammalia
	          Mammalia
	          Amniota
	          Tetrapoda
	          Dipnotetrapodomorpha
	          Sarcopterygii
	          Euteleostomi
	          Teleostomi
	          Gnathostomata vertebrate
	          Vertebrata Metazoa
	          Craniata chordata
	          Chordata
	          Deuterostomia
	          Bilateria
	          Eumetazoa
	          Metazoa
	          Opisthokonta
	          Eukaryota
	          cellular organisms
	          root
	
	```	

	Those results give you an idea of how your taxon of interest is hierarchically organized inside the repeat database. 
	

- You can also see the repeat library available for your species.

	`queryRepeatDatabase.pl -species cat`
	
This will actually print the entire repeat library in fasta format, which is not very practical. To count how many entries exist, you can pipe that output and use grep to count the number of sequences.
	
	`queryRepeatDatabase.pl -species cat | grep -c ">" `

	**Output:**
	
	```
	queryRepeatDatabase
	===================
	RepeatMasker Database: RepeatMaskerLib.embl
	RepeatMasker Combined Database: Dfam_3.0
	Species: cat ( felis catus )
	Warning...unknown stuff <
	>
	782
	```
	
	The thing is: this result doesn't necessarily mean that there are 782 repetitive elements specifically for cat. Let's test it with Feliformia:

	`queryRepeatDatabase.pl -species Feliformia | grep -c ">"`
	
	**Results:**
	
	```
	queryRepeatDatabase
	===================
	RepeatMasker Database: RepeatMaskerLib.embl
	RepeatMasker Combined Database: Dfam_3.0
	Species: Feliformia ( feliformia )
	Warning...unknown stuff <
	>
	782
	```
	
Same number, right?
Here's what this script is giving you: it will use the most closely related library within the taxonomic hierarchy for that taxon. So, there's nothing specific for cat that is not present in the Feliformia library. Dr. Vanessa Gonzalez wrote this cool script that searches for all entries in the taxonomy. The script can be found in `/data/genomics/workshops/STRI2020/Repbase_RepeatQuery_Taxonomy_MTNT.sh` and it will output a file with a list of all taxonomy levels and the number of repeats in each library. To run this script, copy it to your `repmasker` folder and do:

`./Repbase_RepeatQuery_Taxonomy_MTNT.sh cat`

This script will output two files: `cat_tax.txt` (taxonomy) and `cat_RepBase_RepQuery.txt` (the file we actually want). If we use the command `cat` to print the contents of the file `cat_RepBase_RepQuery.txt`, here's what we will see:

**Output:**

```
Taxonomic query = 'Felis catus'
Repeats in RepBase = 782
Taxonomic query = 'Felis'
Repeats in RepBase = 782
Taxonomic query = 'Felinae'
Repeats in RepBase = 782
Taxonomic query = 'Felidae'
Repeats in RepBase = 782
Taxonomic query = 'Feliformia'
Repeats in RepBase = 782
Taxonomic query = 'Carnivora'
Repeats in RepBase = 782
Taxonomic query = 'Laurasiatheria'
Repeats in RepBase = 782
Taxonomic query = 'Boreoeutheria'
Repeats in RepBase = 1883
Taxonomic query = 'Eutheria'
Repeats in RepBase = 1888
Taxonomic query = 'Theria Mammalia'
Repeats in RepBase = 1888
Taxonomic query = 'Mammalia'
Repeats in RepBase = 1888
Taxonomic query = 'Amniota'
Repeats in RepBase = 2023
Taxonomic query = 'Tetrapoda'
Repeats in RepBase = 2023
Taxonomic query = 'Dipnotetrapodomorpha'
Repeats in RepBase = 2023
Taxonomic query = 'Sarcopterygii'
Repeats in RepBase = 2023
Taxonomic query = 'Euteleostomi'
Repeats in RepBase = 3915
Taxonomic query = 'Teleostomi'
Repeats in RepBase = 3915
Taxonomic query = 'Gnathostomata vertebrate'
Repeats in RepBase = 3915
Taxonomic query = 'Vertebrata Metazoa'
Repeats in RepBase = 3915
Taxonomic query = 'Craniata chordata'
Repeats in RepBase = 3915
Taxonomic query = 'Chordata'
Repeats in RepBase = 3915
Taxonomic query = 'Deuterostomia'
Repeats in RepBase = 3915
Taxonomic query = 'Bilateria'
Repeats in RepBase = 6244
Taxonomic query = 'Eumetazoa'
Repeats in RepBase = 6244
Taxonomic query = 'Metazoa'
Repeats in RepBase = 6244
Taxonomic query = 'Opisthokonta'
Repeats in RepBase = 6244
Taxonomic query = 'Eukaryota'
Repeats in RepBase = 6244
Taxonomic query = 'cellular organisms'
Repeats in RepBase = 6244
Taxonomic query = 'root'
Repeats in RepBase = 6244
```

As you can see, the number of repetitive elements in the library is the same until Laurasiatheria (which includes carnivorans, ungulates, shrews, bats, whales, pangolins...). Mammals are very well represented compared to other groups, so this is a good thing to keep in mind when choosing a species for your analysis.

### Run GeMoMa

Gene Model Mapper (GeMoMa) is a homology-based gene prediction program. GeMoMa uses the annotation of protein-coding genes in a reference genome to infer the annotation of protein-coding genes in a target genome. Thus, GeMoMa uses amino acid and intron position conservation to create gen models. In addition, GeMoMa allows to incorporate RNA-seq evidence for splice site prediction. (see more in [GeMoMa](http://www.jstacs.de/index.php/GeMoMa-Docs)).


In this section, we show how to run the main modules and analysis with GeMoMa using real data. We will start the genome annotation for our Clouded Leopard final genome. Since GeMoMa uses the annotation of coding genes from reference genome the first step is to download some reference genomes (taxa_assembly.fasta) with its genome annotation (taxa_annotation.gff). We have done this for you already to speed up the process. The reference genomes that we have are located here ```/data/genomics/workshops/smsc_2023/ref_genomes```. We have a total of 6 reference genomes for closely related taxa:

* Neofelis nebulosa (Clouded Leopard v1)
* Neofelis diardi
* Felis catus (Domestic Cat)
* Canis familiaris (Domestic Dog)
* Homo sapiens (Human)

If you wish to download other genomes you need to search on known genome databases or genome hubs. For instance, you can search on NCBI genomes and get the ftp url and download the files (gff and fasta) using wget. 

<details><summary>SOLUTION</summary>
<p>

```
cd /ref_genomes/
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_XXXX/GCF_00000XXXX_genomic.fna.gz
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_XXXX/GCF_00000XXXX_genomic.gff.gz
cd ..
```
</p>
</details>


In General GeMoMa workflow performs 7 stepts: Extract RNA-seq evidence (ERE), DenoiseIntrons, Extractor, external search (tblastn or mmseqs), Gene Model Mapper (GeMoMa), GeMoMa Annotation Filter (GAF), and AnnnotationFinalizer. You can run this modules one by one or you can use the GeMoMa Pipeline, which is a multi-threaded tool that can uses all compute cores on one machine. 

#### Job file: busco_cloud_leopard.job
- Queue: high
- PE: multi-thread
- Number of CPUs: 40
- Memory: 10G (10G per CPU, 400G total)
- Module: `module load bio/gemoma/1.9`
- Commands:


```
set maxHeapSize 400000
GeMoMa -Xmx400G GeMoMaPipeline threads=$NSLOTS outdir=../gemoma/ GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO p=false o=true t=../assembly/mNeoNeb1.pri.cur.20220520.fasta s=own i=cat a=/pool/genomics/ariasc/SMSC_2023/ref_genomes/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.gff g=/pool/genomics/ariasc/SMSC_2023/ref_genomes/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna s=own i=dog a=/pool/genomics/ariasc/SMSC_2023/ref_genomes/GCF_014441545.1_ROS_Cfam_1.0_genomic.gff.gz g=/pool/genomics/ariasc/SMSC_2023/ref_genomes/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna.gz s=own i=human a=/pool/genomics/ariasc/SMSC_2023/ref_genomes/GCF_000001405.40_GRCh38.p14_genomic.gff.gz g=/pool/genomics/ariasc/SMSC_2023/ref_genomes/GCF_000001405.40_GRCh38.p14_genomic.fna.gz s=own i=n_diardi a=/pool/genomics/ariasc/SMSC_2023/ref_genomes/diardi_annotation.gff  g=/pool/genomics/ariasc/SMSC_2023/ref_genomes/neofelis_diardi_masked.fasta s=own i=n_nebulosa_v1 a=/pool/genomics/ariasc/SMSC_2023/ref_genomes/hic_nebulosa_annotation.gff g=/pool/genomics/ariasc/SMSC_2023/ref_genomes/clouded_leopard_HiC.fasta
```

##### Explanation:
```
threads: nnumber of CPUs
AnnotationFinalizer.r: to rename the predictions (NO/YES)
p: obtain predicted proteines (false/true)
o: keep and safe  individual predictions for each reference species (false/true)
t: target genome (fasta)
outdir: path and prefix of output directory
i: allows to provide an ID for each reference organism
g: path to individual reference genome
a: path to indvidual annotation genome
```

**Tips and Notes about GeMoMa PipeLine:**

* Depending on your default settings and the amount of data, you might have to increase the memory that can be used by setting the VM arguments for initial and maximal heap size, e.g., -Xms200G -Xmx400G.
* GeMoMaPipeline writes its version and all parameters that are no files at the beginning of the predicted annotation. This allows to check parameters at any time point.
* If GeMoMaPipline crashes with an Exception, the parameter restart can be used to restart the latest GeMoMaPipeline run, which was finished without results, with very similar parameters. This allows to avoid time-consuming steps like the search that were successful in the latest GeMoMaPipeline run.
* If you like to use mapped RNA-seq data, you have to use the parameters r and ERE.m:

```
GeMoMa -Xmx400G GeMoMaPipeline  threads=$NSLOTS AnnotationFinalizer.r=NO p=false o=true t=clouded_leopard.fna.gz outdir=output/ r=MAPPED ERE.m=<SAM/BAM> a=NCBI/GCF_XXXX_genomic.gff.gz g=NCBI/GCF_XXXX_genomic.fna.gz

```
* If you like to combine GeMoMa predictions with given external annotation external.gff, e.g., from ab-initio gene prediction, you can use the parameter e:

```
GeMoMa -Xmx400G GeMoMaPipeline threads=$NSLOTS AnnotationFinalizer.r=NO o=true t=clouded_leopard.fna.gz outdir=output/ 
p=true i=<REFERENCE_ID> a=NCBI/GCF_XXXX_genomic.gff.gz g=NCBI/GCF_XXXX_genomic.fna.gz ID=<EXTERNAL_ID> e=external.gff
```
* GeMoMa pipeline is extremely complex with many parameter for each of the module. Here is the full description of the paramenter for the pipeline script but also for each of the modules. [link to GeMoMa documentation](http://www.jstacs.de/index.php/GeMoMa-Docs)

