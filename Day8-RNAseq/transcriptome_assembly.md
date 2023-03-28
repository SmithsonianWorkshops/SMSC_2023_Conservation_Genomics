# Transcriptome Assembly and Quality Analysis

The data we're going to use today live here: `/data/genomics/workshops/smsc_2023/red_siskin_RNAseq/`  
They are from the [Red Siskin Genome Project](https://www.braunlab.umd.edu/red-siskin-conservation/). They come from many different tissue types and are paired-end reads.  
Let's just list in that directory to take a look.  

Technically we could use a reference genome to help with some of these analyses, but that won't always be the case, so today we're going to do a de novo transcriptome assembly so that you know how it works.  

I've made three more directories to work with these RNAseq data: data (for trimming and reformatting the reads), trans_assembly (for doing the actual assembly), and assembly_qc (for quality assessment of the assembly).  



As with a genome assembly, we need to make sure that our reads are free of adaptors before we do the assembly. Rather than test each file ahead of time, we're going to move right to trimming. (I happen to know that they could use some adapter trimming.)  

Choose a tissue type to work with for today. Since we're not using these data for the paper, it's not a big deal if lots of folks want to work on the same tissue type, so follow your heart!  


### Trimming the RNAseq reads  

We'll use Trim Galore again, since you're already familiar with that program.  

CPUs: 1
Memory: 5G
Module: bio/trim_galore  
Basic Command: `trim_galore --paired -q 3 --dont_gzip <FILE_1.fastq> <FILE_2.fastq>` 
Options:
- paired: we're using paired-end reads  
- q: this is the Phred quality score threshold that we're using to trim the reads  
- dont-gzip: if we give the program compressed reads, it will compress the output by default. However, it will be easier to do the next part of the read processing if we let trim galore keep the reads unzipped for us.  
- file_1 and file_2: this is where the read files go, complete with their full paths   

If you like, you can also add a "-o" flag that will let you specify an output directory. It will create this directory for you if you specify one, and if you don't, it will just make an output directory in the directory from which the job is submitted.  


### Fixing the format of the read headers  

Lots of assembly programs are very picky about what the header line of the fastq files looks like. It's always a good idea to make sure there are no special characters or spaces, and that the paired reads are designated correctly.  

Look at the headers of one of the fastq files for your tissue type - all kinds of nonsense happening there. Let's fix it.  

First we'll cut off everything after the space in the header, then tag "/1" or "/2" at the end of each one, and then (because reformat.sh introduces another space when it adds the slashes *<insert eyeroll>*) we'll replace the remaining space with an underscore. This process takes more memory than you'd think, unfortunately, but takes almost no time at all.  

In the same job, we'll go ahead and take a random subsample of the reads in the forward and reverse files - let's subsample down to 15 million reads. This won't really do our assembly any favors, but it will make it finish in time to check out the quality later.  

CPUs: 1  
Memory: 16G  
Modules: bio/bbmap, bio/seqtk/  
Basic Commands:  
```  
reformat.sh in=original_fastq_file1.fq in2=original_fastq_file2.fq out=intermediate_file1.fq out2=intermediate_file2.fq trd  
reformat.sh i=intermediate1_file1.fq in2=intermediate1_file2.fq out=intermediate2_file1.fq out2=intermediate2_file2.fq addslash  
reformat.sh in=intermediate2_file1.fq in2=intermediate2_file2.fq out=formatted_file1.fq out2=formatted_file2.fq underscore  
seqtk sample -s 51 formatted_file1.fq 15000000 > sub_formatted_file1.fq  
seqtk sample -s 51 formatted_file2.fq 15000000 > sub_formatted_file2.fq  
```  
Options for reformat.sh:  
- in and out: input file and output file (make sure that you're using the file from the previous step as the input and that you create a new filename for the output)  
- trd: truncates the header line at the whitespace  
- addslash: adds the "/1" or "/2" to the end of the header lines  
- underscore: replaces whitespaces with underscores  

Options for seqtk:  
- sample: part of the command - we're taking a random subsample  
- s: this is a random seed. It doesn't matter what number you pick, but it needs to be the same one for the forward and reverse reads  
- input file: the original fastq file to be subsampled  
- 15000000: the number of reads that we want in the end  
- output file: our new, subsampled fastq file  

After this job finishes, you can `ls -lthr` in this directory and see that the file sizes of the subsampled files should be much smaller than the originals. You can also `less` into one of them to check out the new header formatting.  


### De novo transcriptome assembly  

Now we can move on to the actual assembly. 

CPUs: 24  
Memory: 8G  
Module: bio/transabyss/2.0.1  
Basic Command: `transabyss --pe <file1> <file2> --threads $NSLOTS --name <output_prefix>`  

Options for transabyss:  
- pe: paired-end reads - make sure to include both full paths  
- threads: num of CPUs  
- name: the first part of all the output files transabyss will generate for this particular assembly  


### Assembly evaluation  

We'll run busco in a really similar way to how we did it for the genome assembly, but it's generally easier and faster.  

CPUs: 12  
Memory: 4G  
Module: bio/busco/5.4.3  
Basic Command: `busco -i input_file -l aves_odb10 -o output_directory -m transcriptome -c $NSLOTS`  

Options for busco:  
- i: path to the input file - a finished transcriptome assembly  
- l: the database we want busco to use - there is a specific one for birds  
- o: the name (or path) of the output directory where busco will put all the output files  
- m: the mode, for us it's "transcriptome", but it could also be "genome" or "protein"  
- c: number of CPUs to use  

