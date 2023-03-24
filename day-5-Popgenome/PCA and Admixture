# Population Genomics - SMSC (Day 5)

<!-- TOC depthFrom:2 -->

 * [Population Structure: Runing A Principal Componet Analisis](#Population-Structure:-Runing-A-Principal-Componet-Analisis)
 
<!-- /TOC -->


### Population Structure: Runing A Principal Componet Analisis

Now we have full filtered VCF file with variants for our seven Clouded Leopards. Nos we will investigate population structure using principal components analysis. Examining population structure can give us insight into origin and history of populations. PCA is a Model-free method used to asses population structure and ancestry.

PCAs identify the main axes of variation in a dataset. In the context of genomic data, PCA summarizes the major axes of variation in allele frequencies and then produces the coordinates of individuals along these axes.
. 
To calculate our PCA is we are going to used a very versatil program Called plink. Plink was design for whole genome association analysis and is powderful and easy to used.

Before running PLink we need to prepared our data to not violate assumptions for PCA analisis. The most important assumptions is that every SNP in independent. This is clearly not the case for genomic data sets since allele frequencies are correlated du to physical linkage. So we need to first prune variants that are in physical linkage.

Let's change directory to our pop_genomics folder and copy the variant calling file (vcf) to this folder. The file is located here: 'pool/genomics/ariasc/SMSC_2023/mapping/raw_c_leopard_vcf/clouded_leopard_data.vcf'

Plink is very fast so we do not need to set up job. We will be using an interactive node.

```
qrsh
module load bio/plink/
plink --vcf clouded_leopard_data.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out clouded_leopard_plink
```

##### Explanation:

```
--vcf: Path to input VCF file.
--double-id: this tells to plink to duplicate the id of our samples its assumes that 
--allow-extra-chr: allow additional chromosomes beyond the human chromosome set. Plink by default expectes that the data is human (i.e. chromosomes 1-22 and X chromosome).
--also necessary to set a variant ID for our SNPs. Human datas sets often have annotated SNP names and so plink will look for these. We do not have them so instead we set ours to default to chromosome:position which can be achieved in plink by setting the option @:#.
--indep-pairwise - This is the command that performes the linkage pruning. The first argument, 50 denotes we have set a window of 50 Kb. The second argument, 10 is our window step size - meaning we move 10 bp each time we calculate linkage. Finally, we set an r2 threshold - i.e. the threshold of linkage we are willing to tolerate. Here we prune any variables that show an r2 of greater than 0.1.

```

Now we run our PCA analysis with plink using the following command.

```
plink --vcf clouded_leopard_data.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract clouded_leopard_plink.prune.in \
--make-bed --pca --out clouded_leopard_pca
```

##### Explanation:

```
--extract: this tells plink to extract only these positions from our VCF.
--make-bed: this is to produce and output file that we will used in our next population analysis (model based approach with admixture).
--pca: this tells plink to calculate a principal components analysis.

```


##### Dowload files to your machine

after the command finish, find the results and download the following two files to your local machine:

```
clouded_leopard_pca.eigenval  
clouded_leopard_pca.eigenvec
```

Remember to use ffsend. Tip: You need to load the module ```bioinformatics/ffsend``` and then used the command ```ffsend upload <FILE>```.

### Visualializing resutls in R and Rstudio

Now that we have the files on our machines lets plot our PCA results. For this we will use R and Rstudio.
The first step is to create and Rstudio project were you will performer your analysis. In your R studio project, open a new rmarkdown document or and R script to write and comment your R commands. You can also move the results files from plink  to your project folder.

First let's load the tydiverse library. This library is very powerful and has many functions to manipulate data in R.

```
library (tydyverse)
```

Now lets read our data using the read table function.

```
pca <- read_table2("clouded_leopard_pca.eigenvec", col_names = FALSE)
eigenval <- scan("clouded_leopard_pca.eigenval")
```
Next, we need to start cleaning our data frame table. Look at your data frame. We will delete duplicated columns (i.e. individual ID) and also we will give a proper name to  the columns in the pca data frame.

```
pca <- pca[,-1]
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
```

Let's now organize individuals by its name and its population/species.
for this you will create some vectors that will hold that information (i. nm for name and pop for population).

```
# name
nm <- rep(NA, length(pca$ind))
nm[grep("NN114296", pca$ind)] <- "pepe"
nm[grep("NN114297", pca$ind)] <- "rosa"
nm[grep("NN114393", pca$ind)] <- "shakiral"
nm[grep("NN115950", pca$ind)] <- "carlos"
nm[grep("NN190240", pca$ind)] <- "mike"
# population
pop <- rep(NA, length(pca$ind))
pop[grep("NN114296", pca$ind)] <- "cl"
pop[grep("NN114297", pca$ind)] <- "cl"
pop[grep("NN114393", pca$ind)] <- "cl"
pop[grep("NN115950", pca$ind)] <- "cl"
pop[grep("NN190240", pca$ind)] <- "cl"
```
if you want to plot later by different colors and separted by both color and population you can create a thir vector do:

```
nm_pop <- paste0(nm, "_", pop)
```

After creating this variables you can reorganize your data frame into a new one.

```
pca2 <- as.tibble(data.frame(pca, nm, pop, nm_pop))
```
To plot the percentage of variance explained by each principal component you need to convert your eigenvalues into percentage  and create a new data frame.

```
pve <- data.frame(PC = 1:5, pve = eigenval/sum(eigenval)*100)
```

Now that we have done our cleaning, we can plot our eignvalues and PCA.  
First we will plot the eigenvalues. 
```
plot_pve <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
plot_pve <- plot_pve  + ylab("%  of variance expl") +
theme_classic()
plot_pve
```
Cumulatively, as you can see the five PC explain 100% of the variance, but PC1 and PC2 explain about 58% of the variance. We could calculate this formaly  with the cumsum function, like this:

```
cumsum_cal<-cumsum(pve$pve)
cumsum_cal
```

Finally let's plot PC1 and PC2. The way to interpret this result is that the closer the dots, more similar the individuals are.

```
pca_plot <- ggplot(pca2, aes(PC1, PC2, col = nm, shape = pop)) + geom_point(size = 3)
pca_plot <- pca_plot + scale_colour_manual(values = c("red", "blue", "green", "purple","black"))
pca_plot <- pca_plot + coord_equal() + theme_light()
pca_plot <- pca_plot + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
pca_plot
```

 * What Happen if you plot other PCs?
 * Play with the colors and change the names?
 
### Challange

Let's play know with another an more complex file. I have put one of my data files, this represents an ample sampling of close to 60 individuals of Brachyhypopomus occidentalis across all of his natural range (From Venezuela, Costarica).


Try the run the same commands but with this data set:


<details><summary>SOLUTION</summary>
<p>

##### Plink Commands
```
# perform linkage pruning - i.e. identify prune sites
plink --vcf ../efish.vcf --double-id --allow-extra-chr  --set-missing-var-ids @:#  --indep-pairwise 50 10 0.1 --out efish

# prune and create pca
plink --vcf ../efish.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract efish.prune.in --make-bed --pca --out efish

```
 
#### R ploting Commands

```
# load tidyverse package
library(tidyverse)
# read in data
pca <- read_table2("./efish.eigenvec", col_names = FALSE)
eigenval <- scan("./efish.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
#population map
pop_map<- read.delim("pop_map4.txt", header = T)

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
pve_plot <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
pve_plot + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

#PCA_joined_popm <- left_join(pca, pop_map, by = c("ind" = "ind"))
PCA_joined_popm$loc <- as.factor(PCA_joined_popm$loc)

#PCA_joined_popm %>% rename(loc1=V2,loc2=V3)
PCA_joined_popm <- merge(pca, pop_map, by.x = "ind",  by.y = "ind", all.x = TRUE, all.y = F)

levels(PCA_joined_popm$loc)

cols <- c("Atrato"= "#009E73","Bayano"="#56B4E9","Bocas" ="#CC79A7","Calovebora"= "brown","CanalZone"="#B6DBFF","Cocle"="#F0E442","Darien"="#882255","GunaYala"="#7F7F7F", "Orinoco"="darkgreen", "SanJuan"= "green", "Lajas"="#999933")  

cols2 <- c("Bocas" ="#E69F00","Calovebora"= "#999933","Santa_Maria"="deepskyblue","Cocle"="red2","Chagres"="khaki4","GunaYala"="#7F7F7F","Bayano"="#CC79A7", "Darien"="darkgreen","P_Col"= "blueviolet", "Venezuela"="red4")  


#scale_color_manual(values=c("#009E73","#D55E00","#F0E442","#CC79A7","#882255","#56B4E9","#E69F00","#0072B2","#7F7F7F","#B6DBFF","#999933"))

L_labels= c("Bocas" ="Bocas (W.Pan.)", "Calovebora"="Calovebora (C.Pan)", "Santa_Maria"="Santa Maria (C.Pan)",  "Cocle"="Cocle Norte (N.Pan)",  "Chagres"="Chagres (C.Pan)","GunaYala"="GunaYala (N.Pan)","Bayano"="Bayano (E.Pan)", "Darien"="Darien (E.Pan)", "P_Col"= "Pacific Col", "Venezuela"="Venezuela") 

PCA_joined_popm$loc <- factor(PCA_joined_popm$loc, levels = c("Bocas","Calovebora","Santa_Maria", "Cocle", "Chagres","GunaYala",   "Bayano", "Darien", "P_Col", "Venezuela"), labels =c("Bocas" ="Bocas (W.Pan.)", "Calovebora"="Calovebora (C.Pan)", "Santa_Maria"="Santa Maria (C.Pan)",  "Cocle"="Cocle Norte (N.Pan)",  "Chagres"="Chagres (C.Pan)","GunaYala"="GunaYala (N.Pan)","Bayano"="Bayano (E.Pan)", "Darien"="Darien (E.Pan)", "P_Col"= "Pacific Col", "Venezuela"="Venezuela"))

# plot pca
library("ggplot2")
b <- ggplot(PCA_joined_popm, aes(PC1, PC2, col = loc)) + geom_point(size = 3, alpha=0.6)
b <- b + scale_colour_manual(values =cols2, name="Localities", labels=L_labels)
b <- b + coord_equal() + theme_light()
b <- b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
b

# Save file as PNG
ggsave(b, file="PCA_efish2.png", dpi = 300)

```
</p>
</details>
 
#### Runing ADMIXTURE

ADMIXTURE is a model base program. That calculates individual and population ancestries.

To run admixture let's keep using the electric fish data set. Admixture can used as input some of the files generated by plink (bed and bim files) in our previous analysis. Admixture reads bed file (A PLINK bed file is a binary biallelic genotype table) and uses information collected on the bim file too. However, since plink was develop to work with human genomes,  we need to do some cleaning. We need to change the first column of the bim file with a O. We will used the next 'awk' command and them rename back our bim file to mach the name on our bed file.

```
awk '{$1=0;print $0}' efish.bim > efish.bim.tmp
mv efish.bim.tmp efish.bim
```
Once we have our new bim corrected file,  we can run admixture let's again an interactive node and call the admixture module. 

```
qrsh
module load bio/admixture
admixture --cv efish.bed 2 > log2.out
```

Previous admixture command test specifically ancestri assuming that there was two genomic cluster. We now need to test admixture and acestri on a range of K values to compare and see which K value fits better our data set. Since this is a very simple command we can use a for loop that iterates K values and gathers the results in differen log files.

```
#Runing a loop for different K

#!/bin/sh
for i in {2..15}
do
admixture --cv efish.bed $i > log${i}.out
done
```
Admixture produced 2 files. A .Q file, which contains cluster assignments for each individual and, a .P file , which contains for each SNP the population allele frequencies. Lets used grep and awk together to identify the best value of K clusters. The best value of K will be the lowest cross-validation error. Thus, we need to collect the cv errors. Below are three different ways to extract the number of K and the CV error for each corresponding K. 

```
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > efish.cv.error
#grep "CV" *out | awk '{print $3,$4}' | cut -c 4,7-20 > efish.cv.error
#awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > efish.cv.error
```

* which is the best K for this data?


We also need to extract a list of the individuals used in the analysis. This will later help when ploting the results.

```
#list of individuals in the analysis
awk '{split($1,name,"."); print $1,name[2]}' efish.nosex > efish.list

```

### Plotting admixture results in R.

Download to your machine que .Q file for the best number of Clusters and used the next script to plot your results. Also dowload the .list file


```
# read .Q file
tbl<-read.table("efish.10.Q")
#Simple plot
barplot(t(as.matrix(tbl)), col=rainbow(10),
               xlab="Individual #", ylab="Ancestry", border=NA)

#pop_map<- pop_map %>% rename(ind=V1,loc=V2, loc2=V3)

pop_map_admix <-read.table("efish.list")

Popadmix_order_joined_popm <- left_join(pop_map_admix, pop_map, 
              by = c("V1" = "ind"))
Popadmix_order_joined_popm <- Popadmix_order_joined_popm %>% rename("ind"="V1")

tbl_pop_map<- cbind(Popadmix_order_joined_popm,tbl)

# Melt (reshape data from wide format to long format).
library(reshape2)
tbl_pop_gather = melt(tbl_pop_map, id.vars=c("ind", "pop", "loc"), variable.name="Ancestry", value.name="Fraction")
#change values within df
#levels(tbl_pop_gather$Ancestry)[levels(tbl_pop_gather$Ancestry)=="V1"] <- "K1"


# Simple stacked bar plot:
#col=c("red", "blue","darkblue", "green", "darkgreen", "black", "orange", "purple", "pink","brown")

tbl_pop_gather$loc <- as.factor(tbl_pop_gather$loc)

levels(tbl_pop_gather$loc)

tbl_pop_gather$loc = factor(tbl_pop_gather$loc,levels=c("Venezuela","P_Col", "Darien","Bayano","GunaYala","Chagres", "Cocle","Santa_Maria","Calovebora","Bocas"),labels=c("V","Col","Da","By","GY","Ch", "CN","SM", "Cv","BT")) 

#cols2 <- c("Bocas" ="#E69F00","Calovebora"= "#999933","Santa_Maria"="deepskyblue","Cocle"="red2","Chagres"="khaki4","GunaYala"="#7F7F7F","Bayano"="#CC79A7", "Darien"="darkgreen","P_Col"= "blueviolet", "Venezuela"="red4")  


library("ggplot2")
p = ggplot(tbl_pop_gather, aes(x=ind, y=Fraction, fill=Ancestry)) +
    geom_bar(stat="identity", position="stack") +
    facet_grid(. ~ loc, drop=TRUE, space="free", scales="free")

p2 = p + theme(panel.grid=element_blank())+theme(axis.title.x=element_blank())+ theme(axis.text.x=element_blank()) + theme(axis.ticks.x=element_blank()) +  
#+ scale_fill_manual(values=col) 
theme(legend.position="none")
p2<- p2 + scale_fill_manual(values=c("blueviolet","green","red2","#E69F00","darkgreen","#CC79A7", "deepskyblue","#7F7F7F","khaki4","#999933"))
p2

ggsave(p2, file="admixture_efish2.png", dpi = 300)

```
