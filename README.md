# RNA Sequencing Project for BIOL 3615 with Prof. Peter Armbruster
### October - December 2024
## Project summary
Using *C. albicans* data, courtesy of Professor Ronda Rolfes at Georgetown University, this project aims to elucidate the differential expression of genes in *C. albicans* replicates grown in the presence and absence of thiamine treatment. *C. albicans* is a yeast species responsible for  urinary tract and other genital infections, and transmission is especially common in hospital settings, where patients are, unfortunately, already susceptible to infection due to weakened immune systems.
## Files used
My particular analysis was performed on two files, WTA2_1.fq.gz, WTA2_2.fq.gz. These files contain the genomic data from one isolate of *C. albicans*. The A2 designation means that the strain was grown in the presence of thiamine. The 1 and 2 designations in the file title represent the read pairs. These file names will be important to keep track of, as they're how I designated the data that I worked on throughout this project.
## Experimental design
In the Rolfes Lab at Georgetown University, researchers grew *C. albicans* replicates in the presence and absence of thiamine, a critical B vitamin for the species. We used 12 data files, but I individually conducted analysis on two of the 12 files, WTA2_1.fq.gz, WTA2_2.fq.gz. The "_1" and "_2" designations represent the forward and reverse reads, while the WTA1/B1/C1 designation indicates that the replicates were grown in the presence of thiamine, WTA2/B2/C2 replicates were grown in the absence of thiamine.
### Data files used by the class:
#### Grown in the presence of thiamine:
WTA1_1.fq.gz, WTA1_2.fq.gz

WTB1_1.fq.gz, WTB1_2.fq.gz

WTC1_1.fq.gz, WTC1_2.fq.gz

#### Grown in the absence of thiamine:
WTA2_1.fq.gz, WTA2_2.fq.gz

WTB2_1.fq.gz, WTB2_2.fq.gz

WTC2_1.fq.gz, WTC2_2.fq.gz
# Project workflow:
The workflow of this project contained many steps. Below, I will describe each step and its utility, as well as provide relevant code and scripts.
## FastQC
I ran FastQC on the initial data to assess its quality. By looking at the initial graphical output on FastQC, it was clear that the data needed to be cleaned up to maximize its quality because several of the reported values were in the "yellow" or "red" zone, suggesting medium-to-poor quality data.
## Trimmomatic cleaning
Based on raw FastQC results, I determined the optimal Trimmomatic parameters (as specified in the Trimmomatic SBATCH script file) and ran Trimmomatic on the raw data. In particular, there was concern about the per-base sequence content, so we used the HEADCROP command to improve that metric.

[This spreadsheet](https://docs.google.com/spreadsheets/d/1AOa-XaTzR_PKMIRQDmu8oDTmawXXnkIwEjKOQkNC7Vs/edit?usp=sharing) was used to record data quality scores before and after cleaning.

See: [Trimmomatic Script](https://github.com/hsw33/RNASeqProject/blob/main/Scripts/Trimmomatic_Script)
## Analysis with FastQC
Then, I used FastQC to assess the cleaned results. The parameters used in Trimmomatic successfully  increased the overall quality of the data. 95.63% of the reads were retained after cleaning, which is an excellent result. 

After I deemed the output to be high-quality, I used bowtie2 to align the data with the *C. albicans* reference genome, linked [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000182965.3/).

## Use bowtie2
I then used bowtie2 to align my experimental reads with the reference genome. After running the program, the bowtie2 results yielded an alignment of 88.51% (of alignments that occurred exactly one time). At this stage in the process, my results were slightly different from Grace's , who is also examining the same WTA2_1 and WTA2_2 files. This is important to note because Grace and I should have identical results. This discrepancy may be due to a small difference in our Trimmomatic scripts, like the order in which we invoked the trimming commands.

See: [bowtie2 Script](https://github.com/hsw33/RNASeqProject/blob/main/Scripts/bowtie2_Script)
## Use samtools and bamtools
The output of running bowtie2 on the data was a file with the extension .sam. SAM files are extremely large files that are too cumbersome to interpret, so I converted the .sam file to a .bam file. The bam file is not human-readable because it is written in binary code, but its size is significantly reduced so I was able to use it as an input file for further analysis. At this point, I also had to create a .bam index file, which is a sorted and indexed version of our .bam file with the extension .bam.bai.
### Code to convert from .sam to .bam & to index the .bam file:
#### Load samtools
$ srun --pty bash
$ module load samtools
#### Convert sam -> bam:
$samtools view -S -b sample.sam > sample.bam
#### View to visualize alignments and header:
$ samtools view sample.bam|more 
$ samtools view –H sample.bam|more
#### Count number of alignments:
$ samtools view sample.bam|wc –l
#### Sort:
$ samtools sort sample.bam –o sample.srt.bam
#### Index:
$samtools index sample.srt.bam (output file is .bam.bai)

## Use conda environment and run HTSeq on the bam index file:
After setting up the conda environment (see initialization instructions below), I utilized the program HTSeq. I wrote an HTSeq SBATCH script as the analysis took several hours to run. HTSeq quantified the read counts of the alignment between the two sets of *C. albicans* data.

See: [HTSeq Script](https://github.com/hsw33/RNASeqProject/blob/main/Scripts/HTSeq_Script)

### Code to set up the Conda environment:
$ module load anaconda3

$ conda create --[name] htseq

$ conda init

(restart shell by logging off HPC, log back on)

$ module load anaconda3

$ conda activate htseq

$ conda install -c bioconda htseq

## Use RStudio and DESeq2 to conduct biological analysis
After the output from HTSeq was available, I downloaded the files onto my desktop computer. Using RStudio, the R package DESeq2, and a script available on this page, I analyzed the output files from HTSeq to extract pertinent biological information. Specifically, I obtained a table that identified 13 genes that are differentially expressed in thiamine-present vs. absent *C. albicans* isolates. DESeq2 was able to identify the differential expression of genes between the transcriptomes. There were only 13 genes that were significantly differentially expressed, and this data was consolidated into a table, as well as a volcano plot. 

Separately, information about the transcriptome as a whole was charted in the Principle Components Analysis (PCA) plot. 

All of this data, including the R script, is available in the [RStudioData](https://github.com/hsw33/RNASeqProject/tree/main/RStudioData) folder. 
## Interpretation of data
We used multiple programs and resources to make sense of our data and interpret its biological relevance.

The output from DESeq2 returned 13 genes that were significantly differentially expressed. Using a file on the HPC that had specific annotations for each gene of significance, I used the following command:

$ grep -wFf signifcant_geneIDs GCF_000182965.3_ASM18296v3_genomic.gtf | grep "protein_coding" | cut -f9| cut -d ";" -f1,3,5 > signif_gene_annotation_info

This command generated a text file that I input into my spreadsheet containing information on the 13 genes. The gene names and gene IDs were now accessible in [this Excel file](https://github.com/hsw33/RNASeqProject/blob/main/RStudioData/Updated%20Biological%20Functions%20of%20Signif_TH-vTH%2B.xlsx) and I began to research the role of these genes in *C. albicans*.

### DESeq2 interpretation with gene function investigation
DESeq2 identified 13 genes that were significantly differentially expressed in the thiamine-absent condition. Using the [Candida Genome Database](http://www.candidagenome.org/) and [UniProt](https://www.uniprot.org/) I researched the proposed function of each gene. Based on my research, many of the differentially expressed genes were involved in thiamine biosynthesis, transmembrane transport, and biofilm synthesis. Of these 13 genes, 4 coded for uncharacterized proteins, and it was more challenging to elucidate their function. Thus, the remaining analysis will focus on the 9 genes with better-understood functions.

Considering the remaining genes that coded for identifiable proteins, I expected to see upregulation of genes involved in pathways related to thiamine synthesis, and indeed, the results confirmed my hypothesis. Thiamine is an essential B vitamin for metabolism in many organisms, including *C. albicans*, and in its environmental absence, the organism is expected to make its own thiamine to survive. This trait is evolutionarily advantageous, thus, it is not surprising that *C. albicans* can synthesize thiamine. 

The upregulation of genes involved in transmembrane transport is also consistent with mechanisms the organism would employ to sustain its life. It is energetically advantageous for cells to utilize as much material as they can from their environment, rather than synthesize it from scratch, which is a much more energetically-costly strategy. As such, I would expect to see the organism increase its rate of cellular transmembrane transport to be energetically efficient and increase its survival odds. 

Lastly, some of these genes were associated with the production of spider biofilm. While this result was not as straightforward as the others, with some more investigation, it became clear that the production of biofilm allows *C. albicans* to create colonies supported by cellular adhesion and aid in the multiplication and reproduction of cells. Thus, the production of biofilm is another evolutionarily advantageous strategy for the species. 

### Gene ontology analysis
As the last step of this project, I used Gene Ontology (GO) and PANTHER software to further interrogate my data. Specifically, the goal of this step was to determine if there more genes with the same GO classification than expected by chance (i.e., by random sampling of genes in genome). *C. albicans* has 6,035 genes, and genome-wide frequency is calculated by dividing the number genes in GO category in genome by the total number of genes in the genome. In this case, that number is 6,035.

According to PANTHER, genes involved in the thiamine biosynthetic process had a greater than 100-fold enrichment. Similarly, genes related to the pyridoxine biosynthetic process also had a greater than 100-fold enrichment. These results suggest that in the thiamine absent-condition, there are more genes with these GO classifications than we would expect to see by random chance in *C. albicans*. This result is also logical: the more genes that can facilitate thiamine biosynthesis (whether by synthesizing thiamine itself or by synthesizing one of its precursors, pyridoxine) the more successful the organism will be at surviving in an vitamin-depleted environment.
