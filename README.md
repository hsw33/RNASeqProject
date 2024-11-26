# RNASeqProject
## Project summary
Using C. albicans data, courtesy of Professor Ronda Rolfes at Georgetown University, this project aims to elucidate the genomic differences of of C. albicans grown in the presence and absence of thiamine treatment. C. albicans is a yeast species that is responsible for causing urinary tract and other genital infections, and transmission is especially common in hospital settings, where patients are, unfortunately, already susceptible to infection due to weakened immune system.
## Files used
My particular analysis was performed on two files, WTA2_1.fq.gz, WTA2_2.fq.gz. These files contain the genomic data from one isolate of C. albicans. The A2 designation means that the strain was grown in the presence of of thiamine. The 1 and 2 designations in the file title represent the read pairs. 
# Project workflow:
## FastQC
Run fastqc on the initial data to assess its quality
## Trimmomatic cleaning
Based on raw fastqc results, determine optimal Trimmomatic parameters (as specified in the Trimmomatic SBATCH script file) and run Trimmomatic on the raw data
https://docs.google.com/spreadsheets/d/1AOa-XaTzR_PKMIRQDmu8oDTmawXXnkIwEjKOQkNC7Vs/edit?usp=sharing
## Analysis with FastQC
Use fastqc to assess the cleaned results. Once the cleaned results are deemed high-quality, use bowtie2 to align the data with the C. albicans reference genome
## Use bowtie2
Assess bowtie2 results. The bowtie2 results yielded an alignment of 88.51% (of alignments that occurred exactly one time). At this stage in the process, my results were slightly different from Grace's results, who is also examining the same WTA2_1 and WTA2_2 files. 
## Use samtools and bamtools
Need to add more explanation here
## Use conda environment and run HTSeq on the bam index file
Need to add more explanation here
## Use RStudio to conduct biological analysis
Using RStudio and a script available on this page, I ran analysis to extract pertinent biological information. Specifically, I obtained a table that identified 13 genes that are differentially expressed in thiamine-present vs absent C. albicans isolates.
