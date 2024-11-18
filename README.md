# RNASeqProject
This project is using C. albicans genomic data to learn more and practice bioinformatics.
# The workflow of this project is as follows:
## FastQC
Run fastqc on the initial data to assess its quality
## Trimmomatic cleaning
Based on raw fastqc results, determine optimal Trimmomatic parameters (as specified in the Trimmomatic SBATCH script file) and run Trimmomatic on the raw data
## Analysis with FastQC
Use fastqc to assess the cleaned results. Once the cleaned results are deemed high-quality, use bowtie2 to align the data with the C. albicans reference genome
## Use bowtie2
Assess bowtie2 results. The bowtie2 results yielded an alignment of 88.51% (of alignments that occurred exactly one time). At this stage in the process, my results were slightly different from Grace's results, who is also examining the same WTA2_1 and WTA2_2 files. 
