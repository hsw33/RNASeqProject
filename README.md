# RNASeqProject
## Project summary
Using C. albicans data, courtesy of Professor Ronda Rolfes at Georgetown University, this project aims to elucidate the genomic differences of of C. albicans grown in the presence and absence of thiamine treatment. C. albicans is a yeast species responsible for causing urinary tract and other genital infections, and transmission is especially common in hospital settings, where patients are, unfortunately, already susceptible to infection due to weakened immune systems.
## Files used
My particular analysis was performed on two files, WTA2_1.fq.gz, WTA2_2.fq.gz. These files contain the genomic data from one isolate of C. albicans. The A2 designation means that the strain was grown in the presence of thiamine. The 1 and 2 designations in the file title represent the read pairs. These file names will be important to keep track of, as they're how I designated the data that I worked on throughout this project.
# Project workflow:
## FastQC
I ran FastQC on the initial data to assess its quality. By looking at the initial graphical output on FastQC, it was clear that the data needed to be cleaned up to maximize its quality because several of the reported values were in the "yellow" or "red" zone, suggesting medium-to-poor quality data.
## Trimmomatic cleaning
Based on raw FastQC results, I determined optimal Trimmomatic parameters (as specified in the Trimmomatic SBATCH script file) and ran Trimmomatic on the raw data. In particular, there was concern about the per base sequence content, so we used the HEADCROP command to improve that metric, in particular.

The following spreadsheet was used to record data quality scores before and after cleaning.
https://docs.google.com/spreadsheets/d/1AOa-XaTzR_PKMIRQDmu8oDTmawXXnkIwEjKOQkNC7Vs/edit?usp=sharing

See: Scripts/Trimmomatic_Script
## Analysis with FastQC
Then, I used FastQC to assess the cleaned results. The parameters used in Trimmomatic were successful in increasing the overall quality of the data. After I deemed the results to be high-quality, I used bowtie2 to align the data with the C. albicans reference genome, which is linked below.
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000182965.3/
## Use bowtie2
I then used bowtie2 to align my experimental reads with the reference genome. After running the program, the bowtie2 results yielded an alignment of 88.51% (of alignments that occurred exactly one time). At this stage in the process, my results were slightly different from Grace's results, who is also examining the same WTA2_1 and WTA2_2 files. This is of importance because we should have exactly identical results. This discrepancy may be due to a small difference in our Trimmomatic scripts, like the order in which we invoked the trimming commands.

See: Scripts/bowtie2_Script
## Use samtools and bamtools
The output of running bowtie2 on the data was a file with the extension .sam. SAM files are extremely large files that would be too cumbersome to interpret, so I converted the .sam file to a .bam file. The bam file is not human-readable because it is written in binary code, but its size is significantly reduced so I was able to use it as an input file for further analysis. At this point, I also had to create a .bam index file, which is a sorted and indexed version of our .bam file.
## Use conda environment and run HTSeq on the bam index file
Using the conda environment, I set up HTSeq. Then, I wrote an HTSeq SBATCH script as the analysis took several hours to run. HTSeq quantified the read counts of the alignment between the experimental genome and the reference genome.
## Use RStudio and DESeq2 to conduct biological analysis
After the output from HTSeq was available, I downloaded the files onto my desktop computer. Using RStudio, DESeq2, and a script available on this page, I ran an analysis on the output files from HTSeq to extract pertinent biological information. Specifically, I obtained a table that identified 13 genes that are differentially expressed in thiamine-present vs absent C. albicans isolates. DESeq2 was able to identify the differential expression of genes between the two genomes. There were only 13 genes that were differentially expressed, and this data was consolidated into a table, as well as into a principle component analysis plot and a Volcano plot.
## Interpretation of data
The output from DESeq2 returned 13 genes that were significantly differentially expressed. Using a file on the HPC that had specific annotations for each gene of significance, I used the following command:

grep -wFf signifcant_geneIDs GCF_000182965.3_ASM18296v3_genomic.gtf | grep "protein_coding" | cut -f9| cut -d ";" -f1,3,5 > signif_gene_annotation_info

This command generated a text file that I input into my spreadsheet containing information on the 13 significant genes. The gene names and gene IDs were now accessible in this Excel file (RStudioData/Updated Biological Functions of Signif_TH-vTH+.xlsx) and I began to research the role of these genes in C. albicans.
