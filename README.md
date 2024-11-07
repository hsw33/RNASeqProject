# RNASeqProject
This project is using C. albicans genomic data to learn more and practice bioinformatics.
The workflow of this project is as follows:
1. Run fastqc on the initial data to assess its quality
2. Based on raw fastqc results, determine optimal Trimmomatic parameters (as specified in the Trimmomatic SBATCH script file) and run Trimmomatic on the raw data
3. Use fastqc to assess the cleaned results. Once the cleaned results are deemed high-quality, use bowtie2 to align the data with the C. albicans reference genome
4. Assess bowtie2 results 
