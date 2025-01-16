#Data Downloading

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR152/056/SRR15299556/SRR15299556_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR152/056/SRR15299556/SRR15299556_2.fastq.gz

#Extraction of fastq files
gunzip SRR15299556_1.fastq.gz SRR15299556_2.fastq.gz
#OR
gunzip *.gz

##Quality control using fastQC

fastqc SRR15299556_1.fastq SRR15299556_2.fastq

#OR
fastqc *.fastq

#There are three parameters we have to check to deicide the quality
#1.Per base sequence quality
#2. Overrepresneted sequence
#3. Adapter content

##Trimming of data using cutadapt(for Paired end data)
cutadapt -b GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG -B GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG -q 30,30 -o Trimmed_SRR15299556_1.fastq -p Trimmed_SRR15299556_2.fastq SRR15299556_1.fastq SRR15299556_2.fastq

#Trimming for single end data
cutadapt -b seq_to_trim -q 30,30 -o Trim_sample.fastq sample.fastq 

#Download refgenome from UCSC
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz

#Again check the quality
fastqc Trimmed_SRR15299556_1.fastq Trimmed_SRR15299556_2.fastq

##Aligning of reads to ref genome using BWA tool
#But alignment of reads to ref genome has two steps
#A. Indexing of ref genome

bwa index -a bwtsw genome.fa

#Aligning of reads to indexed genome(MEM= Maximal Exact Matches)

bwa mem genome.fa Trimmed_SRR15299556_1.fastq Trimmed_SRR15299556_2.fastq > bwa_SRR15299556.sam

#Conversion of SAM into sorted BAM
samtools sort --threads 2 bwa_SRR15299556.sam > sorted_SRR15299556.bam

#Conversion of BAM into SAM(optional)
samtools view sorted_SRR15299556.bam > sorted_SRR15299556.sam

#To remove the duplicates
samtools rmdup -sS sorted_SRR15299556.bam rmdup_SRR15299556.bam

##Calling of variations using GATK
#But GATK needs input in picard-tools format

sudo apt-get install picard-tools

#Convert the ref genome in Picard-tool format

picard-tools CreateSequenceDictionary R=genome.fa O=genome.dict

#Prepare the input unique bam file into picard-tools format

picard-tools AddOrReplaceReadGroups I=rmdup_SRR15299556.bam O=picard_output.bam RGLB=DEMO RGPL=illumina RGPU=run RGSM=SRR15299556 SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT

#Call variation using GATK(HaplotypeCaller= Germline mutation , Mutect2= Somatic variation)
samtools faidx genome.fa

java -jar /mnt/c/Users/bhojr/Desktop/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar HaplotypeCaller -R genome.fa -I picard_output.bam -O Gatk_output.vcf

#Variation filteration using SnpSift
cat Gatk_output.vcf | java -jar /mnt/c/Users/bhojr/Desktop/DNASeq_SRR15299556_NGS_ANalysis/DNAseq/snpEff_latest_core/snpEff/SnpSift.jar filter "((QUAL>=30) & (DP>=10) & (MQ>=30))" > Filterd_GATK.vcf
# path for VCF
/mnt/c/Users/bhojr/Desktop/DNASeq_SRR15299556_NGS ANalysis/DNAseq/snpEff_latest_core/snpEff

# VEP( Variant Effect Predictor)> go to the google page then type variant effect predictor > Ensamble page is open > choose the gene name and upload the GATK vcf or snpff VCF and run it then downlod the result in txt then save into EXCEL file


