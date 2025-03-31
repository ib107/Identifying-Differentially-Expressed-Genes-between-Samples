# Create a directory for Assignment 2
mkdir a2 

# Copy all sequence files from the shared directory to the current working directory
cp /scratch/lukens/Assignment_2_Seqs/* .

# Unzip all FASTQ files
gunzip *.fastq.gz

# Rename paired-end FASTQ files by removing '_1' suffix
for file in *_1.fastq; do
    mv "$file" "${file/_1.fastq/.fastq}"
done

# Load necessary computational modules
module load StdEnv/2023 gcc/12.2.0 r/4.2.2 star/2.7.10a samtools/1.16 fastqc/0.11.9 subread/2.0.3

# Create a directory for quality control (QC) reports
mkdir -p QC_reports

# Allocate computational resources for 2 hours with 32GB memory and 8 CPUs
salloc --time=02:00:00 --mem=32G --cpus-per-task=8

# Run FastQC on all FASTQ files for quality assessment
for file in *.fastq; do 
    echo "Running FastQC on $file..." 
    fastqc -t 8 "$file" -o QC_reports 
done 


# Generate a consolidated QC report using MultiQC
multiqc QC_reports -o QC_reports

# Unzip all FastQC reports for inspection
for file in *_fastqc.zip; do 
    echo "Unzipping $file..."
    unzip -o "$file" -d . 
done 

# Remove zipped FastQC reports to save space
rm -v *_fastqc.zip

# Transfer QC reports to local machine
scp -r ibaxi@graham.computecanada.ca:/home/ibaxi/scratch/genomics/a2/QC_reports .

scp -r ibaxi@graham.computecanada.ca:/home/ibaxi/scratch/genomics/a2/trimmed_fastq/QC_output/multiqc_report.html .

# Organize FASTQ files into separate directories
mkdir unzipped zipped
mv SRR*.fastq unzipped/

# Load Trimmomatic for adapter and quality trimming
module load trimmomatic

# Run Trimmomatic for each FASTQ file to remove low-quality bases
for file in SRR*.fastq; do
    base=$(basename "$file" .fastq)
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 8 \
        "$file" "${base}_trimmed.fastq" HEADCROP:20
    # HEADCROP:20 removes the first 20 bases from each read
    # Single-end (SE) mode is used
    # Uses 8 threads for faster processing
    
done

# Check available STAR aligner versions
module spider star/2.7.9a

# Perform read alignment using STAR for each FASTQ file
for i in *.fastq; do
    j=$(basename "$i")  # Extracts the filename from the path
    STAR --runMode alignReads \
         --runThreadN 8 \
         --genomeDir /scratch/lukens/Assignment_2_Genome \
         --readFilesIn "$i" \
         --outFileNamePrefix "$j"
done

# Create directories to organize output files
mkdir -p fastq_files aligned_sam_files logs_final logs_progress logs_out splice_junctions

# Move files into respective directories
mv SRR*.fastq fastq_files/  # Move FASTQ files
mv SRR*.fastqAligned.out.sam aligned_sam_files/  # Move aligned SAM files
mv SRR*.fastqLog.final.out logs_final/  # Move final log files
mv SRR*.fastqLog.progress.out logs_progress/  # Move progress logs
mv SRR*.fastqLog.out logs_out/  # Move general log files
mv SRR*.fastqSJ.out.tab splice_junctions/  # Move splice junction files

# Convert SAM files to BAM format using samtools
for file in *.sam; do
    samtools view -S -b "$file" > "${file%.sam}.bam"
done

# Sort BAM files for downstream analysis
for file in *.bam; do
    samtools sort "$file" -o "${file%.bam}.sorted.bam"
done

# Load Subread package for feature counting
module load subread

# Perform gene expression quantification using featureCounts
featureCounts -T 8 -a /scratch/lukens/Assignment_2_Genome/genomic.gtf -o gene_counts.txt *.sorted.bam

# Transfer the gene count file to the local machine
scp -r ibaxi@graham.computecanada.ca:/home/ibaxi/scratch/genomics/a2/trimmed_fastq/sam_files/gene_counts.txt .
