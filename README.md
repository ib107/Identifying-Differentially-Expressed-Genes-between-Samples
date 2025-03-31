# Identifying-Differentially-Expressed-Genes-between-Samples 

### Introduction
This project focuses on identifying differentially expressed genes (DEGs) across developmental stages of *Saccharomyces cerevisiae* (yeast) during velum development. The analysis utilizes RNA-Seq technology to quantify gene expression changes across different time points, employing rigorous computational and statistical methods to ensure accuracy and reliability.

### Objectives
- Perform RNA-Seq analysis to identify DEGs.
- Implement quality control and preprocessing of raw sequence data.
- Conduct statistical analysis to determine significant gene expression changes.
- Visualize and interpret results using various plots.

### Methods
#### Software and Tools Used:
- **FastQC**: Quality control assessment of raw RNA-Seq data.
- **Trimmomatic**: Adapter trimming and removal of low-quality reads.
- **STAR**: Read alignment to reference genome.
- **SAMtools**: Processing and manipulation of alignment files.
- **Subread (featureCounts)**: Gene expression quantification.
- **edgeR (RStudio)**: Differential expression analysis and statistical modeling.

#### Workflow Overview:
1. **Data Preparation**: Retrieve RNA-Seq data and preprocess FASTQ files.
2. **Quality Control**: Analyze sequencing quality using FastQC and MultiQC.
3. **Adapter Trimming**: Remove adapters and low-quality bases using Trimmomatic.
4. **Sequence Alignment**: Align reads to the reference genome with STAR, followed by BAM file processing using SAMtools.
5. **Gene Quantification**: Use featureCounts to generate a count matrix for differential expression analysis.
6. **Differential Expression Analysis (DEA)**:
   - Filter low-expression genes
   - Normalize count data
   - Estimate dispersion and fit statistical models in edgeR
7. **Visualization**:
   - Volcano Plots: Show significant DEGs based on log fold change and p-values.
   - Venn Diagrams: Display overlap of DEGs between developmental stages.
   - Smear Plots: Illustrate distribution of gene expression changes.

### Results
- **806 genes** were found to be differentially expressed across all pairwise comparisons
- **Top DEGs** were identified with log fold change (logFC) and statistical significance
- **Key findings** highlight core genes involved in biofilm development and potential regulators of yeast velum formation

### Discussion
The results indicate dynamic gene expression changes during velum development, with core regulatory genes identified as potential candidates for further experimental validation. The study highlights:
- The reliability of RNA-Seq for transcriptomic analysis
- The importance of stringent statistical analysis for DEG identification
- Future directions, including functional validation and pathway enrichment analysis

### Repository Structure
```
├── Data/              # Bash Output Feature Count data
├── Scripts/           # Python/R scripts for analysis
├── Results/           # Output files including DEG lists
├── Figures/           # Plots (volcano, smear, venn diagrams)
├── README.md          # Project documentation
├── Assignment2.pdf    # Full assignment write-up
```

### How to Reproduce
1. Clone this repository:  
   ```bash
   git clone https://github.com/yourusername/DGE-analysis.git
   cd DGE-analysis
   ```
2. Run the analysis:
   ```bash
   bash scripts/DGE_analysis.sh
   ```
