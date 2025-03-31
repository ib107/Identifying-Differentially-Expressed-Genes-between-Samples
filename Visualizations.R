# Load required libraries
install.packages("ggplot2")
install.packages("VennDiagram")
install.packages("EnhancedVolcano")
BiocManager::install("edgeR")

library(ggplot2)
library(VennDiagram)
library(EnhancedVolcano)
library(edgeR)

library(ggplot2)
library(dplyr)
library(gridExtra)

# Load count data
gene_counts <- read.delim("featureGenecount.txt", comment.char="#", check.names=FALSE, stringsAsFactors=FALSE)

# Define sample groups
Stage <- factor(rep(c("Early biofilm", "Thin biofilm", "Mature biofilm"), each=3))
y <- DGEList(counts=gene_counts[,7:ncol(gene_counts)], group=Stage, genes=gene_counts[,1:6])

# Filter and normalize
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")

# Create design matrix
designmatrix <- model.matrix(~0 + Stage)
rownames(designmatrix) <- colnames(y)

# Estimate dispersion
y <- estimateDisp(y, designmatrix, robust=TRUE)

# Fit the model
fit <- glmQLFit(y, designmatrix)

# Differential Expression Analysis
deg_early_vs_thin <- glmQLFTest(fit, contrast = c(1, -1, 0))  
deg_thin_vs_mature <- glmQLFTest(fit, contrast = c(0, 1, -1))  
deg_early_vs_mature <- glmQLFTest(fit, contrast = c(1, 0, -1)) 

# Extract top DEGs
top_early_vs_thin <- topTags(deg_early_vs_thin, n = 10, p.value = 0.05)$table
top_thin_vs_mature <- topTags(deg_thin_vs_mature, n = 10, p.value = 0.05)$table
top_early_vs_mature <- topTags(deg_early_vs_mature, n = 10, p.value = 0.05)$table

# Add regulation status (Upregulated/Downregulated) based on logFC
top_early_vs_thin$Regulation <- ifelse(top_early_vs_thin$logFC > 0, "Upregulated", "Downregulated")
top_thin_vs_mature$Regulation <- ifelse(top_thin_vs_mature$logFC > 0, "Upregulated", "Downregulated")
top_early_vs_mature$Regulation <- ifelse(top_early_vs_mature$logFC > 0, "Upregulated", "Downregulated")

# Combine results
list_anytwotreatments <- data.frame(
  Gene = rownames(top_early_vs_thin),
  Regulation_Early_vs_Thin = top_early_vs_thin$Regulation,
  logFC_Early_vs_Thin = top_early_vs_thin$logFC,
  PValue_Early_vs_Thin = top_early_vs_thin$PValue,
  Regulation_Thin_vs_Mature = top_thin_vs_mature$Regulation,
  logFC_Thin_vs_Mature = top_thin_vs_mature$logFC,
  PValue_Thin_vs_Mature = top_thin_vs_mature$PValue,
  Regulation_Early_vs_Mature = top_early_vs_mature$Regulation,
  logFC_Early_vs_Mature = top_early_vs_mature$logFC,
  PValue_Early_vs_Mature = top_early_vs_mature$PValue
)

# Early vs Later (Thin & Mature)
early_vs_later <- c(1, -0.5, -0.5) 
deg_early_vs_later <- glmQLFTest(fit, contrast=early_vs_later) 

# Extract DEGs
top_early_vs_later <- topTags(deg_early_vs_later, n=10, p.value=0.05)$table

# Add regulation status (Upregulated/Downregulated) based on logFC
top_early_vs_later$Regulation <- ifelse(top_early_vs_later$logFC > 0, "Upregulated", "Downregulated")

# Combine results
list_early_vs_later <- data.frame(
  Gene = rownames(top_early_vs_later),
  Regulation_Early_vs_Later = top_early_vs_later$Regulation,
  logFC_Early_vs_Later = top_early_vs_later$logFC,
  PValue_Early_vs_Later = top_early_vs_later$PValue, 
  FDRValue_Early_vs_Later = top_early_vs_later$FDR
)

# Remove unnecessary columns (chr, start, end, strand, length) from the top 10 lists
#list_anytwotreatments <- list_anytwotreatments %>%
  select(-matches("chr|start|end|strand|length"))

list_early_vs_later <- list_early_vs_later %>%
  select(-matches("chr|start|end|strand|length"))

top_early_vs_mature <- top_early_vs_mature %>%
  select(-matches("chr|start|end|strand|length"))

top_early_vs_thin <- top_early_vs_thin %>%
  select(-matches("chr|start|end|strand|length"))

top_thin_vs_mature <- top_thin_vs_mature %>%
  select(-matches("chr|start|end|strand|length"))

# Save results as text files
write.table(list_anytwotreatments, "DEGs_Pairwise.txt", row.names=FALSE, quote=FALSE, sep="\t")
write.table(list_early_vs_later, "DEGs_Early_vs_Later.txt", row.names=FALSE, quote=FALSE, sep="\t")

### Visualization 1 ----

# Extract all DEGs for each pairwise comparison
all_early_vs_thin <- topTags(deg_early_vs_thin, n = Inf, p.value = 0.05)$table
all_thin_vs_mature <- topTags(deg_thin_vs_mature, n = Inf, p.value = 0.05)$table
all_early_vs_mature <- topTags(deg_early_vs_mature, n = Inf, p.value = 0.05)$table

# Create a Venn diagram for all DEGs
venn_data_all <- list(
  'Early vs Thin' = rownames(all_early_vs_thin),
  'Thin vs Mature' = rownames(all_thin_vs_mature),
  'Early vs Mature' = rownames(all_early_vs_mature)
)

venn.diagram(venn_data_all,
             filename = 'venn_diagram_all_genes.png',
             category.names = c("Early vs Thin", "Thin vs Mature", "Early vs Mature"),
             output = TRUE)

### Visualization #2 ----
# Smear Plot
dt_early_vs_later <- decideTestsDGE(deg_early_vs_later)
plotSmear(deg_early_vs_later, de.tags = rownames(y)[as.logical(dt_early_vs_later)], main = "Smear Plot for Early vs Later Comparison")
abline(h = c(-1, 1), col = "red")

# Visualization #3 -------

# Combine all pairwise results into a single data frame for volcano plots
volcano_data <- data.frame(
  Gene = rownames(y$genes),
  logFC_Early_vs_Thin = deg_early_vs_thin$table$logFC,
  PValue_Early_vs_Thin = deg_early_vs_thin$table$PValue,
  logFC_Thin_vs_Mature = deg_thin_vs_mature$table$logFC,
  PValue_Thin_vs_Mature = deg_thin_vs_mature$table$PValue,
  logFC_Early_vs_Mature = deg_early_vs_mature$table$logFC,
  PValue_Early_vs_Mature = deg_early_vs_mature$table$PValue
)

# Define thresholds for significance
p_threshold <- 0.05
fc_threshold <- 1

# Add columns to classify genes as upregulated, downregulated, or not significant for each comparison
volcano_data <- volcano_data %>%
  mutate(
    Significance_Early_vs_Thin = case_when(
      PValue_Early_vs_Thin < p_threshold & abs(logFC_Early_vs_Thin) > fc_threshold ~ ifelse(logFC_Early_vs_Thin > 0, "Upregulated", "Downregulated"),
      TRUE ~ "Not Significant"
    ),
    Significance_Thin_vs_Mature = case_when(
      PValue_Thin_vs_Mature < p_threshold & abs(logFC_Thin_vs_Mature) > fc_threshold ~ ifelse(logFC_Thin_vs_Mature > 0, "Upregulated", "Downregulated"),
      TRUE ~ "Not Significant"
    ),
    Significance_Early_vs_Mature = case_when(
      PValue_Early_vs_Mature < p_threshold & abs(logFC_Early_vs_Mature) > fc_threshold ~ ifelse(logFC_Early_vs_Mature > 0, "Upregulated", "Downregulated"),
      TRUE ~ "Not Significant"
    )
  )

# Volcano Plot for Early vs Thin
volcano_plot_Early_vs_Thin <- ggplot(volcano_data, aes(x = logFC_Early_vs_Thin, y = -log10(PValue_Early_vs_Thin), color = Significance_Early_vs_Thin)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(fc_threshold, -fc_threshold), linetype = "dashed", color = "black") +
  labs(
    title = "Early vs Thin",
    x = "Log2 Fold Change (Early vs Thin)",
    y = "-Log10(p-value)",
    color = "Significance"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

# Volcano Plot for Thin vs Mature
volcano_plot_Thin_vs_Mature <- ggplot(volcano_data, aes(x = logFC_Thin_vs_Mature, y = -log10(PValue_Thin_vs_Mature), color = Significance_Thin_vs_Mature)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(fc_threshold, -fc_threshold), linetype = "dashed", color = "black") +
  labs(
    title = "Thin vs Mature",
    x = "Log2 Fold Change (Thin vs Mature)",
    y = "-Log10(p-value)",
    color = "Significance"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

# Volcano Plot for Early vs Mature
volcano_plot_Early_vs_Mature <- ggplot(volcano_data, aes(x = logFC_Early_vs_Mature, y = -log10(PValue_Early_vs_Mature), color = Significance_Early_vs_Mature)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(fc_threshold, -fc_threshold), linetype = "dashed", color = "black") +
  labs(
    title = "Early vs Mature",
    x = "Log2 Fold Change (Early vs Mature)",
    y = "-Log10(p-value)",
    color = "Significance"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

# Combine the plots side by side
grid.arrange(volcano_plot_Early_vs_Thin, volcano_plot_Thin_vs_Mature, volcano_plot_Early_vs_Mature, ncol = 3)
