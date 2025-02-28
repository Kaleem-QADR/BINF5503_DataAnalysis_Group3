# DAS Assignment 2 - RNA-Seq Differential Expression Analysis

# Load Required Libraries
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

# ============================
# Step 1: Load and Wrangle Data
# ============================

# Define directory path containing the dataset
read_counts_dir <- "1_CancerTranscriptomics/read_counts"

# Define file paths
file_paths <- list.files(path = read_counts_dir, pattern = "*.tsv", full.names = TRUE)

# Extract sample names from file names
file_names <- list.files(path = read_counts_dir, pattern = "*.tsv", full.names = FALSE)
sample_names <- tools::file_path_sans_ext(file_names)

# Load datasets and assign sample names
data_list <- setNames(lapply(seq_along(file_paths), function(i) {
  read_tsv(file_paths[i], col_names = c("Ensembl_ID", "Gene_Name", "Read_Counts")) %>%
    mutate(Sample = sample_names[i])
}), sample_names)

# Merge datasets into a single dataframe
merged_data <- bind_rows(data_list)

# Check loaded data
print(dim(merged_data))  # Rows x Columns
print(head(merged_data))  # Preview first few rows
print(unique(merged_data$Sample))  # Ensure sample names are correct

# ============================
# Step 2: Prepare Data for DESeq2
# ============================

# Create sample information metadata
sample_info <- data.frame(
  Sample = unique(merged_data$Sample),
  condition = ifelse(grepl("_D0", unique(merged_data$Sample)), "D0", "D6")
)
rownames(sample_info) <- sample_info$Sample

# Transform data into count matrix
count_matrix <- merged_data %>%
  select(Ensembl_ID, Sample, Read_Counts) %>%
  pivot_wider(names_from = Sample, values_from = Read_Counts, values_fill = list(Read_Counts = 0)) %>%
  column_to_rownames(var = "Ensembl_ID")

# Convert to numeric matrix
count_matrix <- as.matrix(count_matrix)

# Check matrix format
print(dim(count_matrix))  # Genes x Samples
print(head(rownames(count_matrix), 10))  # Should return Ensembl IDs

# ============================
# Step 3: Differential Expression Analysis
# ============================

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_info, design = ~ condition)

# Set reference condition
dds$condition <- relevel(dds$condition, ref = "D0")

# Run DESeq2 normalization & analysis
dds_out <- DESeq(dds)
res <- results(dds_out, alpha = 0.05)
summary(res)

# Convert results to a dataframe
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Remove NA values
res_df <- res_df %>% drop_na(log2FoldChange, padj)

# ============================
# Step 4: Visualization
# ============================

## Volcano Plot
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

## Heatmap

### Heatmap

# Extract top 30 significant genes (removing NAs first)
filtered_res_df <- res_df %>% drop_na(padj)
top_genes <- rownames(filtered_res_df[order(filtered_res_df$padj), ])[1:min(30, nrow(filtered_res_df))]

# Ensure the selected genes exist in the dataset
if (exists("dds_out")) {
  vst_data <- vst(dds_out, blind = FALSE)  # Precompute VST data to avoid redundant calls
  top_genes <- top_genes[!is.na(top_genes) & top_genes %in% rownames(vst_data)]
}

if (length(top_genes) > 0) {
  norm_counts <- assay(vst_data)[top_genes, ]
  
  # Replace Ensembl ID with Gene Name, fallback to ID if NA
  gene_names <- merged_data$Gene_Name[match(top_genes, merged_data$Ensembl_ID)]
  rownames(norm_counts) <- ifelse(is.na(gene_names), top_genes, gene_names)
  
  # Fix sample names by removing "_gene"
  colnames(norm_counts) <- gsub("_gene", "", colnames(norm_counts))  
  
  pheatmap(norm_counts,
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           main = "Heatmap of Top 30 DEGs",
           fontsize_row = 8,   
           fontsize_col = 10,   
           angle_col = 45)
} else {
  print("No significant genes found for heatmap.")
}

# ============================
# Step 5: Table of Top DEGs
# ============================

# Select top 10 upregulated and downregulated genes using head()
top_up <- res_df %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange)) %>% head(5)
top_down <- res_df %>% filter(padj < 0.05) %>% arrange(log2FoldChange) %>% head(5)
top_genes_table <- bind_rows(top_up, top_down) %>%
  select(gene, log2FoldChange, padj)

# Replace Ensembl IDs with Gene Names, ensuring no NA mismatches
top_genes_table$gene <- merged_data$Gene_Name[match(top_genes_table$gene, merged_data$Ensembl_ID)]
top_genes_table$gene[is.na(top_genes_table$gene)] <- top_genes_table$gene[is.na(top_genes_table$gene)]

# Print table
print(top_genes_table)

# Save as CSV
write.csv(top_genes_table, "Top_10_DEGs.csv", row.names = FALSE)
