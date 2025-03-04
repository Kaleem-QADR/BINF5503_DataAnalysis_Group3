# DAS Assignment 2 - RNA-Seq Differential Expression Analysis

# Load Required Libraries
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)


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
  dplyr::select(Ensembl_ID, Sample, Read_Counts) %>%
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
top_up <- res_df %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange)) %>% head(10)
top_down <- res_df %>% filter(padj < 0.05) %>% arrange(log2FoldChange) %>% head(10)
top_genes_table <- bind_rows(top_up, top_down) %>%
  dplyr::select(gene, log2FoldChange, padj)

# Replace Ensembl IDs with Gene Names, ensuring no NA mismatches
top_genes_table$gene <- merged_data$Gene_Name[match(top_genes_table$gene, merged_data$Ensembl_ID)]
top_genes_table$gene[is.na(top_genes_table$gene)] <- top_genes_table$gene[is.na(top_genes_table$gene)]

# Print table
print(top_genes_table)

# Save as CSV
write.csv(top_genes_table, "Top_20_DEGs.csv", row.names = FALSE)


# ============================
# Step 6: GO/KEGG Enrichment for Pathway Identification
# ============================

# Remove version numbers from Ensembl IDs
top_genes <- unique(res_df %>% filter(padj < 0.05) %>% pull(gene))
top_genes_clean <- gsub("\\..*", "", top_genes)  

# Convert Ensembl IDs to Entrez IDs
gene_ids <- mapIds(org.Hs.eg.db, 
                   keys = top_genes_clean, 
                   column = "ENTREZID", 
                   keytype = "ENSEMBL", 
                   multiVals = "first")

# If mapping fails, try biomaRt
if (sum(is.na(gene_ids)) > length(gene_ids) * 0.5) { 
  message("⚠ Many genes failed to map. Using biomaRt...")
  
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  converted <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                     filters = "ensembl_gene_id",
                     values = top_genes_clean,
                     mart = mart)
  
  gene_ids <- setNames(converted$entrezgene_id, converted$ensembl_gene_id)
}

# Remove NA values
gene_ids <- na.omit(gene_ids)

# Perform GO enrichment
go_results <- enrichGO(
  gene          = gene_ids, 
  OrgDb         = org.Hs.eg.db, 
  keyType       = "ENTREZID", 
  ont           = "BP",  # Biological Process 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# Print summary
print(head(go_results))

# Save results
write.csv(as.data.frame(go_results), "GO_Enrichment_Results.csv", row.names = FALSE)

# Plot Top 10 GO Terms
if (nrow(as.data.frame(go_results)) > 0) {
  dotplot(go_results, showCategory = 10, title = "Top 10 GO Enriched Biological Processes")
} else {
  message("⚠ No significant GO terms found.")
}

# Perform KEGG pathway enrichment
kegg_results <- enrichKEGG(
  gene          = gene_ids, 
  organism      = "hsa",  # Human
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# Print summary
print(head(kegg_results))

# Save results
write.csv(as.data.frame(kegg_results), "KEGG_Enrichment_Results.csv", row.names = FALSE)

# Plot Top 10 KEGG Pathways
if (nrow(as.data.frame(kegg_results)) > 0) {
  dotplot(kegg_results, showCategory = 10, title = "Top 10 Enriched KEGG Pathways")
} else {
  message("⚠ No significant KEGG pathways found.")
}
