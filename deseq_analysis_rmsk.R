library(tidyverse)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(functions)
library(DESeq2)

setwd("/scratch/rms2jg/plasmidsaurus_de/")

######################################Generate count matrix for DESEQ2
files <- list.files(
  pattern = "M7.*rmsk.*ReadsPerGene.expressed.out.tab$",
  full.names = TRUE
)
sampleNames <- sub("ReadsPerGene.expressed.out.tab$", "", basename(files))

#strandedness
use_col <- 2

read_star_counts <- function(fn, col = 2) {
  tab <- read.delim(fn, header = FALSE, stringsAsFactors = FALSE)
  # V1 = gene ID, V[col] = counts
  counts <- tab[, c(1, col)]
  colnames(counts) <- c("gene_id", "count")
  counts
}

counts_list <- lapply(files, read_star_counts, col = use_col)

genes_master <- data.frame(gene_id = sort(unique(unlist(lapply(counts_list, `[[`, "gene_id")))))
for(i_index in 1:length(counts_list)){
  temp_df <- merge(genes_master, counts_list[[i_index]], by = "gene_id", all = TRUE)
  temp_df[is.na(temp_df$count), "count"] <- 0
  counts_list[[i_index]] <- temp_df
  rm(temp_df)
}

#sanity check that gene_ids are consistent across samples
stopifnot(all(sapply(counts_list, function(x) identical(counts_list[[1]]$gene_id, x$gene_id))))

gene_ids <- counts_list[[1]]$gene_id
cts <- do.call(
  cbind,
  lapply(counts_list, function(x) x$count)
)

rownames(cts) <- gene_ids
colnames(cts) <- sampleNames
 #
#########################################Create DESeq2 colData
possible_conditions <- gsub("^M7L76V_[0-9]{1,2}_sample_[0-9]{1,2}_", "", sampleNames)
possible_conditions <- gsub("_rep[0-9]{1}_RAWSEQ", "", possible_conditions)
sampleTable <- data.frame(
  row.names = sampleNames,
  condition = possible_conditions
)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = sampleTable,
                              design = ~condition)

dds <- DESeq(dds)  

# #HCEC vs tumoroid
# res_hcec_tumoroid <- results(dds, contrast = c("condition", "HCEC", "tumoroid"))
# #log2FoldChange > 0 means genes higher in HCEC than tumoroid
# #log2FoldChange < 0 means genes higher in tumoroid than HCEC
# 
# #HCEC vs DMEM
# res_hcec_dmem <- results(dds, contrast = c("condition", "HCEC", "DMEM"))  
# #log2FoldChange > 0 means genes higher in HCEC than DMEM
# #log2FoldChange < 0 means genes higher in DMEM than HCEC
# 
# #tumoroid vs DMEM
# res_tumoroid_dmem <- results(dds, contrast = c("condition", "tumoroid", "DMEM"))  
# #log2FoldChange > 0 means genes higher in tumoroid than DMEM
# #log2FoldChange < 0 means genes higher in DMEM than tumoroid

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
plotPCA(vsd, intgroup = "condition")

# plotMA(res_hcec_tumoroid)
# plotMA(res_hcec_dmem)
# plotMA(res_tumoroid_dmem)

#plotCounts(dds, gene = "XP_055153856.1", intgroup = "condition")

#erv_library <- readRDS("/scratch/rms2jg/gene_names.RDS")
#housekeepers <- c("ENSG00000089157", "ENSG00000166710", "ENSG00000075624") #RPLP0 (ribosomal protein), B2M, ACTB
#present_hk <- intersect(housekeepers, rownames(dds))

norm_counts <- counts(dds, normalized = TRUE)
# norm_sub <- norm_counts[which(rownames(norm_counts) %in% c(present_hk, erv_library)),]
# erv_library <- rownames(norm_sub)
erv_library <- rownames(norm_counts)

df <- as.data.frame(t(norm_counts))
df$sample <- rownames(df)
meta <- as.data.frame(colData(dds))
meta$sample <- rownames(meta)

df_long <- df %>%
  pivot_longer(
    cols = all_of(erv_library),
    names_to = "gene",
    values_to = "norm_count"
  ) %>%
  left_join(meta, by = "sample") %>%
  mutate(
    norm_count = replace_na(norm_count, 0)
  )

ggplot(df_long, aes(x = gene, y = norm_count, color = condition)) +
  theme_bw()+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  geom_point(size = 3, position = position_jitter(width = 0.1)) +
  scale_y_log10() +  # optional but recommended
  labs(
    x = "Gene",
    y = "Normalized counts (log10)",
    title = "Expression of selected genes across samples"
  )




df_long_expressed <- df_long[which(df_long$norm_count > 0),]
df_long_expressed <- df_long[which(df_long$gene %in% unique(df_long_expressed$gene)),]
# df_long_expressed[which(df_long_expressed$gene == "ENSG00000089157"), "gene"] <- "RPLP0"
# df_long_expressed[which(df_long_expressed$gene == "ENSG00000166710"), "gene"] <- "B2M"
# df_long_expressed[which(df_long_expressed$gene == "ENSG00000075624"), "gene"] <- "ACTB"
# df_long_expressed$gene <- factor(df_long_expressed$gene, levels = c("ACTB", "B2M", "RPLP0", unique(df_long_expressed$gene)[4:length(unique(df_long_expressed$gene))]))
df_long_expressed <- df_long_expressed %>%
  mutate(norm_plot = norm_count + 1) 
dodge <-position_dodge(width = 0.8)
summary(df_long_expressed$norm_plot)
df_long_expressed_min_threshold <- df_long_expressed[which(df_long_expressed$norm_plot > 321),]
df_long_expressed_min_threshold <- df_long_expressed[which(df_long_expressed$gene %in% unique(df_long_expressed_min_threshold$gene)),]


png("M7L76V_rmsk.png", width = 16, height = 21, units = "in", res = 300)
print(
  ggplot(df_long_expressed_min_threshold, aes(x = gene, y = norm_plot, fill = condition)) +
    theme_bw()+
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
          #axis.text.x = element_blank()
          )+
    stat_summary(fun = mean, geom = "col", position = dodge) +
    geom_point(aes(fill = condition), shape = 21, color = "black", position = dodge)+
    scale_y_log10() +
    labs(
      x = "Gene",
      y = "Normalized counts (log10)",
      title = "Expression of ERV genes (M7L76V)"
    )
)
dev.off()

res <- results(dds)
plotCounts(dds, gene=which.min(res$padj), intgroup = "condition")
