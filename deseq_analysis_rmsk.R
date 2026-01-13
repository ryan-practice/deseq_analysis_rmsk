library(tidyverse)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(functions)
library(DESeq2)
library(vsn)
library(pheatmap)

setwd("/scratch/rms2jg/plasmidsaurus_de/")

######################################Generate count matrix for DESEQ2
#awk 'NR==1 || ($2+$3+$4) > 0' counts.tsv > counts.nonzero.tsv

files <- list.files(
  pattern = ".*rmsk.*ReadsPerGene.expressed.out.tab$",
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

#get all possible gene names and ensure each sample has them in the same order for their count file
genes_master <- data.frame(gene_id = sort(unique(unlist(lapply(counts_list, `[[`, "gene_id")))))
for(i_index in 1:length(counts_list)){
  temp_df <- merge(genes_master, counts_list[[i_index]], by = "gene_id", all = TRUE)
  temp_df[is.na(temp_df$count), "count"] <- 0
  temp_df <- temp_df[-which(temp_df$gene_id == "N_unmapped"),]
  temp_df <- temp_df[-which(temp_df$gene_id == "N_multimapping"),]
  temp_df <- temp_df[-which(temp_df$gene_id == "N_noFeature"),]
  counts_list[[i_index]] <- temp_df
  rm(temp_df)
}

#sanity check that gene_ids are consistent across samples
stopifnot(all(sapply(counts_list, function(x) identical(counts_list[[1]]$gene_id, x$gene_id))))

#create count matrix where columns are sample and rows are genes
gene_ids <- counts_list[[1]]$gene_id
cts <- do.call(
  cbind,
  lapply(counts_list, function(x) x$count)
)

rownames(cts) <- gene_ids
colnames(cts) <- sampleNames

#########################################Create DESeq2 colData
#any factor used in the deseq formula should be a column in this table. From documentation:
#"The two factor variables batch and condition should be columns of coldata."
#Here I intend to combine both RNAseqs performed by plasmidsaurus and correct for batch effects
possible_conditions <- gsub("^.*_[0-9]{1,2}_sample_[0-9]{1,2}_", "", sampleNames)
possible_conditions <- gsub("_rep[0-9]{1}_RAWSEQ", "", possible_conditions)
possible_conditions <- gsub("_rep_[0-9]{1}", "", possible_conditions)
possible_conditions[grep(".*CRC.*3.*Tumoroid.*", possible_conditions)] <- "CRC_3_Tumoroid_rmsk_herv"
possible_conditions[grep(".*CRC.*3.*DMEM.*", possible_conditions)] <- "CRC_3_DMEM_rmsk_herv"
possible_conditions[grep(".*CRC.*4.*Tumoroid.*", possible_conditions)] <- "CRC_4_Tumoroid_rmsk_herv"

sampleTable <- data.frame(
  row.names = sampleNames,
  condition = possible_conditions
)

sampleTable$batch <- ""
sampleTable[which(substr(rownames(sampleTable), 1, 6) == "9FC5BD"), "batch"] <- "1"
sampleTable[which(substr(rownames(sampleTable), 1, 6) == "M7L76V"), "batch"] <- "2"


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = sampleTable,
                              design = ~ batch + condition)

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

ntd <- normTransform(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
plotPCA(ntd, intgroup = c("batch", "condition"))
plotPCA(vsd, intgroup = c("batch", "condition"))
plotPCA(rld, intgroup = c("batch", "condition"))
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

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



#This step is redundant now that I filter in bash before reading into R.
#Right now, I filter > 0 but in the DESeq2 documentation, it suggests 10. May be worth considering but for now I want to be permissive
df_long_expressed <- df_long[which(df_long$norm_count > 0),]
df_long_expressed <- df_long[which(df_long$gene %in% unique(df_long_expressed$gene)),]

#Keeping this code for future when I have the human genome included I can re-include housekeeping
# df_long_expressed[which(df_long_expressed$gene == "ENSG00000089157"), "gene"] <- "RPLP0"
# df_long_expressed[which(df_long_expressed$gene == "ENSG00000166710"), "gene"] <- "B2M"
# df_long_expressed[which(df_long_expressed$gene == "ENSG00000075624"), "gene"] <- "ACTB"
# df_long_expressed$gene <- factor(df_long_expressed$gene, levels = c("ACTB", "B2M", "RPLP0", unique(df_long_expressed$gene)[4:length(unique(df_long_expressed$gene))]))

df_long_expressed <- df_long_expressed %>%
  mutate(norm_plot = norm_count + 1) 
dodge <-position_dodge(width = 0.8)
summary(df_long_expressed$norm_plot)
df_long_expressed_min_threshold <- df_long_expressed[which(df_long_expressed$norm_plot > 200),]
df_long_expressed_min_threshold <- df_long_expressed[which(df_long_expressed$gene %in% unique(df_long_expressed_min_threshold$gene)),]


png("batch_by_condition_rmsk.png", width = 16, height = 21, units = "in", res = 300)
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
plotMA(res)

##################################
dds$condition <- relevel(factor(dds$condition), ref = "HCEC_rmsk_herv")
design(dds)
dds <- DESeq(dds)

#see documentation of DESeq for how reduced works. "the full formula with the term(s) of interest removed"
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ batch) 
res_lrt <- results(dds_lrt)
res_lrt["LTR54B__chr19__41771967-41772234__plus(+)_13", ]
res_lrt_df <- as.data.frame(res_lrt)
res_lrt_df$gene_name <- rownames(res_lrt_df)
res_lrt_df_sig <- res_lrt_df[which(res_lrt_df$pvalue <= 0.05),]
v_plotCounts <- Vectorize("plotCounts", vectorize.args = "gene")
v_plotCounts(dds, gene=res_lrt_df_sig$gene_name, intgroup = "condition")

df_long_expressed_sig <- df_long_expressed[which(df_long_expressed$gene %in% res_lrt_df_sig$gene_name),]
png("batch_by_condition_rmsk_sig.png", width = 64, height = 21, units = "in", res = 300)
print(
  ggplot(df_long_expressed_sig, aes(x = gene, y = norm_plot, fill = condition)) +
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