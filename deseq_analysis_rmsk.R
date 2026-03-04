library(tidyverse)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(functions)
library(DESeq2)
library(vsn)
library(pheatmap)
library(openxlsx)
library(Biostrings)
library(stringdist)

setwd("/standard/CookLab/Ryan/plasmidsaurus_de/")

######################################Generate count matrix for DESEQ2
#awk 'NR==1 || ($2+$3+$4) > 0' counts.tsv > counts.nonzero.tsv

files <- list.files(
  pattern = ".*alignment_04.*ReadsPerGene.out.expressed.tab.gz$",
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
possible_conditions <- sampleNames
possible_conditions[c(1:3, 14:16)] <- "HCEC"
possible_conditions[c(4:6,10, 17,18)] <- "CRC3_tumoroid"
possible_conditions[c(7:9, 19:21)] <- "CRC3_dmem"
possible_conditions[11:13] <- "CRC4_tumoroid"
possible_conditions[22:24] <- "CRC2_tumoroid"
possible_conditions

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
dds$condition <- relevel(factor(dds$condition), ref = "HCEC")

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
#plotPCA(dds, intgroup = c("batch", "condition"))
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))
meanSdPlot(assay(dds))

notAllZero <- rowSums(counts(dds)) > 0

meanSdPlot(assay(ntd)[notAllZero, ])
meanSdPlot(assay(vsd)[notAllZero, ])

plotDispEsts(dds)

dds2 <- DESeq(dds, fitType="local")   # or "mean"
vsd2 <- varianceStabilizingTransformation(dds2, blind=FALSE)
meanSdPlot(assay(vsd2)[notAllZero, ])
plotDispEsts(dds2)

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
erv_library <- readDNAStringSet("/standard/CookLab/Ryan/r_projects/HERV_library_work/zenodo_17662456/HERV_internal_sequences_v5.fasta")
erv_library <- names(erv_library)
erv_library <- rownames(norm_counts)[which(rownames(norm_counts) %in% erv_library)]
View(norm_counts[which(rownames(norm_counts) %in% erv_library),])
grep("HERV", rownames(norm_counts), ignore.case = TRUE)
look_for_herv <- rownames(norm_counts)[!grepl("ENSG", rownames(norm_counts))]

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
which(colnames(df_long) %in% erv_library == TRUE)
grep(erv_library[1], colnames(df_long))

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
df_long_expressed_min_threshold <- df_long_expressed[which(df_long_expressed$norm_plot > 0),]
df_long_expressed_min_threshold <- df_long_expressed[which(df_long_expressed$gene %in% unique(df_long_expressed_min_threshold$gene)),]


#png("batch_by_condition_rmsk.png", width = 16, height = 21, units = "in", res = 300)
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
#dev.off()

res <- results(dds)
plotCounts(dds, gene=which.min(res$padj), intgroup = "condition")
plotMA(res)

##################################
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
#png("batch_by_condition_rmsk_sig.png", width = 64, height = 21, units = "in", res = 300)
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
#dev.off()

df_long_expressed_sig$locus <- ""
for(i_row in 1:nrow(df_long_expressed_sig)){
  df_long_expressed_sig[i_row, "locus"] <- substr(df_long_expressed_sig[i_row,"gene"][[1]], 1, as.numeric(str_locate_all(df_long_expressed_sig[i_row, "gene"][[1]], "__")[[1]][2,1])-1)
}
unique(df_long_expressed_sig$locus)
unique(df_long_expressed_sig$gene)
#write.xlsx(df_long_expressed_sig, "df_long_expressed_sig.xlsx")
#write.table(df_long_expressed_sig$gene, "significant_hervs.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)

wide_df <- norm_counts[which(rownames(norm_counts) %in% unique(df_long_expressed_sig$gene)),]
wide_df <- cbind("gene" = rownames(wide_df),
                 wide_df)
#write.xlsx(wide_df, "df_wide_expressed_sig.xlsx")

sig_fa <- fread("significant_reads.fa", header = FALSE)
sig_fa_longform <- data.frame("read_name" = sig_fa[seq(1, nrow(sig_fa),2 ),],
                              "nt_seq" = sig_fa[seq(2, nrow(sig_fa),2 ),])
sig_fa_longform$V1 <- substr(sig_fa_longform$V1,2, nchar(sig_fa_longform$V1))
colnames(sig_fa_longform) <- c("read_name", "nt_seq")
sig_fa_longform$aa_seq <- as.character(translate(DNAStringSet(sig_fa_longform$nt_seq)))  
unique(substr(sig_fa_longform$aa_seq, 1, 1))
string_dist <- stringdistmatrix(sig_fa_longform$aa_seq,sig_fa_longform$aa_seq)
rownames(string_dist) <- sig_fa_longform$aa_seq
colnames(string_dist) <- sig_fa_longform$aa_seq
x_chr <- as.character(sig_fa_longform$aa_seq)
n <- length(x_chr)

D <- as.matrix(stringdistmatrix(x_chr, method = "osa"))
L <- outer(nchar(x_chr), nchar(x_chr), pmax)  # max length per pair

D_norm <- D / L
#rownames(D_norm) <- sig_fa_longform$aa_seq
#colnames(D_norm) <- sig_fa_longform$aa_seq
#diag(D_norm) <- 0
string_dist_long <- as.data.frame(D_norm) |>
  tibble::rownames_to_column("string1") |>
  pivot_longer(
    -string1,
    names_to  = "string2",
    values_to = "distance"
  )

#png("string_dist_long_length_normalized_without_29mers.png", width = 64, height = 21, units = "in", res = 300)
print(
ggplot(string_dist_long, aes(x = string1, y = string2, fill = distance)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "String distance"
  )
)
#dev.off()

View(string_dist_long[which(string_dist_long$string1 == "209" & string_dist_long$string2 == "101"),])
sig_fa_longform$aa_seq[209]
sig_fa_longform$aa_seq[101]
low_values <- string_dist_long[which(string_dist_long$distance <0.15),]
low_values$seq1 <- ""
low_values$seq2 <- ""
low_values$nchar1 <- 0
low_values$nchar2 <- 0
for(i_row in 1:nrow(low_values)){
  low_values[i_row, "seq1"] <- sig_fa_longform$aa_seq[as.numeric(low_values[i_row, "string1"])]
  low_values[i_row, "seq2"] <- sig_fa_longform$aa_seq[as.numeric(low_values[i_row, "string2"])]
  low_values[i_row, "nchar1"] <- nchar(sig_fa_longform$aa_seq[as.numeric(low_values[i_row, "string1"])])
  low_values[i_row, "nchar2"] <- nchar(sig_fa_longform$aa_seq[as.numeric(low_values[i_row, "string2"])])
}
low_values$nchar_diff <- abs(low_values$nchar1 - low_values$nchar2)
low_values <- low_values[which(low_values$nchar_diff > 0),]

string_dist_long$seq1 <- ""
string_dist_long$seq2 <- ""
string_dist_long$nchar1 <- 0
string_dist_long$nchar2 <- 0
for(i_row in 1:nrow(string_dist_long)){
  string_dist_long[i_row, "seq1"] <- sig_fa_longform$aa_seq[as.numeric(string_dist_long[i_row, "string1"])]
  string_dist_long[i_row, "seq2"] <- sig_fa_longform$aa_seq[as.numeric(string_dist_long[i_row, "string2"])]
  string_dist_long[i_row, "nchar1"] <- nchar(sig_fa_longform$aa_seq[as.numeric(string_dist_long[i_row, "string1"])])
  string_dist_long[i_row, "nchar2"] <- nchar(sig_fa_longform$aa_seq[as.numeric(string_dist_long[i_row, "string2"])])
  print(i_row)
}
string_dist_long$nchar_diff <- abs(string_dist_long$nchar1 - string_dist_long$nchar2)
string_dist_long <- string_dist_long[which(string_dist_long$nchar_diff > 0),]

interesting_cases <- string_dist_long[which(string_dist_long$nchar_diff > 32 &
                                              string_dist_long$distance < 0.6),]


chunk_string <- function(x, chunk_length) {
  stopifnot(is.character(x))
  stopifnot(length(chunk_length) == 1, is.numeric(chunk_length), chunk_length > 0)
  
  lapply(x, function(s) {
    chars <- strsplit(s, "", fixed = TRUE)[[1]]
    n <- length(chars)
    
    # If shorter than or equal to chunk length, return as-is
    if (n <= chunk_length) return(s)
    
    # Group character positions into chunk bins: 1..k
    grp <- ceiling(seq_len(n) / chunk_length)
    
    # Paste each group back into a chunk string
    vapply(split(chars, grp), paste0, character(1), collapse = "")
  })
}
class(chunk_string("ABCDE", 3)[[1]])
class(chunk_string("ABC", 3)[[1]])
chunk_string("ABCDEF", 3)
unlist(chunk_string(c("ABC", "DEF"), 3))
chunk_string(c("ABCDE", "DEF"), 3)[[1]][2]

chunked_aa <- unlist(chunk_string(sig_fa_longform$aa_seq, 29))
summary(nchar(chunked_aa))
chunked_aa <- chunked_aa[which(nchar(chunked_aa) >= 9)]

string_dist <- stringdistmatrix(chunked_aa,chunked_aa)
rownames(string_dist) <- 1:367
colnames(string_dist) <- 1:367

string_dist_long <- as.data.frame(string_dist) |>
  tibble::rownames_to_column("string1") |>
  pivot_longer(
    -string1,
    names_to  = "string2",
    values_to = "distance"
  )

png("string_dist_chunked_29mers.png", width = 25, height = 21, units = "in", res = 300)
print(
  ggplot(string_dist_long, aes(x = string1, y = string2, fill = distance)) +
    geom_tile() +
    #scale_fill_viridis_c() +
    coord_fixed() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = "String distance"
    )
)
dev.off()

write.xlsx(sig_fa_longform, "sig_antigens.xlsx")
