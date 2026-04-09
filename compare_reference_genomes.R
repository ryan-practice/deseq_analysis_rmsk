#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biostrings)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript compare_herv_fastas.R mine.fa theirs.fa [output_prefix]")
}

fa1 <- args[1]
fa2 <- args[2]
prefix <- ifelse(length(args) >= 3, args[3], "herv_fasta_compare")

message("Reading FASTA files...")
x1 <- readDNAStringSet(fa1)
x2 <- readDNAStringSet(fa2)

# ----------------------------
# Helper functions
# ----------------------------

safe_median <- function(x) if (length(x) == 0) NA_real_ else median(x)
safe_mean   <- function(x) if (length(x) == 0) NA_real_ else mean(x)
safe_min    <- function(x) if (length(x) == 0) NA_real_ else min(x)
safe_max    <- function(x) if (length(x) == 0) NA_real_ else max(x)

calc_gc <- function(seq_char) {
  chars <- strsplit(seq_char, "", fixed = TRUE)[[1]]
  if (length(chars) == 0) return(NA_real_)
  chars <- toupper(chars)
  sum(chars %in% c("G", "C")) / length(chars) * 100
}

calc_n_pct <- function(seq_char) {
  chars <- strsplit(seq_char, "", fixed = TRUE)[[1]]
  if (length(chars) == 0) return(NA_real_)
  chars <- toupper(chars)
  sum(chars == "N") / length(chars) * 100
}

shannon_complexity <- function(seq_char) {
  chars <- strsplit(toupper(seq_char), "", fixed = TRUE)[[1]]
  chars <- chars[chars %in% c("A", "C", "G", "T", "N")]
  if (length(chars) == 0) return(NA_real_)
  p <- table(chars) / length(chars)
  -sum(p * log2(p))
}

extract_locus_like <- function(header) {
  # tries to pull locus-ish info from headers like:
  # repName__chr__start-end__strand
  # or chr:start-end
  out <- str_extract(header, "chr[[:alnum:]_]+[:_][0-9]+[-_][0-9]+")
  out[is.na(out)] <- str_extract(header, "chr[[:alnum:]_]+")
  out
}

extract_family_like <- function(header) {
  # crude family extraction: before first "__" or first whitespace
  x <- str_extract(header, "^[^[:space:]]+")
  x <- str_replace(x, "__.*$", "")
  x
}

seq_df <- function(x, label) {
  headers <- names(x)
  seqs <- as.character(x)
  lens <- width(x)
  
  tibble(
    library = label,
    header = headers,
    sequence = seqs,
    length = as.integer(lens),
    gc_percent = vapply(seqs, calc_gc, numeric(1)),
    n_percent = vapply(seqs, calc_n_pct, numeric(1)),
    complexity = vapply(seqs, shannon_complexity, numeric(1)),
    locus_like = extract_locus_like(headers),
    family_like = extract_family_like(headers)
  )
}

count_exact_duplicates <- function(seqs) {
  tab <- table(seqs)
  n_dup_entries <- sum(tab > 1)
  n_dup_records <- sum(tab[tab > 1])
  tibble(
    unique_sequences = length(tab),
    duplicate_sequence_groups = n_dup_entries,
    duplicate_records_total = n_dup_records
  )
}

kmer_set <- function(seqs, k = 21, max_seqs = Inf) {
  out <- new.env(hash = TRUE, parent = emptyenv())
  n <- min(length(seqs), max_seqs)
  
  for (i in seq_len(n)) {
    s <- toupper(seqs[i])
    s <- gsub("[^ACGT]", "N", s)
    if (nchar(s) < k) next
    starts <- seq_len(nchar(s) - k + 1)
    km <- substring(s, starts, starts + k - 1)
    km <- km[!grepl("N", km, fixed = TRUE)]
    for (kk in km) out[[kk]] <- TRUE
  }
  
  ls(out, all.names = TRUE)
}

jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  if (union == 0) return(NA_real_)
  inter / union
}

private_kmer_fraction <- function(a, b) {
  if (length(a) == 0) return(NA_real_)
  sum(!(a %in% b)) / length(a)
}

n50_fun <- function(lengths) {
  if (length(lengths) == 0) return(NA_real_)
  x <- sort(lengths, decreasing = TRUE)
  csum <- cumsum(x)
  target <- sum(x) / 2
  x[which(csum >= target)[1]]
}

# ----------------------------
# Build tables
# ----------------------------

d1 <- seq_df(x1, "mine")
d2 <- seq_df(x2, "theirs")
d <- bind_rows(d1, d2)

# ----------------------------
# Basic summary
# ----------------------------

summary_tbl <- d %>%
  group_by(library) %>%
  summarise(
    n_records = n(),
    total_bp = sum(length),
    min_len = safe_min(length),
    median_len = safe_median(length),
    mean_len = safe_mean(length),
    max_len = safe_max(length),
    N50 = n50_fun(length),
    mean_gc_percent = safe_mean(gc_percent),
    mean_n_percent = safe_mean(n_percent),
    mean_complexity = safe_mean(complexity),
    unique_headers = n_distinct(header),
    unique_sequences = n_distinct(sequence),
    unique_locus_like = n_distinct(locus_like, na.rm = TRUE),
    unique_family_like = n_distinct(family_like, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------
# Duplicate stats
# ----------------------------

dup1 <- count_exact_duplicates(d1$sequence) %>% mutate(library = "mine")
dup2 <- count_exact_duplicates(d2$sequence) %>% mutate(library = "theirs")
dup_tbl <- bind_rows(dup1, dup2) %>%
  select(library, everything())

# ----------------------------
# Shared/private exact sequences
# ----------------------------

shared_exact <- tibble(
  metric = c("shared_exact_sequences",
             "mine_private_exact_sequences",
             "their_private_exact_sequences"),
  value = c(
    length(intersect(unique(d1$sequence), unique(d2$sequence))),
    length(setdiff(unique(d1$sequence), unique(d2$sequence))),
    length(setdiff(unique(d2$sequence), unique(d1$sequence)))
  )
)

# ----------------------------
# Header/locus/family summaries
# ----------------------------

family_tbl <- d %>%
  group_by(library, family_like) %>%
  summarise(
    n_records = n(),
    total_bp = sum(length),
    .groups = "drop"
  ) %>%
  arrange(library, desc(total_bp), desc(n_records))

locus_tbl <- d %>%
  group_by(library) %>%
  summarise(
    n_nonmissing_locus_like = sum(!is.na(locus_like)),
    n_unique_locus_like = n_distinct(locus_like, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------
# K-mer overlap
# ----------------------------
# k=21 is a reasonable starting point for RNA-seq-like uniqueness checks.
# If the FASTAs are huge, this may take a while and use memory.

message("Building k-mer sets...")
km1 <- kmer_set(d1$sequence, k = 21)
km2 <- kmer_set(d2$sequence, k = 21)

kmer_tbl <- tibble(
  metric = c(
    "mine_unique_21mers",
    "their_unique_21mers",
    "shared_21mers",
    "jaccard_21mers",
    "fraction_of_mine_21mers_private",
    "fraction_of_their_21mers_private"
  ),
  value = c(
    length(km1),
    length(km2),
    length(intersect(km1, km2)),
    jaccard(km1, km2),
    private_kmer_fraction(km1, km2),
    private_kmer_fraction(km2, km1)
  )
)

# ----------------------------
# Write tables
# ----------------------------

write.csv(summary_tbl, paste0(prefix, "_summary.csv"), row.names = FALSE)
write.csv(dup_tbl, paste0(prefix, "_duplicate_stats.csv"), row.names = FALSE)
write.csv(shared_exact, paste0(prefix, "_shared_exact_sequences.csv"), row.names = FALSE)
write.csv(family_tbl, paste0(prefix, "_family_summary.csv"), row.names = FALSE)
write.csv(locus_tbl, paste0(prefix, "_locus_summary.csv"), row.names = FALSE)
write.csv(kmer_tbl, paste0(prefix, "_kmer_summary.csv"), row.names = FALSE)

# per-sequence table if you want to inspect specific outliers
write.csv(d, paste0(prefix, "_per_sequence_metrics.csv"), row.names = FALSE)

# ----------------------------
# Plots
# ----------------------------

p_len <- ggplot(d, aes(x = length, fill = library)) +
  geom_histogram(bins = 80, alpha = 0.5, position = "identity") +
  scale_x_log10() +
  theme_bw() +
  ggtitle("Sequence length distribution (log10 x-axis)")

ggsave(paste0(prefix, "_length_distribution.png"), p_len, width = 8, height = 5, dpi = 300)

p_gc <- ggplot(d, aes(x = gc_percent, fill = library)) +
  geom_histogram(bins = 60, alpha = 0.5, position = "identity") +
  theme_bw() +
  ggtitle("GC% distribution")

ggsave(paste0(prefix, "_gc_distribution.png"), p_gc, width = 8, height = 5, dpi = 300)

p_complex <- ggplot(d, aes(x = complexity, fill = library)) +
  geom_histogram(bins = 60, alpha = 0.5, position = "identity") +
  theme_bw() +
  ggtitle("Shannon sequence complexity distribution")

ggsave(paste0(prefix, "_complexity_distribution.png"), p_complex, width = 8, height = 5, dpi = 300)

# top families by bp
top_fams <- family_tbl %>%
  group_by(library) %>%
  slice_max(order_by = total_bp, n = 20, with_ties = FALSE) %>%
  ungroup()

p_family <- ggplot(top_fams, aes(x = reorder(family_like, total_bp), y = total_bp, fill = library)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw() +
  ggtitle("Top family-like header groups by total bp") +
  xlab("family_like")

ggsave(paste0(prefix, "_top_family_bp.png"), p_family, width = 9, height = 7, dpi = 300)

message("Done.")
message("Wrote:")
message("  ", paste0(prefix, "_summary.csv"))
message("  ", paste0(prefix, "_duplicate_stats.csv"))
message("  ", paste0(prefix, "_shared_exact_sequences.csv"))
message("  ", paste0(prefix, "_family_summary.csv"))
message("  ", paste0(prefix, "_locus_summary.csv"))
message("  ", paste0(prefix, "_kmer_summary.csv"))
message("  ", paste0(prefix, "_per_sequence_metrics.csv"))
message("  ", paste0(prefix, "_length_distribution.png"))
message("  ", paste0(prefix, "_gc_distribution.png"))
message("  ", paste0(prefix, "_complexity_distribution.png"))
message("  ", paste0(prefix, "_top_family_bp.png"))