## ----install packages---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# remotes::install_github("mojaveazure/seurat-disk")
# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }

# BiocManager::install("GRaNIE")
# BiocManager::install(c("TxDb.Hsapiens.UCSC.hg38.knownGene", "BSgenome.Hsapiens.UCSC.hg38"))
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("ChIPseeker")


## ----load packages------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(readr)
library(GRaNIE)
library(dplyr)
library(org.Hs.eg.db)


## ----load data----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
datafolder <- "/projectnb/paxlab/EnhancerDiscovery/data"
genes_file <- paste0(
  datafolder,
  "/Bulk-Ribo-Depleted-RNAseq-Heart-Lung-Liver_Human_FBs/hFB_riboDepletedBulkRnaSeq/NF_out/star_salmon/salmon.merged.transcript_tpm.tsv"
)
peaks_file <- paste0(
  datafolder,
  "/Bulk-ATAC-Heart-Lung-Liver_Human_FBs/NF_out/bwa/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.featureCounts.txt"
)

counts_rna_df <- read_tsv(genes_file, col_types = cols())
counts_peaks_df <- read_tsv(peaks_file, col_types = cols(), skip = 1)


## ----average data-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
organ_names <- c("Heart", "Lung", "Liver")
condition_names <- c("Unstim", "TGFb", "TGFb+JQ1")
rna_cols <- c("Unstim", "TGFb", "TGFb_._JQ1")
atac_cols <- c("Untim", "TGFB", "TGFB_JQ1")

counts_rna_avg_df <- data.frame(ENSEMBL = counts_rna_df$gene_id)
for (organ in organ_names) {
  for (condition in rna_cols) {
    counts_rna_avg_df[paste(organ, condition_names[match(condition, rna_cols)], sep = "_")] <-
      rowMeans(counts_rna_df %>% dplyr::select(matches(paste0(organ, ".+", condition, "_...$"))))
  }
}

counts_peaks_avg_df <- data.frame(peakID = paste0(counts_peaks_df$Chr, ":", counts_peaks_df$Start, "-", counts_peaks_df$End))
for (organ in organ_names) {
  for (condition in atac_cols) {
    counts_peaks_avg_df[paste(organ, condition_names[match(condition, atac_cols)], sep = "_")] <-
      rowMeans(counts_peaks_df %>% dplyr::select(matches(paste(condition, organ, sep = "_"))))
  }
}

metadata_df <- data.frame(sample_id = colnames(counts_peaks_avg_df)[-1])
metadata_df$organ <- unlist(strsplit(metadata_df$sample_id, "_"))[c(TRUE, FALSE)]
metadata_df$condition <- unlist(strsplit(metadata_df$sample_id, "_"))[c(FALSE, TRUE)]


## ----convert gene ids---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
counts_rna_avg_df$ENSEMBL <- mapIds(org.Hs.eg.db, keys = counts_rna_avg_df$ENSEMBL, column = "ENSEMBL", keytype = "SYMBOL")


## ----initialize GRaNIE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
outputFolder <- "./data"
grn <- initializeGRN(
  outputFolder = outputFolder,
  genomeAssembly = "hg38"
)
grn_file_rds <- paste0(outputFolder, "/GRN.rds")
tfbs_folder <- paste0(outputFolder, "/H12INVIVO")
grn


## ----add data to GRaNIE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
grn <- addData(grn,
  counts_peaks = counts_peaks_avg_df, normalization_peaks = "none",
  counts_rna = counts_rna_avg_df, normalization_rna = "none",
  sampleMetadata = metadata_df, forceRerun = TRUE
)
grn


## ----qc time!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
grn <- plotPCA_all(grn)


## ----add TFBS-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
grn <- addTFBS(grn,
  motifFolder = tfbs_folder, TFs = "all",
  filesTFBSPattern = "_TFBS", fileEnding = ".bed.gz", forceRerun = TRUE
)
nCores <- as.numeric(Sys.getenv("NSLOTS"))
grn <- overlapPeaksAndTFBS(grn, nCores = 1, forceRerun = TRUE)


## ----tf-enhancer connections--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
grn <- addConnections_TF_peak(grn,
  plotDiagnosticPlots = FALSE, connectionTypes = c("expression"),
  corMethod = "pearson", forceRerun = TRUE
)


## ----gene-enhancer connections------------------------------------------------------------------------------------------------------------------------------------------------------------------------
grn <- addConnections_peak_gene(grn)
grn <- plotDiagnosticPlots_peakGene(grn,
  gene.types = list(c("protein_coding", "lincRNA")),
  plotAsPDF = FALSE, pages = 1
)


## ----save grn to disk---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(grn, grn_file_rds)
