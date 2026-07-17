# run with 32 cores

PROJECT_DIR <- "./ST8SIA4_GWAS_molQTL/"

CHR <- 5
START_ST8SIA4_b38 <- 100806933
END_ST8SIA4_b38 <- 100903282

START_ST8SIA4_b37 <- 100142637
END_ST8SIA4_b37 <- 100238986
geneID <- "ENSG00000113532"

library(data.table)
library(dplyr)
library(ggplot2)
library(parallel)


###################################################################################################
######## ST8SIA4 coverage tests ###################################################################
###################################################################################################

EGA_samples <- fread(sprintf("%s/data/02_fastq/EGA_RNA_seq_annotation.txt", PROJECT_DIR))
# EGA_samples_corresp <- fread(sprintf("%s/EvoImmunoPop_EGA_corr.txt", PROJECT_DIR)) #(used internally for correspondece with files using older/defunct sample IDs)
# EGA_samples <- merge(EGA_samples, EGA_samples_corresp[, .(ID = subject_ID, sample_ID = old)], by = "ID")


# coverage_bin <- list()
# junction_list <- list()

process_sample_coverage <- function(sample_id, DIR = PROJECT_DIR) {
  cat("coverage", sample_id, "\n")
  library(data.table)
  CHR <- 5
  START_ST8SIA4_b38 <- 100806933
  END_ST8SIA4_b38 <- 100903282
  geneID <- "ENSG00000113532"
  sample_file_coverage <- sprintf("%s/04_coverage/sashimi_ST8SIA4_%s__density_list_%s_aligned.sortedByCoord.out.ST8SIA4_b38.txt", DIR, sample_id, sample_id)
  coverage <- try(fread(sample_file_coverage))
  if (any(class(coverage) == "try-error")) {
    coverage_bin <- NULL
  } else {
    coverage <- coverage[x > START_ST8SIA4_b38 - 1e4 & x < END_ST8SIA4_b38 + 1e4, ]
    coverage[, bin := cut(x, breaks = unique(c(seq(99806, 101904, by = 2) * 1000)))]
    coverage_bin <- coverage[, .(x = mean(x), y = sum(y), mean_coverage = mean(y), x_min = min(x), x_max = max(x)), by = .(bin)]
  }
  coverage_bin
}


process_sample_junction <- function(sample_id, DIR = PROJECT_DIR) {
  library(data.table)
  cat("coverage", sample_id, "\n")
  sample_file_junction <- sprintf("%s/04_coverage/sashimi_ST8SIA4_%s__junction_list_%s_aligned.sortedByCoord.out.ST8SIA4_b38.txt", DIR, sample_id, sample_id)
  junc <- try(fread(sample_file_junction))
  if (any(class(junc) == "try-error")) {
    junc <- NULL
  }
  junc
}

# Create a cluster with 32 cores
num_cores <- 32
cl <- makeCluster(num_cores)
# Parallelize the loop with parLapply()
coverage_bins <- parLapply(cl, EGA_samples[, ID], process_sample_coverage)
names(coverage_bins) <- EGA_samples[, ID]
coverage_bins <- rbindlist(coverage_bins, idcol = "sample")
coverage_bins[, relative_coverage := mean_coverage / sum(mean_coverage), by = sample]

fwrite(coverage_bins, file = sprintf("%s/04_coverage/aggregate/T8SIA4_coverage_bins_all_samples.tsv.gz", PROJECT_DIR), sep = "\t")
# Parallelize the loop with parLapply()
junctions <- parLapply(cl, EGA_samples[, ID], process_sample_junction)
names(junctions) <- EGA_samples[, ID]
junctions <- rbindlist(junctions, idcol = "sample")

fwrite(junctions, file = sprintf("%s/04_coverage/aggregate/ST8SIA4_junctions_all_samples.tsv.gz", PROJECT_DIR), sep = "\t")
# Stop the cluster after completion
stopCluster(cl)



