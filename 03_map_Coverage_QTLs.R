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


junctions <- fread(sprintf("%s/data/04_coverage/ST8SIA4_junctions_all_samples.tsv.gz", PROJECT_DIR))
coverage_bins <- fread(sprintf("%s/data/04_coverage/ST8SIA4_coverage_bins_all_samples.tsv.gz", PROJECT_DIR))
coverage_bins <- merge(coverage_bins, EGA_samples, by.x = "sample", by.y = "ID")

###################################################################################################
###  extract genotypes (Quach et al 2016, b37) ####################################################
###################################################################################################

library(VariantAnnotation)

# input (obatianed from )
Genotype_path_chr5_200ind_QuachEtAl <- sprintf("%s/data/05_genotype/EvoImmunoPop_imputation_200x19619457_chr5", PROJECT_DIR)

# output with subset of SNPs near ST8SIA4
Genotype_path_ST8SIA4_200ind_QuachEtAl <- sprintf("%s/data/05_genotype/EvoImmunoPop_imputation_200x19619457_ST8SIA4", PROJECT_DIR)


plink_command <- paste(
  "plink",
  "--bfile", genotype_path_chr5_200ind_QuachEtAl,
  "--chr", CHR,
  "--from-bp", START_ST8SIA4_b37 - 1e5,
  "--to-bp", END_ST8SIA4_b37 + 1e5,
  "--make-bed",
  "--out", Genotype_path_ST8SIA4_200ind_QuachEtAl
)

b37_genotypes <- system(plink_command)


# Load the snpStats package
library(snpStats)

# Read the data using read.plink
plink_data <- read.plink(
  bed = paste0(Genotype_path_ST8SIA4_200ind_QuachEtAl, ".bed"),
  bim = paste0(Genotype_path_ST8SIA4_200ind_QuachEtAl, ".bim"),
  fam = paste0(Genotype_path_ST8SIA4_200ind_QuachEtAl, ".fam")
)

# Explore the resulting data
names(plink_data)

# Access the genotype data
genotypes <- 2 - as(plink_data$genotypes, "numeric")
str(genotypes)

# Access subject data
subject_data <- as.data.table(plink_data$fam)
print(subject_data)

# Access SNP map data
snp_data <- as.data.table(plink_data$map)
print(snp_data)
setnames(snp_data, c("allele.1", "allele.2"), c("effect_allele", "ref_allele"))

bin_annot <- unique(coverage_bins[, .(bin, chr = 5, start = x_min, end = x_max)])
snp_annot <- snp_data[, .(rsID = snp.name, chromosome, position)]
library(MatrixEQTL)

coverage_bins[, IID := substr(sample_ID, 1, 6)]


############################################################################
###  map coverage QTLs  ####################################################
############################################################################


QTL_pvalues <- list()
covQTL_sumStats <- list()
for (myCOND in c("NS", "LPS", "PAM3CSK4", "IAV", "R848")) {
  for (PERM in 0:100) {
    relative_coverage_matrix <- dcast(coverage_bins[condition == myCOND, ], bin ~ IID, value.var = "relative_coverage", fill = 0)
    covariate_matrix <- unique(coverage_bins[condition == myCOND, .(IID, GC_pct, Five_Prime_Bias, Pct_European_Ancestry, age)])
    covariate_mat <- as.matrix(covariate_matrix[, -1])
    rownames(covariate_mat) <- covariate_matrix$IID

    # ranktransform gene expression data
    rankTransform <- function(x) {
      require(DescTools)
      notNA <- which(!is.na(x))
      percentile <- rank(x[notNA], ties.method = "random", na.last = NA) / (length(x) + 1)
      mean_level <- mean(DescTools::Winsorize(x[notNA]))
      sd_level <- sd(DescTools::Winsorize(x[notNA]))
      x[notNA] <- qnorm(percentile, mean_level, sd_level)
      x
    }

    set.seed(123)

    # format for MatrixeQTL

    feature_matrix <- t(relative_coverage_matrix[, -1])
    colnames(feature_matrix) <- relative_coverage_matrix$bin
    feature_mat <- t(apply(feature_matrix, 2, rankTransform))

    Feature_Mat <- SlicedData$new()
    Feature_Mat$CreateFromMatrix(feature_mat)
    Feature_Mat$ResliceCombined(sliceSize = 5000)
    # show(Feature_Mat)

    samp <- match(rownames(feature_matrix), rownames(genotypes))
    if (PERM > 0) {
      set.seed(PERM)
      samp <- sample(samp)
    }

    Predictor_Mat <- SlicedData$new()
    Predictor_Mat$CreateFromMatrix(t(genotypes[samp, ]))
    Predictor_Mat$ResliceCombined(sliceSize = 5000)
    # show(Predictor_Mat)

    Cov_Mat <- SlicedData$new()
    Cov_Mat$CreateFromMatrix(t(covariate_mat[match(rownames(feature_matrix), rownames(covariate_mat)), ]))
    # show(Cov_Mat)

    # run MatrixeQTL

    res <- Matrix_eQTL_main(Predictor_Mat,
      Feature_Mat,
      Cov_Mat,
      output_file_name = NULL,
      pvOutputThreshold = 0,
      output_file_name.cis = NULL,
      pvOutputThreshold.cis = 1,
      cisDist = 1e6,
      snpspos = snp_annot,
      genepos = bin_annot,
      min.pv.by.genesnp = TRUE,
      verbose = FALSE
    )

    Pvalues <- data.table(
      featureID = names(res$cis$min.pv.gene),
      pvalue = res$cis$min.pv.gene,
      celltype = "MONO.CD14",
      condition = myCOND,
      population = "ALL",
      perm = PERM
    )
    QTL_pvalues[[paste(PERM, myCOND)]] <- Pvalues[, mget(c("perm", "celltype", "condition", "population", "featureID", "pvalue"))]
    cat(myCOND, PERM, "")
    if (PERM == 0) {
      covQTL_sumStats[[myCOND]] <- as.data.table(res$cis$eqtls)
    }
    cat("\n\n", myCOND, ":")
  }
}
QTL_pvalues <- rbindlist(QTL_pvalues)
covQTL_sumStats <- rbindlist(covQTL_sumStats, idcol = "condition")
fwrite(covQTL_sumStats, file = sprintf("%s/data/06_sum_stats_covQTL/coverage_QTL_ST8SIA4_sumStats.txt.gz", PROJECT_DIR), sep = "\t")
fwrite(QTL_pvalues, file = sprintf("%s/data/06_sum_stats_covQTL/coverage_QTL_ST8SIA4_bestP_perFeature.txt.gz", PROJECT_DIR), sep = "\t")

feature_signif <- QTL_pvalues[, .(P_perm = mean(pvalue[perm > 0] < pvalue[perm == 0]), E0 = mean(pvalue[perm > 0]), V0 = var(pvalue[perm > 0]), obs = pvalue[perm == 0]), by = .(condition, featureID)]
feature_signif[, a_0 := E0 * (E0 * (1 - E0) / V0 - 1)]
feature_signif[, b_0 := (1 - E0) * (E0 * (1 - E0) / V0 - 1)]
feature_signif[, P_param := pbeta(obs, a_0, b_0, lower = T)]
feature_signif[, FDR := p.adjust(P_param, "fdr"), by = .(condition)]
fwrite(feature_signif, file = sprintf("%s/data/06_sum_stats_covQTL/coverage_QTL_ST8SIA4_perFeature_FDR.tsv", PROJECT_DIR), sep = "\t")

# reload output for QC
# QTL_pvalues <- fread(sprintf("%s/sum_stats_covQTL/coverage_QTL_ST8SIA4_bestP_perFeature.txt.gz", PROJECT_DIR))
# covQTL_sumStats <- fread(sprintf("%s/sum_stats_covQTL/coverage_QTL_ST8SIA4_sumStats.txt.gz", PROJECT_DIR))

