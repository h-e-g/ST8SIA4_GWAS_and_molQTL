options(stringsAsFactors = FALSE, max.print = 9999, width = 300, datatable.fread.input.cmd.message = FALSE)

.libPaths(sprintf("%s/single_cell/resources/R_libs/4.1.0", EIP))

library(data.table)
library(tictoc)
library(coloc)
library(susieR)
library(snpStats)

PROJECT_DIR <- "./ST8SIA4_GWAS_molQTL/"

CHR <- 5
START_ST8SIA4_b38 <- 100806933
END_ST8SIA4_b38 <- 100903282
CIS_DIST <- 1e6

START_ST8SIA4_b37 <- 100142639
END_ST8SIA4_b37 <- 100238970


#################### load GWAS sumStats ##########
sumStats_GWAS <- fread(sprintf("%s/data/01_GWAS_results/03_GWAS_summary_statistics/01_sumStats_selected_traits_ST8SIA4_locus.tsv", PROJECT_DIR))

#################### load genotypes ###############
Genotype_path_chr5_200ind_QuachEtAl <- sprintf("%s/data/05_genotype/EvoImmunoPop_imputation_200x19619457_chr5", PROJECT_DIR)
Genotype_path_ST8SIA4_200ind_QuachEtAl <- sprintf("%s/data/05_genotype/EvoImmunoPop_imputation_200x19619457_ST8SIA4", PROJECT_DIR)

# Read the data using read.plink
plink_data <- read.plink(
  bed = paste0(Genotype_path_ST8SIA4_200ind_QuachEtAl, ".bed"),
  bim = paste0(Genotype_path_ST8SIA4_200ind_QuachEtAl, ".bim"),
  fam = paste0(Genotype_path_ST8SIA4_200ind_QuachEtAl, ".fam")
)
# Access the genotype data
genotypes <- 2 - as(plink_data$genotypes, "numeric")
genotypes_char <- 2 - as(plink_data$genotypes, "character")
str(genotypes)

# Access SNP map data
snp_data <- as.data.table(plink_data$map)
print(snp_data)
setnames(snp_data, c("allele.1", "allele.2"), c("effect_allele", "ref_allele"))


########### read covQTL sum stats ###############
covQTL_sumStats <- fread(sprintf("%s/data/06_sum_stats_covQTL/coverage_QTL_ST8SIA4_sumStats.txt.gz", PROJECT_DIR))
feature_signif <- fread(sprintf("%s/data/06_sum_stats_covQTL/coverage_QTL_ST8SIA4_perFeature_FDR.tsv", PROJECT_DIR))
covQTL_sumStats_signif <- covQTL_sumStats[paste(gene, condition) %in% feature_signif[FDR < .01, paste(featureID, condition)], ]
covQTL_sumStats_signif <- merge(covQTL_sumStats_signif, snp_data[, .(snps = snp.name, POS_b37 = position, effect_allele, ref_allele)], by = "snps", all.x = TRUE)

## add GWAS Credible sets info to covQTL sumStats
CS_selected_GWAS_traits <- fread(sprintf("%s/data/01_GWAS_results/05_CS_selected_GWAS_traits.tsv", PROJECT_DIR))
covQTL_sumStats_signif <- merge(covQTL_sumStats_signif, CS_selected_GWAS_traits[, .(hm_rsid, is_CS_ID, hm_variant_id)], by.x = "snps", by.y = "hm_rsid", all.x = TRUE)


############################################################
tic("computing LD matrix") #################################
############################################################

Geno_freq <- apply(genotypes, 2, mean, na.rm = TRUE)
common <- Geno_freq > 0.05 & Geno_freq < 0.95
Geno_res <- apply(genotypes[, which(common)], 2, function(x) {
  y <- rep(NA, length(x))
  y[!is.na(x)] <- lm(x[!is.na(x)] ~ substr(rownames(genotypes)[!is.na(x)], 1, 3))$res
  y
})
rownames(Geno_res) <- rownames(genotypes)
LDMAP_full <- cor(Geno_res, use = "p")
toc()

## select covQTL for common SNPs only
covQTL_sumStats_signif[, Freq_effect := Geno_freq[match(snps, names(Geno_freq))]]
covQTL_sumStats_signif <- covQTL_sumStats_signif[Freq_effect > .05 & Freq_effect < .95, ]


CIS_DIST <- 1e5


## merge SummStats from GWAS and covQTL 
SumStats_both <- merge(covQTL_sumStats_signif[, .(variant_id = snps, POS_b37, p_value = pvalue, beta.eQTL = beta, se.eQTL = beta / statistic, bin = gene, condition, effect_allele, ref_allele)],
  sumStats_GWAS[, .(variant_id, POS_b38 = base_pair_location, effect_allele, other_allele, beta.pheno = beta, se.pheno = standard_error, p_value, studyId)],
  by = "variant_id",
  suffix = c(".eQTL", ".pheno"),
  allow.cartesian = TRUE
)
SumStats_both[, effect_allele.eQTL := toupper(effect_allele.eQTL)]
SumStats_both[, effect_allele.pheno := toupper(effect_allele.pheno)]
SumStats_both[, other_allele.eQTL := toupper(ref_allele)]
SumStats_both[, other_allele.pheno := toupper(other_allele)]
SumStats_both[, ref_allele := NULL]
SumStats_both[, other_allele := NULL]

SumStats_both[effect_allele.eQTL != effect_allele.pheno, beta.eQTL := -beta.eQTL]
SumStats_both[effect_allele.eQTL != effect_allele.pheno, beta.eQTL := -beta.eQTL]
SumStats_both[effect_allele.eQTL != effect_allele.pheno, other_allele.eQTL := other_allele.pheno]
SumStats_both[effect_allele.eQTL != effect_allele.pheno, effect_allele.eQTL := effect_allele.pheno]

SumStats_both <- SumStats_both[!is.na(beta.eQTL) & !is.na(beta.pheno), ]
SumStats_both[, z.eQTL := beta.eQTL / se.eQTL]
SumStats_both[, z.pheno := beta.pheno / se.pheno]
SumStats_both <- SumStats_both[!duplicated(paste(variant_id, bin, condition, studyId)), ]
SumStats_both[, N_snps := .N, by = .(bin, condition, studyId)]
SumStats_both <- SumStats_both[N_snps > 50, ]



################################################################################
################ colocalization and fine mapping functions  ####################
################################################################################

################# get coloclaization results 
get_coloc <- function(sumStat, LDMAP_full) {
  LDMAP <- LDMAP_full[sumStat[, variant_id], sumStat[, variant_id]]

  data_eQTL <- list(
    beta = setNames(sumStat[, beta.eQTL], sumStat[, variant_id]),
    varbeta = setNames(sumStat[, se.eQTL^2], sumStat[, variant_id]),
    snp = sumStat[, variant_id],
    position = sumStat[, POS_b38],
    LD = LDMAP,
    type = "quant",
    sdY = 1
  )

  data_pheno <- list(
    beta = setNames(sumStat[, beta.pheno], sumStat[, variant_id]),
    varbeta = setNames(sumStat[, se.pheno^2], sumStat[, variant_id]),
    snp = sumStat[, variant_id],
    position = sumStat[, POS_b38],
    LD = LDMAP,
    type = "quant",
    sdY = 1
  )

  tic("running coloc")
  coloc.res <- try(coloc.signals(data_eQTL, data_pheno, p12 = 1e-5)$summary)
  toc()
  coloc.res
}

################# get credible sets
get_CS <- function(sumStat, LDMAP_full, SUSIE_COVERAGE = 0.95, N_PRED_MAX = 10, N_IND_POP = 200, zval = "z.eQTL", variant = "variant_id") {
  LDMAP <- LDMAP_full[sumStat[, get(variant)], sumStat[, get(variant)]]
  susie.obj <- susie_rss(sumStat[, get(zval)], LDMAP, coverage = SUSIE_COVERAGE, n = N_IND_POP, L = N_PRED_MAX, niter = 200)
  CSlist <- susie.obj$sets$cs
  credible_set <- rep("", nrow(sumStat))
  for (i in seq_along(CSlist)) {
    credible_set[sumStat[, get(variant)] %in% sumStat[CSlist[[i]], get(variant)]] <- names(CSlist)[i]
  }
  credible_set
}


################# get psoterior probab of inclusion 
get_PIP <- function(sumStat, LDMAP_full, SUSIE_COVERAGE = 0.95, N_PRED_MAX = 10, N_IND_POP = 200, zval = "z.eQTL", variant = "variant_id") {
  LDMAP <- LDMAP_full[sumStat[, get(variant)], sumStat[, get(variant)]]
  susie.obj <- susie_rss(sumStat[, get(zval)], LDMAP, coverage = SUSIE_COVERAGE, n = N_IND_POP, L = N_PRED_MAX, niter = 200)

  top_component <- apply(susie.obj$alpha, 2, which.max)
  pip_top_component <- susie.obj$alpha[cbind(top_component, seq_along(top_component))]
  pip_top_component
}

################################################################################
################################################################################


#############################################################################
############### add finamapiing stats to covQTL summary stats  ##############
#############################################################################

SumStats_both[studyId == "GCST90002316" & bin == "(1.00904e+08,1.00906e+08]" & condition == "NS", get_coloc(.SD, LDMAP_full)]
coloc_results <- SumStats_both[, credible_set := get_CS(.SD, LDMAP_full), by = .(studyId, bin, condition)]
fwrite(coloc_results, file = sprintf("%s/data/06_sum_stats_covQTL/coloc_results_covQTL.tsv", PROJECT_DIR), sep = "\t")


SumStats_both[, credible_set.eQTL := get_CS(.SD, LDMAP_full), by = .(studyId, bin, condition)]
SumStats_both[, pip.eQTL := get_PIP(.SD, LDMAP_full), by = .(studyId, bin, condition)]
fwrite(SumStats_both, file = sprintf("%s/data/06_sum_stats_covQTL/finemapped_covQTL_both.tsv", PROJECT_DIR), sep = "\t")

covQTL_sumStats_signif[, credible_set := get_CS(.SD, LDMAP_full, zval = "statistic", variant = "snps"), by = .(gene, condition)]
covQTL_sumStats_signif[, PIP := get_PIP(.SD, LDMAP_full, zval = "statistic", variant = "snps"), by = .(gene, condition)]

fwrite(covQTL_sumStats_signif, file = sprintf("%s/data/06_sum_stats_covQTL/finemapped_covQTL_signif.tsv", PROJECT_DIR), sep = "\t")

#########################################################################################
############### aggregate covQTL to have same index SNP for overlapping CS ##############
#########################################################################################

Assoc_CS <- covQTL_sumStats_signif[credible_set != "", ]
bestPredictor_list <- NULL
Assoc_CS_remaining <- Assoc_CS
# Assoc_CS: 1 line per feature, celltype , condition  & snp in the CS of a significant component CS for that celltype & condition
# explained_components: 1 line per feature, component, celltype & condition (where component has a significant effect)à
# bestPredictor_list: 1 line per feature & independent eQTL
cur_round <- 0
while (nrow(Assoc_CS_remaining) > 0) {
  cur_round <- cur_round + 1
  cat(nrow(bestPredictor_list), "independent covQTL identified - ", nrow(Assoc_CS_remaining), "predictor-feature association remaining\n")
  predictor_score <- Assoc_CS_remaining[, .(pred_score = sum(abs(PIP))), by = .(snps)]
  bestPredictor <- predictor_score[order(-pred_score), head(.SD, 1)]
  bestPredictor[, round := cur_round]
  explained_components <- merge(Assoc_CS_remaining, bestPredictor, by = c("snps"))[, .(gene, snps, condition, credible_set)]
  predictor_score <- Assoc_CS_remaining[, .(pred_score = sum(abs(PIP))), by = .(snps)]
  Assoc_CS_remaining <- Assoc_CS_remaining[!paste(gene, condition, credible_set, sep = "_") %chin% explained_components[, paste(gene, condition, credible_set, sep = "_")]]
  bestPredictor_list <- rbind(bestPredictor_list, bestPredictor)
}
explained_components <- merge(Assoc_CS, bestPredictor_list, by = c("gene", "snps"))
fwrite(explained_components[order(round, gene, pvalue), ], file = sprintf("%s/data/06_sum_stats_covQTL/finemapped_covQTL_explained.tsv", PROJECT_DIR), sep = "\t")
#explained_components <- fread(sprintf("%s/sum_stats_covQTL/finemapped_covQTL_explained.tsv", PROJECT_DIR))

########################################################################
############### generate Table of colocalization results  ##############
########################################################################

OPENTARGETS_DIR <- sprintf("%s/data/01_GWAS_results/04_openTargets_25_03/", PROJECT_DIR)
list_COLOC_parquet_files <- dir(sprintf("%s/coloc", OPENTARGETS_DIR), pattern = "parquet")

coloc_results <- list()
for (FILE in list_COLOC_parquet_files) {
  coloc_results[[FILE]] <- as.data.table(read_parquet(sprintf("%s/coloc/%s", OPENTARGETS_DIR, FILE)))
}
coloc_results <- rbindlist(coloc_results)

coloc_ST8SIA4 <- coloc_results[leftStudyLocusId %chin% ST8SIA4_assoc_list[,studyLocusId] & rightStudyLocusId %chin% ST8SIA4_assoc_list[,studyLocusId]]
coloc_ST8SIA4 <- merge(coloc_ST8SIA4, ST8SIA4_assoc_list, by.x = "leftStudyLocusId", by.y = "studyLocusId")
coloc_ST8SIA4 <- merge(coloc_ST8SIA4, ST8SIA4_assoc_list, by.x = "rightStudyLocusId", by.y = "studyLocusId",suffix=c(".right",".left"))
coloc_matrix <- as.matrix(dcast(coloc_ST8SIA4, studyId.left ~ studyId.right, value.var = "h4")[, -1])
rownames(coloc_matrix) <- dcast(coloc_ST8SIA4, studyId.left ~ studyId.right, value.var = "h4")[, studyId.left]


##### molQTL from eQTL catalog 
coloc_results_tableS1a=coloc_ST8SIA4[studyId.right%in%c('GCST011096','GCST90002340','GCST90002316'),][order(h4)][,.(gwasID=studyId.right,best_GWAS_variantID=variantId.right,molQTL_phenoID=studyId.left,best_molQTL_variantID=variantId.left, numberColocalisingVariants,h0,h1,h2,h3,h4,colocalisationMethod,betaRatioSignAverage)]
coloc_results_tableS1a=merge(coloc_results_tableS1a,
                            ST8SIA4_finemap_better[,.(molQTL_phenoID=studyId,best_molQTL_variantID=variantId,posteriorProbability,pval=pValueMantissa*10^pValueExponent,beta)],by=c('molQTL_phenoID','best_molQTL_variantID'),all.x=TRUE)
coloc_results_tableS1a=merge(coloc_results_tableS1a,
                            sumStats[,.(gwasID=studyId,best_GWAS_variantID=hm_variant_id,pval=p_value,beta=hm_beta)],by=c('gwasID','best_GWAS_variantID'),suffix=c('.molQTL','.gwas'),all.x=TRUE)

coloc_results_tableS1a <- merge(coloc_results_tableS1a, unique(sumStats[, .(best_molQTL_variantID = hm_variant_id, best_molQTL_rsID = hm_rsid)]), by = c("best_molQTL_variantID"), all.x = TRUE)
coloc_results_tableS1a <- merge(coloc_results_tableS1a, unique(sumStats[, .(best_GWAS_variantID = hm_variant_id, best_GWAS_rsID = hm_rsid)]), by = c("best_GWAS_variantID"), all.x = TRUE)
coloc_results_tableS1a <- coloc_results_tableS1a[order(gwasID, molQTL_phenoID, -h4), .(gwasID, molQTL_phenoID, condition, best_GWAS_variantID, best_GWAS_rsID, best_molQTL_variantID, best_molQTL_rsID, numberColocalisingVariants, h0, h1, h2, h3, h4, colocalisationMethod, betaRatioSignAverage, pval.gwas, beta.gwas, pval.molQTL, beta.molQTL)]
coloc_results_tableS1a[best_molQTL_variantID == "5_100809449_C_CA", best_molQTL_rsID := "rs35431964"] # two bp shortening of polyA 100850954-100850977 (+)
coloc_results_tableS1a[best_molQTL_variantID == "5_100850954_CTT_C", best_molQTL_rsID := "rs70987828"] # two bp shortening of polyA 100850954-100850977 (-)
coloc_results_tableS1a[best_molQTL_variantID == "5_100853985_A_ATGTGTGTGTGTG", best_molQTL_rsID := "rs55821438"]
coloc_results_tableS1a[best_molQTL_variantID == "5_100872875_TAA_T", best_molQTL_rsID := "rs550303574"] # two bp shortening of a polyA 100872875-100872895 (+)
coloc_results_tableS1a[best_molQTL_variantID == "5_100864593_A_G", best_molQTL_rsID := "rs200660720"] # disrupts a polyA 100864574-100864603 (+)
fwrite(coloc_results_tableS1a, file = sprintf("%s/06_sum_stats_covQTL/TableS1a_coloc_results_moQTL_eQTLcatalog.tsv", PROJECT_DIR), sep = "\t")

##### covQTL from our study

coloc_results_tableS1c <- merge(coloc_results, bin_annot[, .(bin, bin_label2)], by = "bin")
beta_ratio <- SumStats_both[p_value.eQTL < .01 & p_value.pheno < 0.01, .(betaRatioSignAverage = sign(mean(beta.eQTL / beta.pheno * ifelse(other_allele.eQTL == other_allele.pheno, 1, -1)))), by = .(bin, condition, studyId)]
coloc_results_tableS1c <- merge(coloc_results_tableS1c, beta_ratio, by = c("bin", "condition", "studyId"))
coloc_results_tableS1c <- coloc_results_tableS1c[studyId %in% c("GCST011096", "GCST90002340", "GCST90002316"), ][order(studyId, bin)][, .(bin, gwasID = studyId, best_GWAS_rsID = best2, molQTL_phenoID = bin_label2, condition, best_molQTL_rsID = best1, numberColocalisingVariants = nsnps, h0 = PP.H0.abf, h1 = PP.H1.abf, h2 = PP.H2.abf, h3 = PP.H3.abf, h4 = PP.H4.abf, colocalisationMethod = "COLOC", betaRatioSignAverage)]
coloc_results_tableS1c <- merge(coloc_results_tableS1c, SumStats_both[, .(bin, condition, gwasID = studyId, best_GWAS_rsID = variant_id, pval.gwas = p_value.pheno, beta.gwas = beta.pheno, pval.molQTL = p_value.eQTL, beta.molQTL = beta.eQTL)], by = c("gwasID", "best_GWAS_rsID", "bin", "condition"), suffix = c(".molQTL", ".gwas"), all.x = TRUE)
# coloc_results_tableS1c=merge(coloc_results_tableS1c,sumStats[,.(gwasID=studyId,best_GWAS_variantID=hm_rsid,pval=p_value,beta=hm_beta)],by=c('gwasID','best_GWAS_variantID'),suffix=c('.molQTL','.gwas'),all.x=TRUE)
coloc_results_tableS1c <- merge(coloc_results_tableS1c, unique(sumStats[, .(best_molQTL_variantID = hm_variant_id, best_molQTL_rsID = hm_rsid)]), by = c("best_molQTL_rsID"), all.x = TRUE)
coloc_results_tableS1c <- merge(coloc_results_tableS1c, unique(sumStats[, .(best_GWAS_variantID = hm_variant_id, best_GWAS_rsID = hm_rsid)]), by = c("best_GWAS_rsID"), all.x = TRUE)
coloc_results_tableS1c <- coloc_results_tableS1c[order(factor(condition, c("NS", "LPS", "PAM3CSK4", "R848", "IAV")), gwasID, molQTL_phenoID), .(gwasID, molQTL_phenoID, condition, best_GWAS_variantID, best_GWAS_rsID, best_molQTL_variantID, best_molQTL_rsID, numberColocalisingVariants, h0, h1, h2, h3, h4, colocalisationMethod, betaRatioSignAverage, pval.gwas, beta.gwas, pval.molQTL, beta.molQTL)]
fwrite(coloc_results_tableS1c[condition == "NS", ], file = sprintf("%s/06_sum_stats_covQTL/TableS1c_coloc_results_covQTL_NSonly.tsv", PROJECT_DIR), sep = "\t")

