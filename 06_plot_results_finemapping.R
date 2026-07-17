options(stringsAsFactors = FALSE, max.print = 9999, width = 300, datatable.fread.input.cmd.message = FALSE)


PROJECT_DIR <- "./ST8SIA4_GWAS_molQTL/"
FIGURE_DIR <- sprintf("%s/figures/", PROJECT_DIR)

CHR <- 5
START_ST8SIA4_b38 <- 100806933
END_ST8SIA4_b38 <- 100903282

START_ST8SIA4_b37 <- 100142637
END_ST8SIA4_b37 <- 100238986
CIS_DIST <- 1e6

geneID <- "ENSG00000113532"

library(data.table)
library(dplyr)
library(ggplot2)
library(tictoc)
library(coloc)
library(susieR)
library(snpStats)
library(corrplot)
library(rtracklayer)


###################################
##### load gwas Credible sets #####
###################################

credible_sets_GWAS <- fread(sprintf("%s/data/01_GWAS_results/05_CS_selected_GWAS_traits.tsv", PROJECT_DIR))

#####################################################
##### load cov_QTL data and fine-mapped variants ####
#####################################################

QTL_pvalues <- fread(sprintf("%s/data/06_sum_stats_covQTL/coverage_QTL_ST8SIA4_bestP_perFeature.txt.gz", PROJECT_DIR))
feature_signif <- fread(sprintf("%s/data/06_sum_stats_covQTL/coverage_QTL_ST8SIA4_perFeature_FDR.tsv", PROJECT_DIR))
covQTL_sumStats <- fread(sprintf("%s/data/06_sum_stats_covQTL/coverage_QTL_ST8SIA4_sumStats.txt.gz", PROJECT_DIR))

covQTL_sumStats <- merge(covQTL_sumStats, credible_sets_GWAS[, .(hm_rsid, is_CS_ID, hm_variant_id)], by.x = "snps", by.y = "hm_rsid", all.x = TRUE)
explained_components <- fread(sprintf("%s/data/06_sum_stats_covQTL/finemapped_covQTL_explained.tsv", PROJECT_DIR))
covQTL_sumStats_signif <- fread(sprintf("%s/data/06_sum_stats_covQTL/finemapped_covQTL_signif.tsv", PROJECT_DIR), sep = "\t")

#############################
### LOAD EXPRESSION DATA  ###
#############################

# sample annot
EGA_samples <- fread(sprintf("%s/data/02_fastq/EGA_RNA_seq_annotation.txt", PROJECT_DIR))

# junction data
junctions <- fread(sprintf("%s/data/04_coverage/aggregate/ST8SIA4_junctions_all_samples.tsv.gz", PROJECT_DIR))
junctions <- merge(junctions, EGA_samples, by.x = "sample", by.y = "ID")
junctions[, IID := substr(sample_ID, 1, 6)]
# coverage data
coverage_bins <- fread(sprintf("%s/data/04_coverage/aggregate/ST8SIA4_coverage_bins_all_samples.tsv.gz", PROJECT_DIR))
coverage_bins <- merge(coverage_bins, EGA_samples, by.x = "sample", by.y = "ID")
coverage_bins[, IID := substr(sample_ID, 1, 6)]

# annotate coverage bins  (start, end , name)
bin_annot <- unique(coverage_bins[, .(bin, chr = 5, start = x_min, end = x_max)])
bin_annot <- merge(bin_annot, coverage_bins[, .(coverage = mean(relative_coverage)), keyby = .(bin, bin_x = x)], by = "bin")
bin_annot <- bin_annot[order(bin_x), ]
bin_annot[, bin_label := paste(format(floor(start / 1000), big.mark = ",", scientific = FALSE), "-", format(ceiling(end / 1000), big.mark = ",", scientific = FALSE), "kb")]
bin_annot[coverage > 0.01, ]

#############################
##### load genotypes    #####
#############################

  # Genotype_path_chr5_200ind_QuachEtAl <- sprintf("%s/data/05_genotypes/EvoImmunoPop_imputation_200x19619457_chr5",PROJECT_DIR)
  Genotype_path_ST8SIA4_200ind_QuachEtAl <- sprintf("%s/data/05_genotypes/EvoImmunoPop_imputation_200x19619457_ST8SIA4", PROJECT_DIR)

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
  genotypes_long <- as.data.table(melt(genotypes, measure.vars = colnames(genotypes), value.name = "genotype"))
  setnames(genotypes_long, c("Var1", "Var2", "genotype"), c("IID", "rsID", "genotype"))
  # Access SNP map data
  snp_data <- as.data.table(plink_data$map)
  print(snp_data)
  setnames(snp_data, c("snp.name", "allele.1", "allele.2"), c("rsID", "effect_allele", "ref_allele"))
  genotypes_long <- merge(genotypes_long, snp_data[, .(rsID, POS_b37 = position, effect_allele, ref_allele)], by = "rsID", all.x = TRUE)

#############################
### LOAD gene structure  ####
#############################
# ST8SIA4_model <- fread("grep -e ENSG00000113532 ./data/07_ST8SIA4_gene_structure/genes_hg38.gtf")
# fwrite(ST8SIA4_model,file=sprintf('%s/data/07_ST8SIA4_gene_structure/ST8SIA4_gene_model_hg38.gtf',PROJECT_DIR),sep='\t')

ST8SIA4_model=fread(sprintf('%s/data/07_ST8SIA4_gene_structure/ST8SIA4_gene_model_hg38.gtf',PROJECT_DIR))
exon_ST8SIA4 <- ST8SIA4_model[V3 == "exon", ]
exon_ST8SIA4[, transcript := str_split(V9, "; ", simplify = T)[, 3]]
exon_ST8SIA4[, transcript_name := str_split(V9, "; ", simplify = T)[, 8]]
exon_ST8SIA4[, exon_nb := str_split(V9, "; ", simplify = T)[, 9]]
exon_ST8SIA4[, transcript := gsub('transcript_id "(.*)"', "\\1", transcript)]
exon_ST8SIA4[, transcript_name := gsub('transcript_name "(.*)"', "\\1", transcript_name)]
exon_ST8SIA4[, exon_nb := gsub('exon_number "(.*)"', "\\1", exon_nb)]
exon_ST8SIA4[, V9 := NULL]
exon_ST8SIA4[, transcript_Nb := cumsum(!duplicated(transcript))]

#pos_UTRs <- fread("/data/07_ST8SIA4_gene_structure/utrome.e30.t5.gc39.pas3.f0.9999.w500.gtf")[grepl("ST8SIA4", V9) & V3 == "exon", ]
## downloaded from https://github.com/Mayrlab/hcl-utrome/releases/download/v1.0.0/utrome.e30.t5.gc39.pas3.f0.9999.w500.tar.gz
# fwrite(pos_UTRs,file=sprintf('%s/data/07_ST8SIA4_gene_structure/ST8SIA4_UTR_hg38.gtf',PROJECT_DIR),sep='\t')

pos_UTRs <- fread("data/07_ST8SIA4_gene_structure/ST8SIA4_UTR_hg38.gtf")
#polyA_db <- fread(sprintf("%s/TableSX2_polyA_sites.csv", PROJECT_DIR))[!is.na(PeakSite)]

# polyA_ID_minus <- import(sprintf("%s/results/polyA/putative_polya_sites.-.bb", PROJECT_DIR))
# polyA_ST8SIA4 <- polyA_ID_minus[grep("ST8SIA4", polyA_ID_minus$name), ]
# polyA_ST8SIA4 <- as.data.table(polyA_ST8SIA4)
# polyA_ST8SIA4[, type := str_split(name, "\\|", simplify = T)[, 2]]
# polyA_ST8SIA4[, polyaIDclassification := as.numeric(gsub("polyaIDclassification_", "", str_split(name, "\\|", simplify = T)[, 3]))]
# polyA_ST8SIA4[, polyaIDcleavageprofile := as.numeric(gsub("polyaIDcleavageprofile_", "", str_split(name, "\\|", simplify = T)[, 4]))]
# polyA_ST8SIA4[, reads := as.numeric(gsub("reads_", "", str_split(name, "\\|", simplify = T)[, 5]))]
# polyA_ST8SIA4[, ru := as.numeric(gsub("ru_", "", str_split(name, "\\|", simplify = T)[, 6]))]
# polyA_ST8SIA4[reads >= 3, ]
# fwrite(polyA_ST8SIA4, file = sprintf("%s/data/07_ST8SIA4_gene_structure/polyA_sites_ST8SIA4.tsv", PROJECT_DIR), sep = "\t")
polyA_ST8SIA4 <- fread(sprintf("%s/data/07_ST8SIA4_gene_structure/polyA_sites_ST8SIA4.tsv", PROJECT_DIR))

# extract junctions from ST8SIA4 gene 
junctions_ST8SIA4 <- junctions[, .(count_per_sample = sum(count) / length(unique(sample))), by = .(x, xend)][order(-count_per_sample)][x > START_ST8SIA4_b38 - 1e4 & xend < END_ST8SIA4_b38 + 1e4, ]
junctions_ST8SIA4 <- junctions_ST8SIA4[order(x), ]
# define y axis position to avoid overlap
junctions_ST8SIA4[, dup_x := 1:.N, by = x]
junctions_ST8SIA4[10:12, dup_x := c(3, 2, 5)]
junctions_ST8SIA4[16, dup_x := 2]
junctions_ST8SIA4[16, dup_x := 2]
junctions_ST8SIA4[18:21, dup_x := c(3:5, 5)]

################################################################################
########## plot coverage for target SNP by condition (Fig 1B) ##################
################################################################################

mySNP <- "rs2548257"
myCOND <- "NS"
myREF <- "C"
myALT <- "A"
mySNP_POS <- 100826613

# plotting all bins and condition
Fig_data <- merge(coverage_bins, genotypes_long[rsID == mySNP], by = c("IID"))
Fig_data[, genotype_char := case_when(
  genotype == 0 ~ paste(myREF, myREF, sep = "/"),
  genotype == 1 ~ paste(myALT, myREF, sep = "/"),
  genotype == 2 ~ paste(myALT, myALT, sep = "/")
)]

Fig_data_agg <- Fig_data[, .(relative_coverage = mean(relative_coverage)), keyby = .(x, condition, genotype_char)]
levels_inverted <- c(paste(myALT, myALT, sep = "/"), paste(myALT, myREF, sep = "/"), paste(myREF, myREF, sep = "/"))

p <- ggplot(Fig_data_agg[condition == "NS", ]) + ylim(-0.12, 0.38)
p <- p + geom_bar(aes(x = x, y = relative_coverage, fill = factor(genotype_char, levels_inverted)), stat = "identity", position = "dodge") + theme_plot(lpos = "right")
p <- p + geom_segment(data = exon_ST8SIA4[, .(V4 = min(V4), V5 = max(V5)), by = .(transcript_Nb)], aes(x = V4, xend = V5, y = -0.02 - 0.06 / 5 * transcript_Nb))
p <- p + geom_rect(data = exon_ST8SIA4, aes(xmin = V4, xmax = V5, ymin = -0.02 - 0.06 / 5 * transcript_Nb - 0.005, ymax = -0.02 - 0.06 / 5 * transcript_Nb + 0.005))
p <- p + xlab("") + ylab("Relative coverage") + theme(text = element_text(size = 12, family = "sans"))

p <- p + geom_point(data = polyA_ST8SIA4[reads >= 10], aes(x = start, y = -0.005), col = "red", size = 0.5)
p <- p + geom_point(data = pos_UTRs, aes(x = V4, y = -0.005), size = 0.7)
p <- p + geom_segment(data = junctions_ST8SIA4, aes(x = x, xend = xend, y = -0.12 + 0.05 / 7 * dup_x, col = log10(count_per_sample)))

pdf(sprintf("%s/Fig1B__covQTL_%s__%s_%s_02_coverage_ST8SIA4_NSonly.pdf", FIGURE_DIR, i, mySNP, myCOND), height = 4, width = 7)
print(p)
dev.off()

#############################################################################
########## plot localization of coverage QTLs (Fig 1C) ######################
#############################################################################

Fig_Data <- merge(covQTL_sumStats_signif, credible_sets_GWAS, by.x = c("snps", "is_CS_ID", "hm_variant_id"), by.y = c("variant_id", "is_CS_ID", "hm_variant_id"))
Fig_Data <- merge(Fig_Data, bin_annot, by.x = "gene", by.y = "bin")
Fig_Data[, .N, by = .(ref_allele == hm_effect_allele & effect_allele == hm_other_allele, ref_allele == hm_other_allele & effect_allele == hm_effect_allele)]
Fig_Data[ref_allele == hm_effect_allele & effect_allele == hm_other_allele, statistic := -statistic]
Fig_Data[ref_allele == hm_effect_allele & effect_allele == hm_other_allele, beta := -beta]
Fig_Data[, ref_allele := NULL]
Fig_Data[, effect_allele := NULL]

Fig_Data[, logP := -log10(pvalue)]
Fig_Data[, sign := sign(statistic)]

bin_scale <- c("#067BC2", "#84BCDA", "#ADBE87", "#ecb40b", "#F37748", "#c94b4d", "#732728")
names(bin_scale) <- unique(Fig_Data[condition == "NS", ][order(bin_x), gsub("coverage ", "", bin_label)])

p <- ggplot(Fig_Data[condition == "NS", ], aes(x = hm_pos / 1000, y = PIP, col = gsub("coverage ", "", bin_label), fill = gsub("coverage ", "", bin_label), shape = as.character(sign)), alpha = ifelse(credible_set != "", 1, .1)) +
  rasterize(geom_point(), dpi = 200)
# p <- p + geom_rect(data = exon_ST8SIA4, aes(xmin = V4, xmax = V5, ymin = -0.02 - 0.06 / 5 * transcript_Nb - 0.005, ymax = -0.02 - 0.06 / 5 * transcript_Nb + 0.005))
p <- p + scale_shape_manual(values = c("-1" = 25, "1" = 24)) + scale_fill_manual(values = bin_scale) + scale_color_manual(values = bin_scale)
p <- p + xlab("position (kb)") + ylab("PIP") + theme_plot(lpos = "right") + theme(text = element_text(size = 9, family = "sans"))
p <- p + theme(legend.spacing.y = unit(0.3, "cm")) + guides(colour = guide_legend(override.aes = list(alpha = 1), byrow = TRUE), fill = "none", shape = "none")
p <- p + theme(legend.key.size = unit(.1, "cm"))

pdf(sprintf("%s/Fig1C__covQTL_localization_ST8SIA4_PIP.pdf", FIGURE_DIR), width = 5.2, height = 2.2)
print(p)
dev.off()


#############################################################################
########## plot colocalization of coverage QTLs with GWAS (Fig 1D) ##########
#############################################################################

coloc_results_tableS3 <- fread(file = sprintf("%s/TableS1c_coloc_results_covQTL_NSonly.tsv", PROJECT_DIR), sep = "\t")

coloc_matrix_wide <- dcast(coloc_results_tableS3, gwasID ~ molQTL_phenoID, value.var = "h4")
coloc_matrix <- as.matrix(coloc_matrix_wide[, -1])
rownames(coloc_matrix) <- coloc_matrix_wide[, gwasID]
colnames(coloc_matrix) <- gsub("coverage ", "", colnames(coloc_matrix))

pdf(sprintf("%s/Fig1D__coloc_covQTL_GWAS.pdf", FIGURE_DIR), height = 4, width = 7)
corrplot(coloc_matrix[c(2, 3, 1), ], col.lim = c(0, 1), tl.col = "black", tl.cex = 1.5, cl.cex = 1.5, cl.length = 6, cl.align.text = "l")
dev.off()
