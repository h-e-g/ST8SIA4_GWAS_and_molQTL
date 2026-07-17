PROJECT_DIR <- "./ST8SIA4_GWAS_molQTL/"
OPENTARGETS_DIR <- sprintf("%s/data/01_GWAS_results/02_openTargets_25_03/", PROJECT_DIR)

chr <- 5
start_b38 <- 100806933
end_b38 <- 100903282


library(data.table)
library(ggplot2)
library(arrow)
library(tictoc)
library(janitor)

# load GWAS catalog meta data on traits
GWAS_catalog_traits <- fread(sprintf("%s/data/01_GWAS_results/02_gwas-catalog-v1.0.3.1-studies-r2025-04-01.tsv", PROJECT_DIR))
GWAS_catalog_traits <- clean_names(GWAS_catalog_traits)

# define target GWAS ids to plot and add trait name
GWAS_ids <- c("GCST90002316", "GCST90002340", "GCST011096")
name_corresp <- GWAS_catalog_traits[match(GWAS_ids, study_accession), paste(GWAS_ids, mapped_trait, sep = " - ")]
names(name_corresp) <- GWAS_ids

###########################################################################
#####################  Extract summary stats ##############################
###########################################################################

FILES <- dir(sprintf("%s/data/01_GWAS_results/03_GWAS_summary_statistics", PROJECT_DIR), pattern = ".gz")
sumStats <- list()
for (FILE in FILES) {
  if (FILE != "GCST011096.h.tsv.gz") {
    sumStats[[FILE]] <- fread(sprintf("%s/SummStats_GWAS/%s", PROJECT_DIR, FILE))[hm_chrom == 5 & hm_pos > 100.8e6 - 1e6 & hm_pos < 100.8e6 + 1e6, ]
  } else {
    sumStats_FILE <- fread(sprintf("%s/SummStats_GWAS/%s", PROJECT_DIR, FILE))
    sumStats_FILE <- sumStats_FILE[, .(variant_id, rsid, chromosome, base_pair_location, other_allele, effect_allele, beta, odds_ratio, ci_lower, ci_upper, effect_allele_frequency, hm_code, chromosome, variant_id, base_pair_location, effect_allele, other_allele, p_value, ci_lower, standard_error, effect_allele_frequency, odds_ratio, beta, ci_upper)]
    colnames(sumStats_FILE) <- c("hm_variant_id", "hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele", "hm_beta", "hm_odds_ratio", "hm_ci_lower", "hm_ci_upper", "hm_effect_allele_frequency", "hm_code", "chromosome", "variant_id", "base_pair_location", "effect_allele", "other_allele", "p_value", "ci_lower", "standard_error", "effect_allele_frequency", "odds_ratio", "beta", "ci_upper")
    sumStats_FILE[, hm_variant_id := paste(chromosome, base_pair_location, other_allele, effect_allele, sep = "_")]
    sumStats[[FILE]] <- sumStats_FILE[hm_chrom == 5 & hm_pos > 100.8e6 - 1e6 & hm_pos < 100.8e6 + 1e6, ]
  }
}
ii <- intersect(colnames(sumStats[[1]]), colnames(sumStats[[2]]))

sumStats <- rbindlist(lapply(sumStats, function(x) {
  x[, mget(ii)]
}), idcol = "studyId")
sumStats[studyId == "GCST011096.h.tsv.gz", studyId := "XXXXXXXX-GCST011096-EFO_000XXXX.h.tsv.gz"]

sumStats[, PMID := gsub("([X0-9]+)-(GCST[0-9]+)-(EFO_[X0-9]+).h.tsv.gz", "\\1", studyId)]
sumStats[, EFO_number := gsub("([X0-9]+)-(GCST[0-9]+)-(EFO_[X0-9]+).h.tsv.gz", "\\3", studyId)]
sumStats[, studyId := gsub("([X0-9]+)-(GCST[0-9]+)-(EFO_[X0-9]+).h.tsv.gz", "\\2", studyId)]
sumStats[, p_value := as.numeric(p_value)]

sumStats[, full_name := paste0(GWAS_catalog_traits[match(studyId, study_accession), mapped_trait], "\n(", studyId, ")")]
sumStats[, is_CS := paste(studyId, hm_variant_id) %chin% ST8SIA4_finemap[, paste(studyId, variantId)]]
fwrite(sumStats, file = sprintf("%s/data/01_GWAS_results/03_GWAS_summary_statistics/01_sumStats_selected_traits_ST8SIA4_locus.tsv", PROJECT_DIR), sep = "\t")
# sumStats=fread(sprintf("%s/data/01_GWAS_results/03_GWAS_summary_statistics/01_sumStats_selected_traits_ST8SIA4_locus.tsv", PROJECT_DIR))

CS_selected_GWAS_traits=sumStats[order(p_value),.(sumLogP=sum(-log10(p_value)),is_CS_any=sum(is_CS),is_CS_ID=paste(studyId[is_CS],collapse='/')),by=.(variant_id,hm_variant_id,hm_rsid , hm_chrom ,hm_pos, hm_other_allele, hm_effect_allele)][order(-sumLogP),]
fwrite(CS_selected_GWAS_traits,file=sprintf('%s/data/01_GWAS_results/05_CS_selected_GWAS_traits.tsv',PROJECT_DIR),sep='\t')
# CS_selected_GWAS_traits=fread(sprintf('%s/data/01_GWAS_results/05_CS_selected_GWAS_traits.tsv',PROJECT_DIR))
###########################################################################
#####################  Extract credible sets ##############################
###########################################################################

# read credible sets results from OpenTargets
list_CS_parquet_files <- dir(sprintf("%s/credibleSets", OPENTARGETS_DIR), pattern = "parquet")
CS_results <- list()
for (FILE in list_CS_parquet_files) {
  tic(FILE)
  CS_results[[FILE]] <- as.data.table(read_parquet(sprintf("%s/credibleSets/%s", OPENTARGETS_DIR, FILE)))
  toc()
}
CS_results <- rbindlist(CS_results)

### extract all credible sets at the ST8SIA4 locus
ST8SIA4_results <- CS_results[chromosome == chr & position > (start_b38 - 1e4) & position < (end_b38 + 1e4), ]
ST8SIA4_assoc_list <- ST8SIA4_results[, .(studyType, studyLocusId, studyId, variantId)]
ST8SIA4_gwas_list <- ST8SIA4_results[studyType == "gwas", studyLocusId]
ST8SIA4_molQTL_list <- ST8SIA4_results[studyType != "gwas", studyLocusId]

#### focus on the 3 GWAS traits we selected.
ST8SIA4_finemap <- ST8SIA4_results[!grepl("PICS fine-mapped credible set", confidence), locus[[1]], by = .(studyLocusId, studyId, studyType, chromosome, position, bestVariant = variantId, subStudyDescription, finemappingMethod)]
ST8SIA4_finemap <- ST8SIA4_finemap[studyType != "gwas" | studyId %in% GWAS_ids, ]
ST8SIA4_finemap[, full_name := studyId]
ST8SIA4_finemap[studyType == "gwas", full_name := name_corresp[studyId]]
ST8SIA4_finemap[, Zval := beta / standardError]
ST8SIA4_finemap[, direction := sign(beta) / 2]
fwrite(ST8SIA4_finemap, file = sprintf("%s/data/01_GWAS_results/03_GWAS_summary_statistics/02_ST8SIA4_finemapping_selected_traits.tsv", PROJECT_DIR), sep = "\t")
#ST8SIA4_finemap <- fread(sprintf("%s/data/01_GWAS_results/03_GWAS_summary_statistics/02_ST8SIA4_finemapping_selected_traits.tsv", PROJECT_DIR))

###########################################################################
#####################  plot utilities  ####################################
###########################################################################

# custom plot theme
theme_plot <- function(lpos = "bottom", rotate.x = 0) {
  ptheme <- theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      text = element_text(size = 7, family = "sans"),
      legend.position = lpos,
      legend.background = element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      strip.background = element_rect(fill = "transparent", colour = NA),
      panel.spacing = unit(0, "pt"),
      strip.text = element_text(color = "black")
    )
  if (rotate.x > 0) {
    ptheme <- ptheme + theme(axis.text.x = element_text(angle = rotate.x, hjust = 1, vjust = 0.5))
  }
  ptheme
}


#####################################################################################
#####################################################################################





#####################################################################################
##################### .  Figure 1a (top panel) #######################################
#####################################################################################

# reload SumStats (if needed)
#sumStats <- fread( sprintf("%s/data/01_GWAS_results/03_GWAS_summary_statistics/01_sumStats_selected_traits_ST8SIA4_locus.tsv", PROJECT_DIR), sep = "\t")
FigData <- sumStats[hm_chrom == 5 & hm_pos > start_b38 - 5e5 & hm_pos < end_b38 + 5e5 & studyId %in% c("GCST90002316", "GCST90002340", "GCST011096"), ]
FigData[, full_name := gsub("lupus ", "lupus\n", full_name)]
fwrite(FigData, file = sprintf("%s/final_figures/Fig1A_top_GWAS_locusZoom_plot_3TRAITS_source_Data.tsv", PROJECT_DIR), sep = "\t")

pdf(sprintf("%s/final_figures/Fig1A_top_GWAS_locusZoom_plot_3TRAITS.pdf", PROJECT_DIR), height = 4, width = 4)
p <- ggplot(FigData, aes(x = base_pair_location / 1e6, y = -log10(p_value), col = is_CS, alpha = is_CS, size = is_CS)) +
  theme_plot() +
  geom_hline(yintercept = 0, col = "darkgrey")
p <- p + geom_point() + facet_grid(rows = vars(full_name)) + guides(color = "none", size = "none", alpha = "none") + scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"))
p <- p + theme(strip.text.x = element_text(size = 7)) + scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1)) + scale_size_manual(values = c("FALSE" = 0.1, "TRUE" = 1))
p <- p + geom_rect(xmin = start_b38 / 1e6, xmax = end_b38 / 1e6, ymin = -1, ymax = -0.5, fill = "darkgrey", col = "black", alpha = 0.1)
p <- p + ylim(c(-1, 12)) + xlab("Chromosome 5 position (Mb)") + geom_hline(yintercept = -log10(5e-8), col = "lightgrey", lty = 2) + geom_text(x = end_b38 / 1e6 + 0.02, y = -0.75, label = "ST8SIA4", col = "black", hjust = 0, size = 2, fontface = "italic")
print(p)
dev.off()

#####################################################################################
#####################################################################################



#####################################################################################
##################### .  Figure 1a (bottom panel) ####################################
#####################################################################################

# plot Credible sets a the locus (Fig1 A)
library(corrplot)
ST8SIA4_gwas <- dcast(ST8SIA4_finemap[studyType == "gwas", ], full_name ~ variantId, value.var = "direction")
matrix_gwas <- as.matrix(ST8SIA4_gwas[, -1])
rownames(matrix_gwas) <- ST8SIA4_gwas[, full_name]
matrix_gwas[is.na(matrix_gwas)] <- 0
fwrite(data.table(gwas_id=rownames(matrix_gwas),matrix_gwas), file = sprintf("%s/final_figures/Fig1A_bottom_gwas_finemap_allTraits_ordered__source_data.tsv", PROJECT_DIR), sep='\t')

pdf(sprintf("%s/final_figures/Fig1a_bottom_gwas_finemap_allTraits_ordered.pdf", PROJECT_DIR), width = 12)
hr <- hclust(as.dist(1 - cor(t(abs(matrix_gwas)))))
# hc=hclust(as.dist(1-cor(abs(Z_matrix))))
corrplot(matrix_gwas[c(2, 3, 1), ], is.corr = FALSE, tl.cex = 0.6, tl.col = "black", na.label = "square", na.label.col = "white", col.lim = c(-.8, 0.8))
dev.off()

#####################################################################################
#####################################################################################
