this folder should contain 
- a list of GWAS association for ST8SIA4 (downloaded from GWAS catalog)
    ./data/01_GWAS_results/01_genes_ST8SIA4-associations-2022-12-2.csv
- metadata on GWAS studies from GWAS catalog (used for annotation purpose)
    ./data/01_GWAS_results/02_gwas-catalog-v1.0.3.1-studies-r2025-04-01.tsv
- full summary statistics from studies GCST90002316, GCST90002340 and GCST011096.
   ./data/01_GWAS_results/03_GWAS_summary_statistics/
- full summary statistics from studies GCST90002316, GCST90002340 and GCST011096.
   ./data/01_GWAS_results/03_GWAS_summary_statistics/
    These are not provided and should be downladed from GWAS catalog
- parquet files containing credible sets and colocalization results,
These are not provided, but can be downloaded from:
  http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.03/output/credible_set
  http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.03/output/coloc
- a list of SNPs overlapping with SuSiE 95% credibles sets from GWAS studies GCST90002316, GCST90002340 and GCST011096.
		this can be generated from script 00
		./data/01_GWAS_results/05_CS_selected_GWAS_traits.tsv

