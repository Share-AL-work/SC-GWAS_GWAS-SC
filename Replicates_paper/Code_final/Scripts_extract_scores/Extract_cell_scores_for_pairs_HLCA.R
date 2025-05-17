library(data.table)
library(dplyr)

# Read the cell metadata dataframe for HLCA
cell_df = fread("/scratch/project_mnt/S0007/uqali4/HLCA_core_healthy_LP_pc.cell_id_cell_type.csv")

# Process each dataset individually

# FEV1
dis_score_FEV1_HLCA = fread("/QRISdata/Q5729/scDRS_res_step1/HLCA_core_healthy_pconly/LP/FEV1_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = V1)

dis_score_FEV1_HLCA_bronchial_smooth_muscle_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "tracheobronchial smooth muscle cell"),
  dis_score_FEV1_HLCA,
  by = "cell_id"
) %>% mutate(
  Organ = "Lung",
  file = "FEV1_tidy.ma_10kb.full_score.gz"
)

# Celiac Disease
dis_score_Celiac_disease_HLCA = fread("/QRISdata/Q5729/scDRS_res_step1/HLCA_core_healthy_pconly/LP/Celiac_disease_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = V1)

dis_score_Celiac_disease_HLCA_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "CD4-positive, alpha-beta T cell"),
  dis_score_Celiac_disease_HLCA,
  by = "cell_id"
) %>% mutate(
  Organ = "Lung",
  file = "Celiac_disease_tidy.ma_10kb.full_score.gz"
)

# Crohn's Disease
dis_score_Crohns_disease_HLCA = fread("/QRISdata/Q5729/scDRS_res_step1/HLCA_core_healthy_pconly/LP/cd_build37_40266_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = V1)

dis_score_Crohns_disease_HLCA_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "CD4-positive, alpha-beta T cell"),
  dis_score_Crohns_disease_HLCA,
  by = "cell_id"
) %>% mutate(
  Organ = "Lung",
  file = "cd_build37_40266_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz"
)

# Inflammatory Bowel Disease (IBD)
dis_score_IBD_HLCA = fread("/QRISdata/Q5729/scDRS_res_step1/HLCA_core_healthy_pconly/LP/ibd_build37_59957_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = V1)

dis_score_IBD_HLCA_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "CD4-positive, alpha-beta T cell"),
  dis_score_IBD_HLCA,
  by = "cell_id"
) %>% mutate(
  Organ = "Lung",
  file = "ibd_build37_59957_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz"
)

# Lymphocyte Count
dis_score_Lymphocyte_count_HLCA = fread("/QRISdata/Q5729/scDRS_res_step1/HLCA_core_healthy_pconly/LP/30120_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = V1)

dis_score_Lymphocyte_count_HLCA_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "CD4-positive, alpha-beta T cell"),
  dis_score_Lymphocyte_count_HLCA,
  by = "cell_id"
) %>% mutate(
  Organ = "Lung",
  file = "30120_tidy.ma_10kb.full_score.gz"
)

# Vitiligo
dis_score_Vitiligo_HLCA = fread("/QRISdata/Q5729/scDRS_res_step1/HLCA_core_healthy_pconly/LP/Vitiligo_Jin_2016_NG_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = V1)

dis_score_Vitiligo_HLCA_cd8_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "CD8-positive, alpha-beta T cell"),
  dis_score_Vitiligo_HLCA,
  by = "cell_id"
) %>% mutate(
  Organ = "Lung",
  file = "Vitiligo_Jin_2016_NG_tidy.ma_10kb.full_score.gz"
)

# Now, merge all the individual dataframes into a single dataframe
list_of_dfs = list(
  dis_score_FEV1_HLCA_bronchial_smooth_muscle_cell,
  dis_score_Celiac_disease_HLCA_cd4_positive__alpha_beta_t_cell,
  dis_score_Crohns_disease_HLCA_cd4_positive__alpha_beta_t_cell,
  dis_score_IBD_HLCA_cd4_positive__alpha_beta_t_cell,
  dis_score_Lymphocyte_count_HLCA_cd4_positive__alpha_beta_t_cell,
  dis_score_Vitiligo_HLCA_cd8_positive__alpha_beta_t_cell
)

# Merge all dataframes using bind_rows
merged_df_HLCA = bind_rows(list_of_dfs)

# Save the merged dataframe to a gzipped CSV file
fwrite(merged_df_HLCA, "/scratch/project_mnt/S0007/uqali4/NM_dt_results/merged_disease_scores_HLCA.csv.gz", row.names = FALSE, compress = "gzip")









