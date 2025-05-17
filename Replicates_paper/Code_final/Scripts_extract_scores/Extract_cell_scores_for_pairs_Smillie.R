library(data.table)
library(dplyr)

# Read the cell metadata dataframe for Smillie
cell_df = fread("/scratch/project_mnt/S0007/uqali4/2019_Smillie_normal_cellxgene_pc.cell_id_cell_type.csv")

# Process each dataset individually

# Crohn's Disease
dis_score_Crohns_disease_smillie = fread("/QRISdata/Q3999/tmp/scDRS_res_step1/2019_Smillie_normal_cellxgene_pc/mBAT/cd_build37_40266_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = index)

dis_score_Crohns_disease_smillie_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "activated_CD4_positive__alpha_beta_T_cell"),
  dis_score_Crohns_disease_smillie
) %>% mutate(
  Organ = "Colon",
  file = "cd_build37_40266_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz"
)

# Celiac Disease
dis_score_Celiac_disease_smillie = fread("/QRISdata/Q3999/tmp/scDRS_res_step1/2019_Smillie_normal_cellxgene_pc/mBAT/Celiac_disease_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = index)

dis_score_Celiac_disease_smillie_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "activated_CD4_positive__alpha_beta_T_cell"),
  dis_score_Celiac_disease_smillie
) %>% mutate(
  Organ = "Colon",
  file = "Celiac_disease_tidy.ma_10kb.full_score.gz"
)

# Inflammatory Bowel Disease (IBD)
dis_score_IBD_smillie = fread("/QRISdata/Q3999/tmp/scDRS_res_step1/2019_Smillie_normal_cellxgene_pc/mBAT/ibd_build37_59957_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = index)

dis_score_IBD_smillie_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "activated_CD4_positive__alpha_beta_T_cell"),
  dis_score_IBD_smillie
) %>% mutate(
  Organ = "Colon",
  file = "ibd_build37_59957_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz"
)

# Lymphocyte Count
dis_score_Lymphocyte_count_smillie = fread("/QRISdata/Q3999/tmp/scDRS_res_step1/2019_Smillie_normal_cellxgene_pc/mBAT/30120_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = index)

dis_score_Lymphocyte_count_smillie_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "activated_CD4_positive__alpha_beta_T_cell"),
  dis_score_Lymphocyte_count_smillie
) %>% mutate(
  Organ = "Colon",
  file = "30120_tidy.ma_10kb.full_score.gz"
)

# Now, merge all the individual dataframes into a single dataframe
list_of_dfs = list(
  dis_score_Crohns_disease_smillie_cd4_positive__alpha_beta_t_cell,
  dis_score_Celiac_disease_smillie_cd4_positive__alpha_beta_t_cell,
  dis_score_IBD_smillie_cd4_positive__alpha_beta_t_cell,
  dis_score_Lymphocyte_count_smillie_cd4_positive__alpha_beta_t_cell
)

# Merge all dataframes using bind_rows
merged_df_smillie = bind_rows(list_of_dfs)


# Save the merged dataframe to a gzipped CSV file
fwrite(merged_df_smillie, "/scratch/project_mnt/S0007/uqali4/NM_dt_results/merged_disease_scores_smillie.csv.gz", row.names = FALSE, compress = "gzip")
