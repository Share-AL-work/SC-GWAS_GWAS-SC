library(data.table)
library(dplyr)

# Read the cell metadata dataframe for FBM
cell_df = fread("/scratch/project_mnt/S0007/uqali4/FBM_Jardine_2021_Nature_pc.cell_id_cell_type.csv")

# Process each dataset individually
# Platelet count
dis_score_PltC_FBM = fread("/QRISdata/Q3999/tmp/scDRS_res_step1/FBM_Jardine_2021_Nature_pc/mBAT/30080_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = index)

dis_score_PltC_FBM = inner_join(
  cell_df %>% filter(cell_ontology_class == "megakaryocyte"),
  dis_score_Red_blood_cell_count_FBM
) %>% mutate(
  Organ = "Bone Marrow and blood",
  file = "30080_tidy.ma_10kb.full_score.gz"
)


# Red Blood Cell Count
dis_score_Red_blood_cell_count_FBM = fread("/QRISdata/Q3999/tmp/scDRS_res_step1/FBM_Jardine_2021_Nature_pc/mBAT/30010_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = index)

dis_score_Red_blood_cell_count_FBM_erythrocyte = inner_join(
  cell_df %>% filter(cell_ontology_class == "erythrocyte"),
  dis_score_Red_blood_cell_count_FBM
) %>% mutate(
  Organ = "Bone Marrow and blood",
  file = "30010_tidy.ma_10kb.full_score.gz"
)

# Inflammatory Bowel Disease (IBD)
dis_score_IBD_FBM = fread("/QRISdata/Q3999/tmp/scDRS_res_step1/FBM_Jardine_2021_Nature_pc/mBAT/ibd_build37_59957_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = index)

dis_score_IBD_FBM_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "CD4_positive__alpha_beta_T_cell"),
  dis_score_IBD_FBM
) %>% mutate(
  Organ = "Bone Marrow and blood",
  file = "ibd_build37_59957_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz"
)

# Celiac Disease
dis_score_Celiac_disease_FBM = fread("/QRISdata/Q3999/tmp/scDRS_res_step1/FBM_Jardine_2021_Nature_pc/mBAT/Celiac_disease_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = index)

dis_score_Celiac_disease_FBM_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "CD4_positive__alpha_beta_T_cell"),
  dis_score_Celiac_disease_FBM
) %>% mutate(
  Organ = "Bone Marrow and blood",
  file = "Celiac_disease_tidy.ma_10kb.full_score.gz"
)

# Crohn's Disease
dis_score_Crohns_disease_FBM = fread("/QRISdata/Q3999/tmp/scDRS_res_step1/FBM_Jardine_2021_Nature_pc/mBAT/cd_build37_40266_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = index)

dis_score_Crohns_disease_FBM_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "CD4_positive__alpha_beta_T_cell"),
  dis_score_Crohns_disease_FBM
) %>% mutate(
  Organ = "Bone Marrow and blood",
  file = "cd_build37_40266_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz"
)

# Vitiligo
dis_score_Vitiligo_FBM = fread("/QRISdata/Q3999/tmp/scDRS_res_step1/FBM_Jardine_2021_Nature_pc/mBAT/Vitiligo_Jin_2016_NG_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = index)

dis_score_Vitiligo_FBM_cd8_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "CD8_positive__alpha_beta_T_cell"),
  dis_score_Vitiligo_FBM
) %>% mutate(
  Organ = "Bone Marrow and blood",
  file = "Vitiligo_Jin_2016_NG_tidy.ma_10kb.full_score.gz"
)

# Lymphocyte Count
dis_score_Lymphocyte_count_FBM = fread("/QRISdata/Q3999/tmp/scDRS_res_step1/FBM_Jardine_2021_Nature_pc/mBAT/30120_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = index)

dis_score_Lymphocyte_count_FBM_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "CD4_positive__alpha_beta_T_cell"),
  dis_score_Lymphocyte_count_FBM
) %>% mutate(
  Organ = "Bone Marrow and blood",
  file = "30120_tidy.ma_10kb.full_score.gz"
)

# Now, merge all the individual dataframes into a single dataframe
list_of_dfs = list(
  dis_score_Red_blood_cell_count_FBM_erythrocyte,
  dis_score_IBD_FBM_cd4_positive__alpha_beta_t_cell,
  dis_score_Celiac_disease_FBM_cd4_positive__alpha_beta_t_cell,
  dis_score_Crohns_disease_FBM_cd4_positive__alpha_beta_t_cell,
  dis_score_Vitiligo_FBM_cd8_positive__alpha_beta_t_cell,
  dis_score_Lymphocyte_count_FBM_cd4_positive__alpha_beta_t_cell,
  dis_score_PltC_FBM
)

# Merge all dataframes using bind_rows
merged_df_FBM = bind_rows(list_of_dfs)

# Save the merged dataframe to a gzipped CSV file
fwrite(merged_df_FBM, "/scratch/project_mnt/S0007/uqali4/NM_dt_results/merged_disease_scores_FBM.csv.gz", row.names = FALSE, compress = "gzip")
































