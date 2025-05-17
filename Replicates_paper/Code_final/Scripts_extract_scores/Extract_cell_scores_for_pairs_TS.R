library(data.table)
library(dplyr)

# Read the cell metadata dataframe for Tabula Sapiens
cell_df = fread("/scratch/project_mnt/S0007/uqali4/TabulaSapiens_pc_ortholog_with_TMS_minCell20.cell_id_cell_ontology_class.csv")

# Process each dataset individually

# 1. Red Blood Cell Count
dis_score_plt_count_TS = fread("/QRISdata/Q5059/scDRS_res_step1/TabulaSapiens_pc_ortholog/30080_tidy.ma_10kb.full_score.gz")  # 475392 cells
dis_score_plt_count_TS_platelet = inner_join(
  cell_df %>% filter(cell_ontology_class == "platelet"),
  dis_score_plt_count_TS
) %>% mutate(
  Organ = "Atlas - Human",
  file = "30080_tidy.ma_10kb.full_score.gz"
)  # 233 cells

dis_score_Red_blood_cell_count_TS = fread("/QRISdata/Q5059/scDRS_res_step1/TabulaSapiens_pc_ortholog/30010_tidy.ma_10kb.full_score.gz")  # 475392 cells
dis_score_Red_blood_cell_count_TS_erythrocyte = inner_join(
  cell_df %>% filter(cell_ontology_class == "erythrocyte"),
  dis_score_Red_blood_cell_count_TS
) %>% mutate(
  Organ = "Atlas - Human",
  file = "30010_tidy.ma_10kb.full_score.gz"
) 

# 2. Systolic Blood Pressure (SBP)
dis_score_SBP_TS = fread("/QRISdata/Q5059/scDRS_res_step1/TabulaSapiens_pc_ortholog/Systolic_blood_pressure_tidy.ma_10kb.full_score.gz")  # 475392 cells
dis_score_SBP_TS_smooth_muscle_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "vascular_associated_smooth_muscle_cell"),
  dis_score_SBP_TS
) %>% mutate(
  Organ = "Atlas - Human",
  file = "Systolic_blood_pressure_tidy.ma_10kb.full_score.gz"
)  # 2669 cells

# 3. Celiac Disease
dis_score_Celiac_disease_TS = fread("/QRISdata/Q5059/scDRS_res_step1/TabulaSapiens_pc_ortholog/Celiac_disease_tidy.ma_10kb.full_score.gz")
dis_score_Celiac_disease_TS_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "cd4_positive__alpha_beta_t_cell"),
  dis_score_Celiac_disease_TS
) %>% mutate(
  Organ = "Atlas - Human",
  file = "Celiac_disease_tidy.ma_10kb.full_score.gz"
)  # 13921 cells

list_of_dfs = list(
  dis_score_Celiac_disease_TS_cd4_positive__alpha_beta_t_cell,
  dis_score_SBP_TS_smooth_muscle_cell,
  dis_score_Red_blood_cell_count_TS_erythrocyte
)

# Merge all dataframes using bind_rows
merged_df_TS = bind_rows(list_of_dfs)



# Crohn's Disease
dis_score_Crohns_disease_TS = fread("/QRISdata/Q5059/scDRS_res_step1/TabulaSapiens_pc_ortholog/cd_build37_40266_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz")
dis_score_Crohns_disease_TS_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "cd4_positive__alpha_beta_t_cell"),
  dis_score_Crohns_disease_TS
) %>% mutate(
  Organ = "Atlas - Human",
  file = "cd_build37_40266_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz"
)

# Inflammatory Bowel Disease (IBD)
dis_score_IBD_TS = fread("/QRISdata/Q5059/scDRS_res_step1/TabulaSapiens_pc_ortholog/ibd_build37_59957_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz")
dis_score_IBD_TS_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "cd4_positive__alpha_beta_t_cell"),
  dis_score_IBD_TS
) %>% mutate(
  Organ = "Atlas - Human",
  file = "ibd_build37_59957_20161107.txt_splitted.tab_tidy.ma_10kb.full_score.gz"
)

# Vitiligo
dis_score_Vitiligo_TS = fread("/QRISdata/Q5059/scDRS_res_step1/TabulaSapiens_pc_ortholog/Vitiligo_Jin_2016_NG_tidy.ma_10kb.full_score.gz")
dis_score_Vitiligo_TS_cd8_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "cd8_positive__alpha_beta_t_cell"),
  dis_score_Vitiligo_TS
) %>% mutate(
  Organ = "Atlas - Human",
  file = "Vitiligo_Jin_2016_NG_tidy.ma_10kb.full_score.gz"
)

# Lymphocyte Count
dis_score_Lymphocyte_count_TS = fread("/QRISdata/Q5059/scDRS_res_step1/TabulaSapiens_pc_ortholog/30120_tidy.ma_10kb.full_score.gz")
dis_score_Lymphocyte_count_TS_cd4_positive__alpha_beta_t_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "cd4_positive__alpha_beta_t_cell"),
  dis_score_Lymphocyte_count_TS
) %>% mutate(
  Organ = "Atlas - Human",
  file = "30120_tidy.ma_10kb.full_score.gz"
)

# Atrial Fibrillation (AF)
dis_score_AF_TS = fread("/QRISdata/Q5059/scDRS_res_step1/TabulaSapiens_pc_ortholog/nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl_tidy.ma_10kb.full_score.gz")
dis_score_AF_TS_cardiac_muscle_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "cardiac_muscle_cell"),
  dis_score_AF_TS
) %>% mutate(
  Organ = "Atlas - Human",
  file = "nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl_tidy.ma_10kb.full_score.gz"
)

# Melanoma
dis_score_Melanoma_TS = fread("/QRISdata/Q5059/scDRS_res_step1/TabulaSapiens_pc_ortholog/Melanoma_tidy.ma_10kb.full_score.gz")
dis_score_Melanoma_TS_melanocyte = inner_join(
  cell_df %>% filter(cell_ontology_class == "melanocyte"),
  dis_score_Melanoma_TS
) %>% mutate(
  Organ = "Atlas - Human",
  file = "Melanoma_tidy.ma_10kb.full_score.gz"
)

# Heart Failure (HF)
dis_score_HF_TS = fread("/QRISdata/Q5059/scDRS_res_step1/TabulaSapiens_pc_ortholog/HERMES_Jan2019_HeartFailure_summary_data.txt_tidy.ma_10kb.full_score.gz")
dis_score_HF_TS_cardiac_muscle_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "cardiac_muscle_cell"),
  dis_score_HF_TS
) %>% mutate(
  Organ = "Atlas - Human",
  file = "HERMES_Jan2019_HeartFailure_summary_data.txt_tidy.ma_10kb.full_score.gz"
)

# Vitamin D Levels
dis_score_VitD_TS = fread("/QRISdata/Q5059/scDRS_res_step1/TabulaSapiens_pc_ortholog/Revezetal2020_vitD_tidy.ma_10kb.full_score.gz")
dis_score_VitD_TS_keratinocyte = inner_join(
  cell_df %>% filter(cell_ontology_class == "keratinocyte"),
  dis_score_VitD_TS
) %>% mutate(
  Organ = "Atlas - Human",
  file = "Revezetal2020_vitD_tidy.ma_10kb.full_score.gz"
)

# Forced Expiratory Volume in 1 Second (FEV1)
dis_score_FEV1_TS = fread("/QRISdata/Q5059/scDRS_res_step1/TabulaSapiens_pc_ortholog/FEV1_tidy.ma_10kb.full_score.gz")
dis_score_FEV1_TS_bronchial_smooth_muscle_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "bronchial_smooth_muscle_cell"),
  dis_score_FEV1_TS
) %>% mutate(
  Organ = "Atlas - Human",
  file = "FEV1_tidy.ma_10kb.full_score.gz"
)


list_of_dfs = list(
  dis_score_Celiac_disease_TS_cd4_positive__alpha_beta_t_cell,
  dis_score_SBP_TS_smooth_muscle_cell,
  dis_score_Red_blood_cell_count_TS_erythrocyte,
  dis_score_Crohns_disease_TS_cd4_positive__alpha_beta_t_cell,
  dis_score_IBD_TS_cd4_positive__alpha_beta_t_cell,
  dis_score_Vitiligo_TS_cd8_positive__alpha_beta_t_cell,
  dis_score_Lymphocyte_count_TS_cd4_positive__alpha_beta_t_cell,
  dis_score_AF_TS_cardiac_muscle_cell,
  dis_score_Melanoma_TS_melanocyte,
  dis_score_HF_TS_cardiac_muscle_cell,
  dis_score_VitD_TS_keratinocyte,
  dis_score_FEV1_TS_bronchial_smooth_muscle_cell,
  dis_score_plt_count_TS_platelet
)

# Merge all dataframes using bind_rows
merged_df_TS = bind_rows(list_of_dfs)

# Save the merged dataframe to a gzipped CSV file
fwrite(merged_df_TS, "/scratch/project_mnt/S0007/uqali4/NM_dt_results/merged_disease_scores_TS.csv.gz", row.names = FALSE, compress = "gzip")
































