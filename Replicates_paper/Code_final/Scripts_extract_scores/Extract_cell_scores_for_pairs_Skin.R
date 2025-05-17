library(data.table)
library(dplyr)

# Read the cell metadata dataframe for Skin
cell_df = fread("/scratch/project_mnt/S0007/uqali4/Cheng_2018_Cell_Reports_pcc.cell_id_cell_type.csv")

# Process each dataset individually

# Vitamin D Levels
dis_score_VitD_Skin = fread("/QRISdata/Q3999/tmp/scDRS_res_step1/Cheng_2018_Cell_Reports_pc/mBAT/Revezetal2020_vitD_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = index)

dis_score_VitD_Skin_keratinocyte = inner_join(
  cell_df %>% filter(cell_ontology_class == "keratinocyte"),
  dis_score_VitD_Skin,
  by = "cell_id"
) %>% mutate(
  Organ = "Skin",
  file = "Revezetal2020_vitD_tidy.ma_10kb.full_score.gz"
)

# Melanoma
dis_score_Melanoma_Skin = fread("/QRISdata/Q3999/tmp/scDRS_res_step1/Cheng_2018_Cell_Reports_pc/mBAT/Melanoma_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = index)

dis_score_Melanoma_Skin_melanocyte = inner_join(
  cell_df %>% filter(cell_ontology_class == "melanocyte"),
  dis_score_Melanoma_Skin,
  by = "cell_id"
) %>% mutate(
  Organ = "Skin",
  file = "Melanoma_tidy.ma_10kb.full_score.gz"
)

# Now, merge all the individual dataframes into a single dataframe
list_of_dfs = list(
  dis_score_VitD_Skin_keratinocyte,
  dis_score_Melanoma_Skin_melanocyte
)

# Merge all dataframes using bind_rows
merged_df_Skin = bind_rows(list_of_dfs)

# Save the merged dataframe to a gzipped CSV file
fwrite(merged_df_Skin, "/scratch/project_mnt/S0007/uqali4/NM_dt_results/merged_disease_scores_Skin.csv.gz", row.names = FALSE, compress = "gzip")
