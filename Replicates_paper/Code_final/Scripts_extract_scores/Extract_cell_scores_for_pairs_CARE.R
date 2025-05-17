library(data.table)
library(dplyr)

# Read the cell metadata dataframe for CARE
cell_df = fread("/scratch/project_mnt/S0007/uqali4/CARE_snRNA_Heart_pc.cell_id_celltype.csv")

# Process each dataset individually

# Atrial Fibrillation (AF)
dis_score_AF_CARE = fread("/QRISdata/Q5059/scDRS_res_step1/Example_atlas_vs_singleorgan/CARE_pconly/nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl_tidy.ma_10kb.full_score.gz") %>%
  dplyr::rename(cell_id = barcode) %>%
  select(-matches("^ctrl_norm_score_[0-9]+$"))

# Check updated column

dis_score_AF_CARE_cardiac_muscle_cell = inner_join(
  cell_df %>% filter(cell_ontology_class == "Atrial_cardiomyocyte"),
  dis_score_AF_CARE,
  by = "cell_id"
) %>% mutate(
  Organ = "Heart",
  file = "nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl_tidy.ma_10kb.full_score.gz"
)


# Now, merge all the individual dataframes into a single dataframe
list_of_dfs = list(
  dis_score_AF_CARE_cardiac_muscle_cell
)

# Merge all dataframes using bind_rows
merged_df_CARE = bind_rows(list_of_dfs)

# Save the merged dataframe to a gzipped CSV file
fwrite(merged_df_CARE, "/scratch/project_mnt/S0007/uqali4/NM_dt_results_new_score/merged_disease_scores_CARE.csv.gz", row.names = FALSE, compress = "gzip")
