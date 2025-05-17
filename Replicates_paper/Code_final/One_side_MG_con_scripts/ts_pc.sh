#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --job-name=ts_pc_mg
#SBATCH --error=slurm_%A_%a.error
#SBATCH --output=slurm_%A_%a.out
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH --array=1-39

i=$SLURM_ARRAY_TASK_ID

TASKFILE="/QRISdata/Q5729/job_list.txt"
trait=$(sed -n "${i}p" $TASKFILE)
trait=$(echo "$trait" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

echo "Trimmed Trait: '$trait'"
echo ${i}
trait_magma_dir="/scratch/project/genetic_data_analysis/uqali4/MAGMA_10kb_UKB20k_56778genes_pc_tms_ts_consist"
metrics_dirs=("TDEP" "cepo_norm" "det.esw_s" "ep.esw_s" "esmu" "ges.esw_s" "nsi.esw_s" "sclinker") 
for dir_name in "${metrics_dirs[@]}"; do
  COVAR_FILE="/scratch/project/genetic_data_analysis/uqali4/MG_set_conAnnot/ts_pc.${dir_name}.csv"
  mkdir -p /scratch/user/uqali4/geneset_res/Con_MG_set_oneside/ts_pc/${dir_name}
  new_res_dir="/scratch/user/uqali4/geneset_res/Con_MG_set_oneside/ts_pc/${dir_name}"
  /scratch/project/genetic_data_analysis/uqali4/software/magma --gene-results $trait_magma_dir/magma.${trait}.genes.raw --gene-covar $COVAR_FILE --model direction=pos --out $new_res_dir/${trait}
done
