#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --job-name=Kamath_2022_normal_expr_gene_magma_geneset_benchmark_100kb
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
trait_magma_dir="/scratch/project/genetic_data_analysis/uqali4/MAGMA_10kb_UKB20k_56778genes"
metrics_dirs=(TDEP cepo_norm det_esw_s ep_esw_s esmu ges_esw_s nsi_esw_s sclinker) 
for dir_name in "${metrics_dirs[@]}"; do
  COVAR_FILE="/scratch/project/genetic_data_analysis/uqali4/MG_set_conAnnot/Kamath_2022_normal_expr_gene.${dir_name}.csv"
  mkdir -p /scratch/user/uqali4/geneset_res/Con_MG_set/Kamath_2022_normal_expr_gene/${dir_name}
  new_res_dir="/scratch/user/uqali4/geneset_res/Con_MG_set/Kamath_2022_normal_expr_gene/${dir_name}"
  /scratch/project/genetic_data_analysis/uqali4/software/magma --gene-results $trait_magma_dir/magma.${trait}.genes.raw --gene-covar $COVAR_FILE --out $new_res_dir/${trait}
done

