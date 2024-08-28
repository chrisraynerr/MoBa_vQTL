#!/bin/bash
#SBATCH --job-name=vQTL_Regenie_GwaQt
#SBATCH --output=/cluster/projects/p805/crayner/logs/%x.%j.out
#SBATCH --account=p805_tsd
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-22

source /cluster/bin/jobsetup
cd /cluster/projects/p805/crayner/vQTL

PHE_FILE=${1}
OUT_FILE=${PHE_FILE%.*}
mkdir RegenieGwa/${OUT_FILE}

/cluster/projects/p805/crayner/software/regenie-master/regenie \
  --step 2 \
  --bed /cluster/projects/p805/data/genetics/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc \
  --bsize 200 \
  --covarFile /cluster/projects/p805/crayner/data/Cov/20SnpPc3GenoBatchCentrePcN207283.tsv \
  --phenoFile Data/${PHE_FILE} \
  --chr ${SLURM_ARRAY_TASK_ID} \
  --force-qt \
  --pred RegenieLocoPred/ChPsyGenQrsAsr_pred.list \
  --out RegenieGwa/${OUT_FILE}/Chr${SLURM_ARRAY_TASK_ID}

  # --pred RegenieLocoPred/${OUT_FILE}_pred.list \
  # --chrList 1,2,3,4,5,6,7,8,9,10,11 \

# USAGE:
# SCRIPT=Regenie_GwaQtArray.sh
# PHE_FILE=ChPsyGenQrsAsr.tsv
# sh DataTransferCluster.sh ${TSD}/vQTL/${SCRIPT} ${CLUSTER}/vQTL/${SCRIPT}
# sbatch Regenie_GwaQtArray.sh ${PHE_FILE}
# 
# sbatch Regenie_GwaQtArray.sh ChPsyQrxAsr.tsv
# sbatch Regenie_GwaQtArray.sh ChPsyFsxAsr.tsv

