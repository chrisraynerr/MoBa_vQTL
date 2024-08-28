#!/bin/bash
#SBATCH --job-name=vQTLRegenie_LocoPredQt
#SBATCH --output=/cluster/projects/p805/crayner/logs/%x.%j.out
#SBATCH --account=p805
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8GB

source /cluster/bin/jobsetup
cd /cluster/projects/p805/crayner/vQTL

REGENIE=/cluster/projects/p805/crayner/software/regenie-master/regenie
CL_GDIR=/cluster/projects/p805/data/genetics/MoBaPsychGen_v1

echo INPUTS
PHE_FILE=${1}
OUT_FILE=${PHE_FILE%.*}
OUT_DIR=RegenieLocoPred
mkdir ${OUT_DIR}
mkdir ${OUT_DIR}/${OUT_FILE}/


${REGENIE} \
  --step 1 \
  --bed ${CL_GDIR}/MoBaPsychGen_v1-ec-eur-batch-basic-qc \
  --bsize 1000 \
  --extract ${CL_GDIR}/regenie_500k_snpsTEST.txt \
  --phenoFile Data/${PHE_FILE} \
  --covarFile /cluster/projects/p805/crayner/data/Cov/20SnpPc3GenoBatchCentrePcN207283.tsv \
  --force-qt \
  --lowmem \
  --lowmem-prefix ${OUT_DIR}/${OUT_FILE}/tmp \
  --out ${OUT_DIR}/${OUT_FILE}/

# USAGE:
# cd /cluster/projects/p805/crayner/vQTL
# SCRIPT=Regenie_LocoPredQt.sh
# sh DataTransferCluster.sh ${TSD}/vQTL/${SCRIPT} ${CLUSTER}/vQTL/${SCRIPT}
# PHE_FILE=ChPsyGenQrsAsr.tsv # ChPsyFsxAsr.tsv
# sh DataTransferCluster.sh ${TSD}/vQTL/Data/${PHE_FILE} ${CLUSTER}/vQTL/Data/${PHE_FILE}
# sbatch Regenie_LocoPredQt.sh ${PHE_FILE}
# sbatch Regenie_LocoPredQt.sh ChPsyGenFscAsr.tsv

