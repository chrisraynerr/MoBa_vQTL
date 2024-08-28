#!/bin/sh
#SBATCH --output=/cluster/projects/p805/crayner/logs/%x.%j.out
#SBATCH --account=p805_tsd
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2GB

#################################################
# USAGE:
# cd /cluster/projects/p805/crayner/vQTL
# FILE=ChAnth
# ANALYSIS=RegenieMQt
# sbatch Regenie_ReformatGwaSumstats.sh ${FILE} ${ANALYSIS}

#################################################

source /cluster/bin/jobsetup

cd /cluster/projects/p805/crayner/vQTL
FILE=${1}
FILE=${FILE%.*}
ANALYSIS=${2}
# FILE=ChPsyQrxAsr
# FILE=ChPsyFsxAsr
# ANALYSIS=RegenieMQt

cd RegenieGwa/${FILE}/

ls Chr1_* > LIST
sed -i 's/Chr1_//g' LIST
sed -i 's/.regenie//g' LIST

cd /cluster/projects/p805/crayner/vQTL

mkdir GwaSumStats/${FILE}/

while read LIST 
do 
  if [[ ! -f  GwaSumStats/${FILE}/${LIST}_${ANALYSIS}.txt.gz ]]; then
  echo -e 'chr pos rsid a0 a1 af n test beta se chi2 log10p extra p' > GwaSumStats/${FILE}/${LIST}.head
  echo "combine all chromosomes into one file..."
  for i in {1..22} 
  do
    awk 'FNR>1 {print $0, 10^(-1*$12)}' RegenieGwa/${FILE}/Chr${i}_${LIST}.regenie > RegenieGwa/${FILE}/Chr${i}_${LIST}.tmp
  echo ${i}
  done
  cat GwaSumStats/${FILE}/${LIST}.head RegenieGwa/${FILE}/Chr*_${LIST}.tmp >> GwaSumStats/${FILE}/${LIST}_${ANALYSIS}.txt
  Nsnps=$(wc -l GwaSumStats/${FILE}/${LIST}_${ANALYSIS}.txt)
  echo ${Nsnps}
  echo "gzip processed sumstats file..."
  gzip -f GwaSumStats/${FILE}/${LIST}_${ANALYSIS}.txt
  rm RegenieGwa/${FILE}/Chr*_${LIST}.tmp
  rm GwaSumStats/${FILE}/${LIST}.head
  fi
done < RegenieGwa/${FILE}/LIST

while read LIST 
do 
  if [[ -f  GwaSumStats/${FILE}/${LIST}_${ANALYSIS}.txt.gz ]]; then
  rm RegenieGwa/${FILE}/Chr*_${LIST}.regenie
  fi
done < RegenieGwa/${FILE}/LIST


# sh DataTransferCluster.sh $TSD/vQTL/Regenie_ReformatGwaSumstats.sh ./Regenie_ReformatGwaSumstats.sh
# sbatch ./Regenie_ReformatGwaSumstats.sh ChPsyGenQrsAsr RegenieVQt
# sbatch ./Regenie_ReformatGwaSumstats.sh ChPsyGenFscAsr RegenieMQt

# for i in $(ls RegenieGwa/Chr1_); do sbatch Regenie_ReformatGwaSumstats.sh $i; done;
