#!/bin/bash
#SBATCH --job-name=Plot_GwaQQ
#SBATCH --output=/cluster/projects/p805/crayner/logs/%x.%j.out
#SBATCH --account=p805_tsd
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G

source /cluster/bin/jobsetup
cd /cluster/projects/p805/crayner/vQTL
module load R/4.0.0-foss-2020a

Rscript Plot_GwaQQ.R --Phe=${1} --PheDir=${2}

# sh DataTransferCluster.sh $TSD/vQTL/Plot_GwaQQ.R Plot_GwaQQ.R
# sh DataTransferCluster.sh $TSD/vQTL/Plot_GwaQQ.sh Plot_GwaQQ.sh
# 
# # ls GwaSumStats/ChAnthQrs/*.txt.gz > PlotTraitList
# # sed -i 's/GwaSumStats\/ChAnthQrs\///g' PlotTraitList
# # sed -i 's/_RegenieVQt.txt.gz//g' PlotTraitList
# 
# while read PlotTraitList
# do
# sbatch Plot_GwaQQ.sh $PlotTraitList ChAnth
# done < PlotTraitList