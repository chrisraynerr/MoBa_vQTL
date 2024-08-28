#!/bin/bash
#SBATCH --job-name=Plot_TopSnps
#SBATCH --output=/cluster/projects/p805/crayner/logs/%x.%j.out
#SBATCH --account=p805_tsd
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G

source /cluster/bin/jobsetup
cd /cluster/projects/p805/crayner/vQTL
module load R/4.0.0-foss-2020a

rm /cluster/projects/p805/crayner/vQTL/Data/${1}*TopSnps*

Rscript Plot_TopSnps.R --PheDir ${1}

mkdir /tsd/p805/data/durable/projects/crayner/vQTL/Plots_Snps/${1}
scp /cluster/projects/p805/crayner/vQTL/Plots_Snps/${1}/* \
/tsd/p805/data/durable/projects/crayner/vQTL/Plots_Snps/${1}/

scp /cluster/projects/p805/crayner/vQTL/Tables_Snps/${1}/* \
/tsd/p805/data/durable/projects/crayner/vQTL/Tables_Snps/${1}/

# sh DataTransferCluster.sh $TSD/vQTL/Plot_TopSnps.R Plot_TopSnps.R
# sh DataTransferCluster.sh $TSD/vQTL/Plot_TopSnps.sh Plot_TopSnps.sh
# sbatch Plot_TopSnps.sh ChAnth