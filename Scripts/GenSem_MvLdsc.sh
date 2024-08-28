#!/bin/bash
#SBATCH --job-name=GenSemMvLdsc
#SBATCH --output=/cluster/projects/p805/crayner/logs/%x.%j.out
#SBATCH --account=p805_tsd
#SBATCH --time=47:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16GB

source /cluster/bin/jobsetup
module load R/4.0.0-foss-2020a

cd /cluster/projects/p805/crayner/vQTL
Rscript GenSem_MvLdsc.R

# cd /cluster/projects/p805/crayner/vQTL
# sh DataTransferCluster.sh $TSD/vQTL/GenSem_MvLdsc.R GenSem_MvLdsc.R
# sh DataTransferCluster.sh $TSD/vQTL/GenSem_MvLdsc.sh GenSem_MvLdsc.sh
# sbatch GenSem_MvLdsc.sh