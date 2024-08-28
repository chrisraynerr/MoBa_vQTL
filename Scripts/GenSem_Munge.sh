#!/bin/bash
#SBATCH --job-name=GenSemMunge
#SBATCH --output=/cluster/projects/p805/crayner/logs/%x.%j.out
#SBATCH --account=p805_tsd
#SBATCH --time=47:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16GB

source /cluster/bin/jobsetup
module load R/4.0.0-foss-2020a

cd /cluster/projects/p805/crayner/vQTL
Rscript GenSem_Munge.R --sumstats=${1}

# sh DataTransferCluster.sh $TSD/vQTL/GenSem_Munge.R GenSem_Munge.R
# sh DataTransferCluster.sh $TSD/vQTL/GenSem_Munge.sh GenSem_Munge.sh
# sbatch GenSem_Munge.sh /cluster/projects/p805/crayner/vQTL/GwaSumStats/ChAnth