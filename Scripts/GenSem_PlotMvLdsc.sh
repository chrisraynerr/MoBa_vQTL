#!/bin/bash
#SBATCH --job-name=GenSemPlot
#SBATCH --output=/cluster/projects/p805/crayner/logs/%x.%j.out
#SBATCH --account=p805_tsd
#SBATCH --time=47:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16GB

source /cluster/bin/jobsetup
module load R/4.0.0-foss-2020a

cd /cluster/projects/p805/crayner/vQTL/GenSem

SsSet=ChPsyFsxAsr
SsSet=ChPsyQrxAsr

for SsSet in ChPsyFsxAsr ChPsyQrxAsr
  do
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  grep '^Heritability Results for trait'  mvLDSC_${SsSet}.log > mvLDSC_${SsSet}_Traits.txt
  sed -i 's/Heritability Results for trait: \/cluster\/projects\/p805\/crayner\/vQTL\/GwaSumStats\///g' mvLDSC_${SsSet}_Traits.txt
  sed -i 's/_Munged.sumstats.gz//g' mvLDSC_${SsSet}_Traits.txt
  sed -i 's/.Munged.sumstats.gz//g' mvLDSC_${SsSet}_Traits.txt
  sed -i 's/^.*\///g' mvLDSC_${SsSet}_Traits.txt
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  grep '^Intercept'  mvLDSC_${SsSet}.log > mvLDSC_${SsSet}_h2Int.txt
  sed -i 's/Intercept //g' mvLDSC_${SsSet}_h2Int.txt
  sed -i 's/Intercept: //g' mvLDSC_${SsSet}_h2Int.txt
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  grep '^Total Observed Scale h2'  mvLDSC_${SsSet}.log > mvLDSC_${SsSet}_h2.txt
  sed -i 's/Total Observed Scale h2: //g' mvLDSC_${SsSet}_h2.txt
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  grep '^Ratio:'  mvLDSC_${SsSet}.log > mvLDSC_${SsSet}_h2R.txt
  sed -i 's/Ratio: //g' mvLDSC_${SsSet}_h2R.txt
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  grep '^h2 Z:'  mvLDSC_${SsSet}.log > mvLDSC_${SsSet}_h2Z.txt
  sed -i 's/h2 Z: //g' mvLDSC_${SsSet}_h2Z.txt
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  paste mvLDSC_${SsSet}_Traits.txt mvLDSC_${SsSet}_h2.txt mvLDSC_${SsSet}_h2Int.txt mvLDSC_${SsSet}_h2R.txt mvLDSC_${SsSet}_h2Z.txt > mvLDSC_${SsSet}_h2FullTable.txt
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sed -i 's/[)(]//g'  mvLDSC_${SsSet}_h2FullTable.txt
  sed -i 's/ChFPsy/F_/g' mvLDSC_${SsSet}_h2FullTable.txt
  sed -i 's/ChMPsy/M_/g' mvLDSC_${SsSet}_h2FullTable.txt
  sed -i 's/_Ch/_/g' mvLDSC_${SsSet}_h2FullTable.txt
  sed -i 's/_Ch/_/g' mvLDSC_${SsSet}_h2FullTable.txt
  sed -i 's/TscQ/_Q/g' mvLDSC_${SsSet}_h2FullTable.txt
  sed -i 's/_RegenieMQt//g' mvLDSC_${SsSet}_h2FullTable.txt
  sed -i 's/[)(]//g' mvLDSC_${SsSet}_h2FullTable.txt
  sed -i 's/ /\t/g' mvLDSC_${SsSet}_h2FullTable.txt
  sed -i 's/  */ /g' mvLDSC_${SsSet}_h2FullTable.txt
  sed -i.bak $'s/\t/ /g'  mvLDSC_${SsSet}_h2FullTable.txt
  sed -i 's/\/cluster\/projects\/p805\/crayner\/vQTL\/GwaSumStats\// /g' mvLDSC_${SsSet}_h2FullTable.txt
  sed -i 's/_Munged.sumstats.gz//g' mvLDSC_${SsSet}_h2FullTable.txt
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  scp mvLDSC_${SsSet}_h2FullTable.txt NonNegMvLDSC_${SsSet}_h2FullTable.txt
  sed -i '/-/d' NonNegMvLDSC_${SsSet}_h2FullTable.txt
done

# for SsSet in ChPsyFsxAsr
for SsSet in ChPsyGenFscAsr ChPsyGenQrsAsr
  do
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  grep '^Genetic Correlation between' mvLDSC_${SsSet}.log > mvLDSC_${SsSet}_rg.txt
  sed -i 's/Genetic Correlation between //g' mvLDSC_${SsSet}_rg.txt
  sed -i 's/and //g' mvLDSC_${SsSet}_rg.txt
  sed -i 's/\://g' mvLDSC_${SsSet}_rg.txt
  sed -i 's/_RegenieMQt//g' mvLDSC_${SsSet}_rg.txt
  sed -i 's/.Munged.sumstats.gz//g' mvLDSC_${SsSet}_rg.txt
  sed -i 's/[)(]//g' mvLDSC_${SsSet}_rg.txt
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  grep '^Results for genetic covariance between:' mvLDSC_${SsSet}.log > mvLDSC_${SsSet}_pval1.txt
  sed -i 's/Results for genetic covariance between: //g' mvLDSC_${SsSet}_pval1.txt
  sed -i 's/\/cluster\/projects\/p805\/crayner\/vQTL\/GwaSumStats\// /g' mvLDSC_${SsSet}_pval1.txt
  sed -i 's/_RegenieMQt.Munged.sumstats.gz//g' mvLDSC_${SsSet}_pval1.txt
  sed -i 's/and //g' mvLDSC_${SsSet}_pval1.txt
  sed -i 's/\// /g' mvLDSC_${SsSet}_pval1.txt
  sed -i 's/^ //g' mvLDSC_${SsSet}_pval1.txt
  sed -i 's/\t/ /g' mvLDSC_${SsSet}_pval1.txt
  awk '{print $2,$4}' mvLDSC_${SsSet}_pval1.txt > mvLDSC_${SsSet}_pval2.txt
  grep '^g_cov P-value: ' mvLDSC_${SsSet}.log > mvLDSC_${SsSet}_pval3.txt
  sed -i 's/g_cov P-value: //g' mvLDSC_${SsSet}_pval3.txt
  paste mvLDSC_${SsSet}_pval2.txt mvLDSC_${SsSet}_pval3.txt > MvLDSC_${SsSet}_rG_pval.txt
  sed -i 's/\t/ /g' MvLDSC_${SsSet}_rG_pval.txt
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
done 

rm *pval1*
rm *pval2*
rm *pval3*

# for SsSet in ChPsyQrxAsr
for SsSet in ChPsyGenQrsAsr
  do
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  grep '^Genetic Correlation between' NonNegMvLDSC_${SsSet}.log > NonNegMvLDSC_${SsSet}_rg.txt
  sed -i 's/Genetic Correlation between //g' NonNegMvLDSC_${SsSet}_rg.txt
  sed -i 's/and //g' NonNegMvLDSC_${SsSet}_rg.txt
  sed -i 's/\://g' NonNegMvLDSC_${SsSet}_rg.txt
  sed -i 's/_RegenieMQt//g' NonNegMvLDSC_${SsSet}_rg.txt
  sed -i 's/.Munged.sumstats.gz//g' NonNegMvLDSC_${SsSet}_rg.txt
  sed -i 's/[)(]//g' NonNegMvLDSC_${SsSet}_rg.txt
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  grep '^Results for genetic covariance between:' NonNegMvLDSC_${SsSet}.log > NonNegMvLDSC_${SsSet}_pval1.txt
  sed -i 's/Results for genetic covariance between: //g' NonNegMvLDSC_${SsSet}_pval1.txt
  sed -i 's/\/cluster\/projects\/p805\/crayner\/vQTL\/GwaSumStats\// /g' NonNegMvLDSC_${SsSet}_pval1.txt
  sed -i 's/_Munged.sumstats.gz//g' NonNegMvLDSC_${SsSet}_pval1.txt
  sed -i 's/and //g' NonNegMvLDSC_${SsSet}_pval1.txt
  sed -i 's/\// /g' NonNegMvLDSC_${SsSet}_pval1.txt
  sed -i 's/^ //g' NonNegMvLDSC_${SsSet}_pval1.txt
  sed -i.bak $'s/\t/ /g' NonNegMvLDSC_${SsSet}_pval1.txt
  awk '{print $2,$4}' NonNegMvLDSC_${SsSet}_pval1.txt > NonNegMvLDSC_${SsSet}_pval2.txt
  grep '^g_cov P-value: ' NonNegMvLDSC_${SsSet}.log > NonNegMvLDSC_${SsSet}_pval3.txt
  sed -i 's/g_cov P-value: //g' NonNegMvLDSC_${SsSet}_pval3.txt
  paste NonNegMvLDSC_${SsSet}_pval2.txt NonNegMvLDSC_${SsSet}_pval3.txt > NonNegMvLDSC_${SsSet}_pval.txt
  sed -i.bak $'s/\t/ /g' NonNegMvLDSC_${SsSet}_pval.txt
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
done
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd GwaSumStats/${SsSet}/
ls *.Munged.sum* > TRAIT_LIST
while read TRAIT_LIST; do zcat $TRAIT_LIST | awk '{ total += $2; count++ } END { print total/count }' >> MEAN_N; done < TRAIT_LIST
paste TRAIT_LIST MEAN_N > SampleSizes.txt
sed -i 's/_RegenieMQt.Munged.sumstats.gz//g' SampleSizes.txt
sed -i 's/Q/ Q/g' SampleSizes.txt

scp mvLDSC_${SsSet}_h2FullTable.txt $TSD/vQTL/Tables/
scp SampleSizes.txt $TSD/vQTL/Tables/

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rscript ../GenSem_PlotMvLdsc.R ${SsSet}

done

# sh DataTransferCluster.sh $TSD/vQTL/GenSem_PlotMvLdsc.R $CLUSTER/vQTL/GenSem_PlotMvLdsc.R
# sh DataTransferCluster.sh $TSD/vQTL/GenSem_PlotMvLdsc.sh $CLUSTER/vQTL/GenSem_PlotMvLdsc.sh


