#!/usr/bin/Rscript
cat("\n### Preparing workspace \n")

remove(list = ls());start.time <- Sys.time()
local({
  r = getOption("repos")
  r["CRAN"] <- "https://cran.tsd.usit.no"
  r["BIOCONDUCTOR"] <- "https://bioconductor.tsd.usit.no"
  r["BIOCONDUCTOR-annotation"] <- "https://bioconductor.tsd.usit.no/data/annotation"
  r["BIOCONDUCTOR-experiment"] <- "https://bioconductor.tsd.usit.no/data/experiment"
  r["BioC_mirror"] <- "https://bioconductor.tsd.usit.no"
  options(repos=r)
  options(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = FALSE)
  options(download.file.method = "libcurl")
})

.libPaths(c("/cluster/projects/p805/crayner/software/R/4.0",
           "/cluster/projects/p805/software/R/4.0"))

cat("\n... loading packages \n")
packages <- c("optparse","R.utils","GenomicSEM", "bigsnpr", "data.table", "parallel", "doParallel",
              "dplyr", "stringr")

installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages==F)){install.packages(packages[!installed_packages],dependencies = T)}
suppressPackageStartupMessages(invisible(lapply(packages,library,character.only=T)))


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

option_list = list(
  make_option(
    c("-s", "--sumstats"), 
    type="character", 
    default=NA,
    help="GWAS summary statistics directory path",
    metavar="character")
)

opt_parser = OptionParser(option_list=option_list); opt = parse_args(opt_parser);
print_help(opt_parser); print(opt)

SsSet <- "ChAnthQrs"
SsSet <- opt$sumstats
SsDir <- paste0("/cluster/projects/p805/crayner/vQTL/GwaSumStats/",SsSet)

# SsDir <- "/cluster/projects/p805/crayner/vQTL/GwaSumStats/ChPsyGenQrsAsr"
# SsDir <- "/cluster/projects/p805/crayner/vQTL/GwaSumStats/ChPsyGenFscAsr"
# SsDir <- "/cluster/projects/p805/crayner/vQTL/GwaSumStats/ChAnth"

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

T0 <- Sys.time()
cat(paste0("\n*** Beginning analysis at ",T0," ***\n"))

NCORES    <- nb_cores()
cat(paste0("\n... using ",NCORES," cores \n"))

cat("\n ... registering cluster \n")
Operating <- Sys.info()[['sysname']]
if (Operating != "Windows") {
  cl <- makeCluster(NCORES, type = "FORK", outfile = "")
} else {
  cl <- makeCluster(NCORES, type = "PSOCK", outfile = "")
}

registerDoParallel(cl)
on.exit(parallel::stopCluster(cl), add = TRUE)

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to create all pairwise comparisons between all variables in a character vector
pairwise_comparisons <- function(vars) {
  n <- length(vars)
  var1 <- c();  var2 <- c()
  for (i in 1:n) {
    for (j in i:n) {
      var1 <- c(var1, vars[i])
      var2 <- c(var2, vars[j])
    }}
  return(list(var1, var2))
}

compute_SE_of_rg <- function(GsObj) {
  S <- GsObj$S
  V <- GsObj$V
  ratio <- tcrossprod(1 / sqrt(diag(S)))
  S_Stand <- S * ratio
  scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)
  #rescale the sampling correlation matrix by the appropriate diagonals
  V_Stand <- V * tcrossprod(scaleO)
  #enter SEs from diagonal of standardized V
  r <- nrow(S)
  SE_Stand <- matrix(0, r, r)
  SE_Stand[lower.tri(SE_Stand, diag = TRUE)] <- sqrt(diag(V_Stand))
  as.data.frame(as.matrix(Matrix::forceSymmetric(SE_Stand, "L")))
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Multivariable LDSC
# for(SsSet in c("ChPsyGenFscAsr")){ #,"ChPsyFsxAsr","ChFPsyQrs","ChFPsyTsc","ChMPsyQrs","ChMPsyTsc")){
cat(paste0(SsSet, "\n"))
# SsDir <- paste0("/cluster/projects/p805/crayner/vQTL/GwaSumStats/",SsSet)
  ## vector of munged summary statisitcs
  files <- list.files(SsDir,full.names=T,pattern="Munged.sumstats.gz")
  ## name the traits 
  names <- gsub("_Munged.sumstats.gz","",basename(files))
  #define the reference file being used to allign alleles across summary stats
  hm3 <- "/cluster/projects/p805/crayner/data/RefPanels/eur_w_ld_chr/w_hm3.snplist"
  #the folder of LD scores
  ld  <- "/cluster/projects/p805/crayner/data/RefPanels/eur_w_ld_chr"
  #the folder of LD weights [typically the same as folder of LD scores]
  wld <- "/cluster/projects/p805/crayner/data/RefPanels/eur_w_ld_chr"
  hd  <- "/cluster/projects/p805/crayner/data/RefPanels/UKB_imputed_SVD_eigen99_extraction"
  ## enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
  sample.prev     <- c(rep(NA,length(names)))
  ## vector of population prevalences
  population.prev <- c(rep(NA,length(names)))
  
  cat("\n*** Running MV LDSC ***\n")
  
  if(!file.exists(paste0("/cluster/projects/p805/crayner/vQTL/GenSem/mvLDSC_",SsSet,".rds"))){
  
  LDSCoutput <-
    ldsc(
      traits          = files,
      trait.names     = names,
      sample.prev     = sample.prev,
      population.prev = population.prev,
      ld              = ld,
      wld             = wld,
      ldsc.log        = paste0("/cluster/projects/p805/crayner/vQTL/GenSem/mvLDSC_",SsSet)
    )
  
  saveRDS(LDSCoutput,paste0("/cluster/projects/p805/crayner/vQTL/GenSem/mvLDSC_",SsSet,".rds"))
  
  } else {
      LDSCoutput<-readRDS(paste0("/cluster/projects/p805/crayner/vQTL/GenSem/mvLDSC_",SsSet,".rds"))
    }

    
  LOG <- paste0("/cluster/projects/p805/crayner/vQTL/GenSem/mvLDSC_",SsSet,"_ldsc.log")
  h2  <- paste0("/cluster/projects/p805/crayner/vQTL/GenSem/mvLDSC_",SsSet)
  
  if(!file.exists(paste0(h2,"_h2NonNeg.txt"))){
  
  system(paste0("grep '^Heritability Results for trait' ",LOG, " > ",h2,"_Traits.txt"))
  system(paste0("sed -i 's/Heritability Results for trait: //g' ",h2,"_Traits.txt"))
  system(paste0("grep '^Total Observed Scale h2' ",LOG, " > ",h2,"_h2.txt"))
  system(paste0("sed -i 's/Total Observed Scale h2: //g' ",h2,"_h2.txt"))
  system(paste0("grep '^Intercept' ",LOG, " > ",h2,"_h2I.txt"))
  system(paste0("sed -i 's/Intercept: //g' ",h2,"_h2I.txt"))  
  system(paste0("sed -i 's/Intercept //g' ",h2,"_h2I.txt"))  
  system(paste0("grep '^Lambda' ",LOG, " > ",h2,"_h2L.txt"))
  system(paste0("sed -i 's/Lambda GC: //g' ",h2,"_h2L.txt"))
  system(paste0("grep '^Mean Chi' ",LOG, " > ",h2,"_h2C.txt"))
  system(paste0("sed -i 's/Mean Chi^2 across remaining SNPs: //g' ",h2,"_h2C.txt"))
  system(paste0("grep '^Ratio:' ",LOG, " > ",h2,"_h2R.txt"))
  system(paste0("sed -i 's/Ratio: //g' ",h2,"_h2R.txt"))
  system(paste0("grep '^h2 Z:' ",LOG, " > ",h2,"_h2Z.txt"))
  system(paste0("sed -i 's/h2 Z: //g'  ",h2,"_h2Z.txt"))
  system(paste0("paste ",h2,"_Traits.txt ",h2,"_h2.txt ",h2,"_h2I.txt ",h2,"_h2C.txt ",h2,"_h2L.txt ",h2,"_h2R.txt ",h2,"_h2Z.txt > ",h2,"_h2FullTable.txt"))
  system(paste0("sed -i 's/[)(]//g' ",h2,"_h2FullTable.txt"))
  system(paste0("sed -i 's/ /\t/g' ",h2,"_h2FullTable.txt"))
  system(paste0("sed -i 's/  */ /g'  ",h2,"_h2FullTable.txt"))
  system(paste0("sed -i 's/\t/ /g'  ",h2,"_h2FullTable.txt"))
  system(paste0("sed -i 's/\t/ /g'  ",h2,"_h2FullTable.txt"))
  system(paste0("sed -i 's/\\/cluster\\/projects\\/p805\\/crayner\\/vQTL\\/GwaSumStats\\/",SsSet,"\\///g'  ",h2,"_h2FullTable.txt"))
  system(paste0("sed -i 's/_RegenieVQt.Munged.sumstats.gz//g'  ",h2,"_h2FullTable.txt"))
  system(paste0("sed -i 's/_RegenieMQt.Munged.sumstats.gz//g'  ",h2,"_h2FullTable.txt"))
  
  FullTab        <- fread(paste0(h2,"_h2FullTable.txt"),header=F,sep=" ")
  # FullTab        <- fread(paste0(h2,"_h2FullTable.txt"),header=T,sep=",")
  # FullTab        <- fread(paste0(h2,"_h2FullTable.txt"),quote=F)
  
  names(FullTab) <- c("Trait","h2","SE","Int","IntSE","MeanChi2","Lambda","Ratio","RatioSE", "Z")

  fwrite(FullTab,paste0(h2,"_h2FullTable.txt"),sep="\t",quote=F)
  
  NonNeg <- FullTab %>% filter(h2>0 & Z>2)
  
  fwrite(NonNeg,paste0(h2,"_h2NonNeg.txt"),sep="\t",quote=F)
  }
  
  FullTab <- fread(paste0(h2,"_h2FullTable.txt"),header=T,sep="\t")
  
  NonNeg  <- fread(paste0(h2,"_h2NonNeg.txt"),header=T,sep="\t")
  
  varVec <- colnames(LDSCoutput$S)
  varVec <- str_remove_all(varVec,"_RegenieVQt.Munged.sumstats.gz")
  varVec <- str_remove_all(varVec,"_RegenieMQt.Munged.sumstats.gz")
  
  varMat <- pairwise_comparisons(varVec)
  keepVec<- which(FullTab$Trait %in% NonNeg$Trait)
  varKeep<- NonNeg$Trait

  keepMat <- NULL
  for(i in seq_along(varMat[[1]])){
    keepMat[i] <- 
      ifelse(varMat[[1]][i] %in% varKeep & varMat[[2]][i] %in% varKeep,T,F)
  }
  
  Y <- list()
  Y[["V"]] <- as.matrix(LDSCoutput[["V"]][keepMat,keepMat])
  Y[["S"]] <- as.matrix(LDSCoutput[["S"]][keepVec,keepVec])
  Y[["I"]] <- as.matrix(LDSCoutput[["I"]][keepVec,keepVec])
  Y[["N"]] <- as.matrix(LDSCoutput[["N"]][1,keepMat])
  Y[["m"]] <- LDSCoutput[["m"]]

  saveRDS(Y,paste0("/cluster/projects/p805/crayner/vQTL/GenSem/NonNegMvLDSC_",SsSet,".rds"))
  
  rGmat <- cov2cor(Y$S)
  SEmat <- compute_SE_of_rg(Y); colnames(SEmat) <- colnames(rGmat)
  fwrite(rGmat,paste0("/cluster/projects/p805/crayner/vQTL/GenSem/MvLDSC_",SsSet,"_rG_est.csv"))
  fwrite(SEmat,paste0("/cluster/projects/p805/crayner/vQTL/GenSem/MvLDSC_",SsSet,"_rG_err.csv"))

  system(paste0("grep '^Results for genetic covariance between:' ",LOG," > mvLDSC_",SsSet,"_pval1.txt"))
  system(paste0("sed -i 's/Results for genetic covariance between: //g' mvLDSC_",SsSet,"_pval1.txt"))
  system(paste0("sed -i 's/\/cluster\/projects\/p805\/crayner\/vQTL\/GwaSumStats\// /g' mvLDSC_",SsSet,"_pval1.txt"))
  system(paste0("sed -i 's/_RegenieMQt.Munged.sumstats.gz//g' mvLDSC_",SsSet,"_pval1.txt"))
  system(paste0("sed -i 's/_RegenieVQt.Munged.sumstats.gz//g' mvLDSC_",SsSet,"_pval1.txt"))
  system(paste0("sed -i 's/and //g' mvLDSC_",SsSet,"_pval1.txt"))
  system(paste0("sed -i 's/\// /g' mvLDSC_",SsSet,"_pval1.txt"))
  system(paste0("sed -i 's/^ //g' mvLDSC_",SsSet,"_pval1.txt"))
  system(paste0("sed -i 's/\t/ /g' mvLDSC_",SsSet,"_pval1.txt"))
  system(paste0("awk '{print $2,$4}' mvLDSC_",SsSet,"_pval1.txt > mvLDSC_",SsSet,"_pval2.txt"))
  system(paste0("grep '^g_cov P-value: ' ",LOG,"> mvLDSC_",SsSet,"_pval3.txt"))
  system(paste0("sed -i 's/g_cov P-value: //g' mvLDSC_",SsSet,"_pval3.txt"))
  system(paste0("paste mvLDSC_",SsSet,"_pval2.txt mvLDSC_",SsSet,"_pval3.txt > MvLDSC_",SsSet,"_rG_pval.txt"))
  system(paste0("sed -i 's/\t/ /g' MvLDSC_",SsSet,"_rG_pval.txt"))
  
  system(paste0("ls GwaSumStats/",SsSet,"/*.Munged.sumstats.gz > /cluster/projects/p805/crayner/vQTL/",SsSet,"SampleSizes_traits.txt"))
  system(paste0("while read SST; do zcat $SST | tail -n1 | awk 's+=$2{print s/NR}' >> /cluster/projects/p805/crayner/vQTL/",SsSet,"SampleSizes_N.txt; done < /cluster/projects/p805/crayner/vQTL/",SsSet,"SampleSizes_traits.txt"))
  system(paste0("paste /cluster/projects/p805/crayner/vQTL/",SsSet,"SampleSizes_traits.txt /cluster/projects/p805/crayner/vQTL/",SsSet,"SampleSizes_N.txt > /cluster/projects/p805/crayner/vQTL/",SsSet,"_SampleSizes.txt"))

  system(paste0("sed -i 's/\/cluster\/projects\/p805\/crayner\/vQTL\/GwaSumStats\// /g' mvLDSC_",SsSet,"_pval1.txt"))
  system(paste0("sed -i 's/_RegenieMQt.Munged.sumstats.gz//g' mvLDSC_",SsSet,"_pval1.txt"))
  system(paste0("sed -i 's/_RegenieVQt.Munged.sumstats.gz//g' mvLDSC_",SsSet,"_pval1.txt"))
  
  
  
  C <- commonfactor(Y,estimation="DWLS")
  saveRDS(C,paste0("/cluster/projects/p805/crayner/vQTL/GenSem/MvLDSC_",SsSet,"_CommonfactorModels.rds"))
  
  require(Matrix);require(stats)
  Ssmooth  <- as.matrix((nearPD(Y$S, corr = FALSE))$mat)
  
  Z   <- list()
  O   <- list()
  MAX <- floor(length(varKeep)/2)
  MAX <- min(MAX,14)
  
  for(i in c(2:MAX)){
    print(i)
    EFA  <- tryCatch(
      factanal(covmat=Ssmooth, factors=i, rotation="promax", lower = 0.01)
      , error = function(e) NULL)
    
    if(!is.null(EFA)){
      
    loads <-
      as.data.frame(EFA$loadings[,1:i]) %>%
      mutate(across(matches("Factor"), function(x) ifelse(x<abs(.35),"",round(x,2))))

    # print(EFA$criteria[1]); # print(EFA$dof[1]); # print(loads)
    
    f <- list()
    for(j in c(1:i)){
      if(length(rownames(loads)[loads[colnames(loads)[j]]!=''])<2) {
        f[[j]] <- paste0("F",j," =~ ", 
                         rownames(loads)[loads[colnames(loads)[j]]!=''], "; ",
                         rownames(loads)[loads[colnames(loads)[j]]!=''], "~~0*",
                         rownames(loads)[loads[colnames(loads)[j]]!=''])
      } else if(length(rownames(loads)[loads[colnames(loads)[j]]!=''])<3) {
        f[[j]] <- paste0("F",j," =~ ",paste0("a*", rownames(loads)[loads[colnames(loads)[j]]!=''],collapse=" + "))
      } else {
        f[[j]] <- paste0("F",j," =~ ",paste0(rownames(loads)[loads[colnames(loads)[j]]!=''],collapse=" + "))
      }
    }
    
    #Specify the Genomic confirmatory factor model
    CFAofEFA <- paste0(f, collapse="; ")
    
    Z[[i]] <-
      tryCatch(
        usermodel(
          Y, 
          estimation = "DWLS", 
          model = CFAofEFA, 
          CFIcalc = TRUE, 
          std.lv = TRUE, 
          imp_cov = FALSE,
        ), error = function(e) NULL)
  
    O[[i]]  <- tryCatch(print(Z[[i]]$modelfit), error = function(e) NULL)
    
    } else {
      
     Z[[i]] <- NULL
     O[[i]] <- NULL
     
    }
  }
  
  saveRDS(list(c(C,Z,O)),paste0("/cluster/projects/p805/crayner/vQTL/GenSem/MvLDSC_",SsSet,"_factorModels.rds"))
  saveRDS(list(c(C,Z,O)),paste0("/cluster/projects/p805/crayner/vQTL/GenSem/MvLDSC_",SsSet,"_factorModels.rds"))
  
  runHDL <- FALSE 
  if(
    !file.exists(paste0("/cluster/projects/p805/crayner/vQTL/GenSem/mvHDL_",SsSet,".rds"))
    && runHDL==TRUE ){
  
    HDLoutput <-
      hdl(
        traits          = files,
        trait.names     = names,
        sample.prev     = sample.prev,
        population.prev = population.prev,
        LD.path         = hd,
        method          = "piecewise"
      )
    saveRDS(HDLoutput,paste0("/cluster/projects/p805/crayner/vQTL/GenSem/mvHDL_",SsSet,".rds"))
  }
# }

cat("\n*** Done ***\n")
stopCluster(cl) 
cat(paste0("\n *** Analysis started at: ",as.character(T0)," \n"))
cat(paste0("\n *** Analysis completed at: ",as.character(Sys.time())," \n"))
quit(save="no")
