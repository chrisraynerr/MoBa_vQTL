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

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n ### CLUSTER library \n")

.libPaths(c("/cluster/projects/p805/crayner/software/R/4.0",
            "/cluster/projects/p805/software/R/4.0"))

cat("\n... loading packages \n")
packages <- c("optparse","R.utils","GenomicSEM", "bigsnpr", "data.table", "parallel", "doParallel",
              "dplyr", "stringr","foreach")

installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages==F)){install.packages(packages[!installed_packages],dependencies = T)}
suppressPackageStartupMessages(invisible(lapply(packages,library,character.only=T)))


T0 <- Sys.time()
cat(paste0("\n*** Beginning analysis at ",T0," ***\n"))

NCORES    <- nb_cores()
cat(paste0("\n... using ",NCORES," cores \n"))

cat("\n ... registering cluster \n")
Operating <- Sys.info()[['sysname']]
if (Operating != "Windows") {
  cl <- makeCluster(NCORES, type="FORK", outfile="")
} else {
  cl <- makeCluster(NCORES, type="PSOCK", outfile="")
}
registerDoParallel(cl)
on.exit(parallel::stopCluster(cl), add = TRUE)


## ~~~~~~~~~~~~~~~~~~

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

SsDir <- opt$sumstats

# SsDir <- "/cluster/projects/p805/crayner/vQTL/GwaSumStats/ChPsyGenQrsAsr"
# SsDir <- "/cluster/projects/p805/crayner/vQTL/GwaSumStats/ChPsyGenFscAsr"
# SsDir <- "/cluster/projects/p805/crayner/vQTL/GwaSumStats/ChAnth"

## ~~~~~~~~~~~~~~~~~~
## Munge sumstats
## ~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~
## create vector of the summary statistics files
## ~~~~~~~~~~~~~~~~~~
files <- list.files(SsDir,pattern="_RegenieVQt.txt.gz|_RegenieMQt.txt.gz",full.names=T)
# cat(print(files))
## ~~~~~~~~~~~~~~~~~~
## define the reference file being used to allign alleles across summary stats
## ~~~~~~~~~~~~~~~~~~
hm3   <- "/cluster/projects/p805/crayner/data/RefPanels/eur_w_ld_chr/w_hm3.snplist"
## ~~~~~~~~~~~~~~~~~~
## name the traits 
## ~~~~~~~~~~~~~~~~~~
Qcfiles <- files
Qcfiles <- gsub("RegenieVQt.txt.gz","RegenieVQt.Qc",Qcfiles)
Qcfiles <- gsub("RegenieMQt.txt.gz","RegenieMQt.Qc",Qcfiles)

names <- files
names <- gsub("RegenieVQt.txt.gz","RegenieVQt.Munged",names)
names <- gsub("RegenieMQt.txt.gz","RegenieMQt.Munged",names)
names <- gsub("Xirt","",names)

ref   <- fread("/cluster/projects/p805/crayner/data/RefPanels/children_POUN_32413_hm3_ldscores_qc.tsv",
         data.table=F) %>% select(SNP=rsid, MAFSD=sd_val)



for(i in c(1:length(files))){
  if(!file.exists(paste0(names[[i]],".sumstats.gz"))){
    GwaSs <- fread(files[[i]],data.table=F)
    Nrow <- nrow(GwaSs)
    names(GwaSs) <- toupper(names(GwaSs))
    GwaSs <- GwaSs %>%
      dplyr::select(
        CHR,BP=POS,SNP=RSID,A2=A0,A1,MAF=AF,BETA,SE,P,N
      ) %>%
      na.omit() %>%
      inner_join(ref)

      cat("\n ... QC step (checking for unusual effect sizes) \n")
      # Bt
      # sd_ss <- with(GwaSs, 2 / sqrt(N * SE^2 + BETA^2)) #sumstats sd

      # Qt
      sd_y_est <- median(with(GwaSs, MAFSD * SE * sqrt(N)),na.rm = T)
      sd_ss    <- with(GwaSs, sd_y_est / sqrt(N * SE^2))

      is_bad   <- sd_ss<(0.5*GwaSs$MAFSD)|sd_ss>(GwaSs$MAFSD+0.1)|sd_ss<0.1|GwaSs$MAFSD < 0.05

      if(sum(is_bad, na.rm=T) > length(is_bad)*0.5){
        N <- (2/GwaSs$MAFSD)^2 / (GwaSs$SE^2)
        GwaSs$N <- median(N, na.rm = T)
        sd_ss = with(GwaSs, 2 / sqrt(N * SE^2))
        is_bad <- sd_ss < (0.5 * GwaSs$MAFSD) | sd_ss > (GwaSs$MAFSD + 0.1) | sd_ss < 0.1 | GwaSs$MAFSD < 0.05
        cat("More than 50% variants have discordant SD \n")
        cat("Median imputed N: ", median(N, na.rm = T), "\n")
      }
    GwaSs <- GwaSs[!is_bad,] %>% na.omit()
    Loss  <- (Nrow - nrow(GwaSs)) / Nrow
    if(Loss <0.1){
    fwrite(GwaSs, paste0(SsDir,"/",basename(files[[i]])),quote=F,col.names=T,row.names=F,sep="\t",na="NA")
    } else {
    cat(basename(files[[i]]), " have more than 20% variants have discordant SD \n")
      fwrite(GwaSs, paste0(SsDir,"/",basename(Qcfiles[[i]]),".sumstats.gz"),quote=F,col.names=T,row.names=F,sep="\t",na="NA")
    }
  }
}

cat("\n*** Munging summary statistics ***\n")
## ~~~~~~~~~~~~~~~~~~
## run munge
## ~~~~~~~~~~~~~~~~~~
foreach(i = seq_along(1:length(files))) %dopar% {
  if(!file.exists(paste0(names[[i]],".sumstats.gz")) &&
    file.exists(paste0(SsDir,"/",basename(files[[i]])))
  ){
    if(file.exists(paste0(SsDir,"/",basename(Qcfiles[[i]]),".sumstats.gz"))
    ){
    munge(
      files       = paste0(SsDir,"/",basename(Qcfiles[[i]]),".sumstats.gz"),
      hm3         = hm3,
      trait.names = names[[i]],
      )
    } else {
      munge(
        files       = paste0(SsDir,"/",basename(files[[i]])),
        hm3         = hm3,
        trait.names = names[[i]],
      )
    }
  }
}
   
cat("\n*** Done ***\n")

stopCluster(cl) 

cat(paste0("\n *** Analysis started at: ",as.character(T0)," \n"))
cat(paste0("\n *** Analysis completed at: ",as.character(Sys.time())," \n"))

quit(save="no")
