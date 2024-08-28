#!/usr/bin/Rscript
cat("\n*** Script to obtain the quantile integrated rank score\n")
cat("\n### Preparing workspace \n")
remove(list = ls())
start.time <- Sys.time()

suppressMessages(library(data.table))
suppressMessages(library(quantreg))
suppressMessages(library(parallel))
suppressMessages(library(optparse))
suppressMessages(library(quadprog))
suppressMessages(library(progress))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(lmerTest))
suppressMessages(library(lme4))
suppressMessages(library(tidyr))

source("Functions_PhenoData.R")

cat("\n### Defining functions ###\n")
# Step1: fit the null model for k quantiles
GetaTauQi <- function(Qi) {
  if (Qi %% 100 == 0) {
    cat(paste0("... current quantile levels: ", Qi, "/", Qtot, " ...\n"))
  } 
  TauQi  <- Qi / (Qtot + 1)
  # Fit quantile regression between phenotype and covarites at quantile i/k
  Qreg   <- rq(Y ~ DRv, tau = TauQi, method = "fn")
  # Construct the quantile rank score for each quantile
  coeff  <- summary(Qreg, se = "ker")$coefficients
  DRvSE  <- coeff[nrow(coeff), 2]
  aTauQi <- (TauQi - ifelse(residuals(Qreg) < 0, 1, 0)) * DRvSE / (sqrt(-TauQi^2 + TauQi))
  aTau   <- aTauQi
  return(aTau)
}

FitaTauQi <- function(start = 1, end = Qtot) {
  DfQ <- lapply(start:end, GetaTauQi)
  return(DfQ)
}

RankAllPhenos <- function(Dat,Phe) {
  Tmp <- as.data.frame(Dat)
  Tmp <- Tmp[,c("PrID",paste0(Phe))]
  Tmp <- Tmp[!is.na(Tmp[2]), ]
  ID  <- data.frame(PrID=Tmp[,1])
  Y   <- Tmp[,2]
  DRv <- rnorm(length(Y))
  assign( "ID" , ID,  env = .GlobalEnv )
  assign( "Y" ,  Y,   env = .GlobalEnv )
  assign( "DRv", DRv, env = .GlobalEnv )
  # To suppress some Warnings that are harmless: https://github.com/HenrikBengtsson/future/issues/218
  ListQ <- suppressWarnings(FitaTauQi(start=1, end=Qtot))
  IntRankScore <- 0
  for (Qi in 1:Qtot) {
    ## upper quantile weight = 1, lower quantile weight = -1
    weight <- ifelse((Qi > (Qtot) / 2), 1, -1)
    IntRankScore <- IntRankScore + weight * ListQ[[Qi]]
  }
  IntRankScore      <- IntRankScore / ((Qtot) / 2)
  ID[[paste0(Phe)]] <- scale_it(IntRankScore)
  # ID[[paste0(Phe)]] <- IntRankScore
  # ID[[paste0(Phe,"_sc")]] <- scale_it(IntRankScore)
  return(ID)
}

ZSq <- function(Dat, Phe) {
  Tmp <- as.data.frame(Dat)
  Tmp <- Tmp[,c("PrID",paste0(Phe))]
  Tmp <- Tmp[!is.na(Tmp[2]), ]
  ID  <- data.frame(PrID=Tmp[,1])
  Y   <- Tmp[,2]
  Z   <- ((Y - mean(x=Y,na.rm=T))/sd(x=Y,na.rm=T))^2
  ID[[paste0(Phe)]] <- Z
  return(ID)
}

cat("\n### Loading data ###\n")

FILES <- list.files("Data/","xIRTFsc.tsv",full.name=T)
FILES <- FILES[!grepl(x=FILES,"Data/MoRapiFull_xIRTFsc.tsv")]
NAMES <- str_replace_all(FILES,"Fsc","Qrs")
# NAMES <- str_replace_all(NAMES,"Qrs","QrsTest")
TABLES<- str_replace_all(FILES,"Fsc","Cor")
TABLES<- str_replace_all(TABLES,"Data","Tables")
MODELS<- str_replace_all(TABLES,"Cor","Lmm")
MODELS<- str_replace_all(MODELS,".tsv",".doc")

Qtot  <- 1000

for(X in seq_along(FILES[1:2])){
  cat("\n***",FILES[[X]], "\n")
  if(!file.exists(paste0(NAMES[[X]]))){
    Dat    <- data.table::fread(paste0(FILES[[X]]))
    Dat    <- Dat %>% dplyr::select(PrID,matches("^Q")) %>% dplyr::select(!matches("_SE"))
    AllPhe <- names(Dat)[-1]
    Qtot   <- 1000
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cat("### Fitting the quantile regression to obtain quantile rank score ###\n")
    RankList <- NULL
    RankList <- lapply(seq_along(AllPhe), function(i) RankAllPhenos(Dat, AllPhe[i]))
    DatQ     <- plyr::join_all(RankList, by=c("PrID"))
    
    fwrite(DatQ,paste0(NAMES[[X]]),quote=F,col.names=T,row.names=F,sep="\t",na="NA")
    
  cat("*** DONE ! \n\n")
  }
}

