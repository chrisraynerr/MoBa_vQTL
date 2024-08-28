#!/usr/bin/Rscript
cat("\n### Preparing workspace \n")
remove(list = ls()); start.time <- Sys.time()
setwd("N:/durable/projects/crayner/PipePsy/")
source("Setup_LocalRprofile.R")
############################################################
cat("\n ... loading packages \n")
packages <- c("data.table","dplyr","purrr","stringr","ggplot2","ggstatsplot","tidyr")
installed_packages <-  packages %in% rownames(installed.packages())
if(any(installed_packages==F)){install.packages(packages[!installed_packages],
                                                dependencies = T)}
suppressMessages( suppressWarnings( invisible(
  lapply(packages, library, character.only=T) )))


cat("\n ... loading covariates \n")

Cov <- 
  data.table::fread("N:/durable/projects/crayner/Data/Participation/MoBaAgeVarsLongLinked.csv") 

gID <- Cov %>% mutate(PrID=paste0(PrID,"_",ChN)) %>% select(FID,IID=ChIID,PrID) %>% distinct() %>% filter(IID!=""&IID!=" "&!is.na(IID))
  
Cov <- 
  Cov %>% 
  select(MoLNR,FaLNR,PrID,ChN, ChSex=ChSEX, ChAge, MoAge, FaAge, Wave=Q) %>% 
  mutate(Wave=ifelse(Wave=="QF1","Q01",Wave)) %>% 
  mutate(Wave=ifelse(Wave=="QF2","Q12",Wave)) %>% 
  distinct() %>% filter(ChAge>0) %>% filter(Wave!="Q01"&Wave!="Q03") %>% 
  mutate(PrID2=paste0(PrID,"_",ChN), Wave=droplevels(as.factor(Wave))) %>% 
  pivot_wider(id_cols=c(MoLNR,FaLNR,PrID2,PrID,ChN,ChSex),names_from=Wave,values_from=c(ChAge,MoAge,FaAge))%>% 
  mutate( ChAge_Q11M = ChAge_Q09 + 6,
          MoAge_Q11M = MoAge_Q09 + 6,
          FaAge_Q11M = FaAge_Q09 + 6,
          ChAge_Q11C = ChAge_Q09 + 6,
          MoAge_Q11C = MoAge_Q09 + 6,
          FaAge_Q11C = FaAge_Q09 + 6
  )

cat("\n ... loading data \n")

FILES <- list.files("Data/","xIRTFsc.tsv",full.names=T)
FILES <- FILES[grep(pattern="Ch",x=FILES)]
NAMES <- basename(gsub(x=FILES,pattern="_xIRTFsc.tsv",replacement=""))

DfLs  <- lapply(FILES,fread)

selectCols <- 
  function(Df){ 
    Df <- Df %>% 
      select(PrID,starts_with("Q")) %>% 
      select(!matches("SE")) 
  }

DfLs <- lapply(DfLs, selectCols)
# Define a function to rename columns in a data frame
renameCols <- function(Df, NAMES) {
  Qcols <- grepl("^Q", colnames(Df)) # Identify columns that start with "Q"
  colnames(Df)[Qcols] <- paste0(NAMES, colnames(Df)[Qcols]) # Rename columns
  return(Df)
}
# Use lapply to apply the rename_cols function to each data frame in DfLs
DfLs <- lapply(seq_along(DfLs), function(i) renameCols(DfLs[[i]], NAMES[i]))
# Define the two ID columns to merge on
IDcols <- c("PrID")
# Use Reduce() function to merge all data frames in DfLs_renamed
Df  <- plyr::join_all(DfLs,"full",by=IDcols)
fwrite(Df,"Data/ChPsyFsc.tsv",quote=F,col.names=T,row.names=F,sep="\t",na="NA")
# 
# gDf <- Df %>% inner_join(gID) %>% select(FID,IID,everything()) %>% select(-PrID)
# fwrite(gDf,"Data/ChPsyGenFsc.tsv",quote=F,col.names=T,row.names=F,sep="\t",na="NA")

# fwrite(Df,"Data/FamPsyFsc.tsv",quote=F,col.names=T,row.names=F,sep="\t",na="NA")

DfLs<-NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n ... covariate adjustment \n")

DfList2 <- list()
for ( Q in c(paste0("Q0",4:9),"Q11M","Q11C")){
  DfList2[[paste0(Q)]] <-
    Df %>%
    dplyr::select(PrID2=PrID,matches(paste0(Q))) %>% # FID,IID
    dplyr::inner_join(Cov %>% select(MoLNR,FaLNR,PrID2,ChSex,matches(paste0(Q)))) %>%
    # tidyr::separate(PrID2,c("PrID","ChN")) %>%
    dplyr::mutate(FamID=as.numeric(as.factor(paste0(MoLNR,"_",FaLNR)))) %>%  #,PrID=as.numeric(as.factor(PrID))) %>%
    dplyr::select(-MoLNR,-FaLNR) %>% #,-ChN)
    dplyr::select(PrID=PrID2,FamID,everything())
}

for ( Q in c(paste0("Q0",4:9),"Q11M","Q11C")){
  for ( X in c(names( DfList2[[paste0(Q)]] %>% dplyr::select(!matches("ID|Age|Sex"))) )){
    DfList2[[paste0(Q)]][paste0(X)] <-
      residuals(lm(as.formula(paste0(X, "~ ChSex + ChAge_",Q)),data=DfList2[[paste0(Q)]],na.action = "na.exclude"))
  }
  DfList2[[paste0(Q)]] <- DfList2[[paste0(Q)]] %>% dplyr::select(!matches("Age|Sex"))
}

IDcols <- c("PrID","FamID")

Df2  <- plyr::join_all(DfList2,"full",by=IDcols)
fwrite(Df2,"Data/ChPsyFscAsr.tsv",quote=F,col.names=T,row.names=F,sep="\t",na="NA")

gDf <- Df2 %>% inner_join(gID) %>% select(FID,IID,everything()) %>% select(-PrID,-FamID)
fwrite(gDf,"Data/ChPsyGenFscAsr.tsv",quote=F,col.names=T,row.names=F,sep="\t",na="NA")


# library(lmerTest)
# ModList <- list()
# for ( Q in c(paste0("Q0",4:9))){
#   ModList[[paste0(Q)]] <- list()
#   for ( X in c(names( DfList2[[paste0(Q)]] %>% select(!matches("ID|ChN|Age|Sex|LNR"))) )){
#     ModList[[paste0(Q)]][[paste0(X)]] <-
#       lme4::lmer(as.formula(
#         paste0(X, "~ ChSex + ChAge_",Q," + MoAge_",Q," + FaAge_",Q, "+ (1|FamID)")
#         # + (1|PrID) + (1|MoLNR) + (1|FaLNR)")
#         ),data=DfList2[[paste0(Q)]],na.action="na.exclude")
#   }
# }

cat("\n ... prep datasets for multiple imputation \n")

names(Df) <- str_replace_all(names(Df),"Q","_Q")

Df2 <- 
  Df %>% 
  select(PrID2=PrID,matches(paste0("_Q"))) %>% # FID,IID
  inner_join(Cov %>% select(MoLNR,FaLNR,PrID2,ChSex,matches(paste0("_Q"))), "PrID2") %>% 
  separate(PrID2,c("PrID","ChN")) %>% 
  select(MoLNR,FaLNR,PrID,ChN,ChSex,matches("Age"),everything())

fwrite(Df2,"Data/FamPsyFsc_PreImp.tsv",quote=F,col.names=T,row.names=F,sep="\t",na="NA")

Df3 <- 
  haven::read_sav("N:/durable/data/moba/Original files/Delivery 05 2023/PDB2601_20230513/PDB2601_MFR_v12.sav") %>%
  haven::as_factor() %>% 
  filter(is.na(DAAR)) %>%
  select(!matches("VERSJON"),!matches("ALDER"),!DMND,!DAAR,FODT_MFR,!KJONN) %>% 
  select(PrID=PREG_ID_2601,ChN=BARN_NR,everything()) %>% 
  mutate(PrID=as.character(PrID),ChN=as.character(ChN),FLERFODSEL=ifelse(is.na(FLERFODSEL),"Single birth",FLERFODSEL)) %>% 
  inner_join(Df2, by=c("PrID","ChN")) #  %>% select(PrID,ChN))

fwrite(Df3,"Data/FamPsyFsc_CovPreImp.tsv",quote=F,col.names=T,row.names=F,sep="\t",na="NA")
