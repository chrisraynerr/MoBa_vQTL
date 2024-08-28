#!/usr/bin/Rscript
cat("\n ### Preparing workspace \n")
remove(list = ls()); start.time <- Sys.time()

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

.libPaths(c("/cluster/projects/p805/crayner/software/R/4.0","/cluster/projects/p805/software/R/4.0"))
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... loading packages \n")
packages <- c("optparse","data.table","purrr","bigsnpr","dplyr","R.utils","runonce","stringr",
              "parallel","doParallel","foreach","ggplot2","ggbeeswarm")

installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages==F)){install.packages(packages[!installed_packages],dependencies=T)}
suppressPackageStartupMessages(invisible(lapply(packages, library, character.only=T)))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n... loading options  \n")

options(bitmapType="cairo")

option_list = list(
  make_option("--PheDir",action="store",default=NA,type="character")
  # make_option("--Phe",action="store",default=NA,type="character")
)
opt = parse_args(OptionParser(option_list=option_list))

# PHE <- "ChAnth"
# PHE <- "ChPsy"
PHE  <- opt$PheDir

outDir  <- paste0("Plots_Snps/",PHE)
dir.create(paste0(outDir))

tabDir  <- paste0("Tables_Snps/",PHE)
dir.create(paste0(tabDir))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# TabPath <- list.files(paste0("Plots_QQ/",PHE,"/"), pattern="TopSnps.tsv",full.names=T)

TabPath <- list.files(tabDir, pattern="TopSnps.tsv",full.names=T)
TabPath <- TabPath[!grepl(x=TabPath,"ChWeight_[0-9][0-9][a-z][a-z]_TopSnps.tsv")]
SnpTab  <- lapply(TabPath,fread,data.table=F)
# Trait   <- str_remove_all(basename(TabPath),"_TopSnps.tsv")
# for(i in seq_along(TabPath)){SnpTab[[i]]$Trait}
SnpTab  <- Filter(function(x) dim(x)[1] > 0,SnpTab)
SnpTab  <- do.call("bind_rows",SnpTab)

TopSnps <- SnpTab %>% filter(vP<.00001|mP<.00001)%>%rowwise() %>% 
           mutate(P=-log10(min(vP,mP))) %>% ungroup()%>% 
           select(rsid,P,vP,mP)%>%distinct(rsid,.keep_all=T)

GenObj  <- snp_attach(paste0("/cluster/projects/p805/crayner/data/Gen/MoBaTriosNindv207283.rds"))
SnpInd  <- which(GenObj$map$marker.ID %in% TopSnps$rsid)
GenObj2 <- snp_subset(GenObj,ind.col=SnpInd,backingfile=paste0("/cluster/projects/p805/crayner/vQTL/Data/",PHE,"_TopSnps.GenObj"))
GenObj2 <- snp_attach(paste0("/cluster/projects/p805/crayner/vQTL/Data/",PHE,"_TopSnps.GenObj.rds"))

vClumped <- snp_clumping(GenObj2$genotypes,S=-log10(TopSnps$vP),infos.chr=GenObj2$map$chromosome,infos.pos=GenObj2$map$physical.pos,thr.r2=0.1)
mClumped <- snp_clumping(GenObj2$genotypes,S=-log10(TopSnps$mP),infos.chr=GenObj2$map$chromosome,infos.pos=GenObj2$map$physical.pos,thr.r2=0.1)
Clumped  <- c(vClumped,mClumped)
Snps    <- GenObj2$map$marker.ID[Clumped]
CluSnps <- TopSnps[Clumped,]

TopSnps <- SnpTab %>% filter(vP<.00001|mP<.00001) %>% rowwise() %>% 
  mutate(P=-log10(min(vP,mP))) %>% ungroup()%>% filter(rsid %in% CluSnps$rsid)

fwrite(TopSnps, paste0("/cluster/projects/p805/crayner/vQTL/Tables_Snps/",PHE,"/TopSnps_Clumped.tsv"), sep="\t")

# TopSnps$Trait<- str_replace_all(TopSnps$Trait,"([a-z])([A-Z])","\\1 \\2")
# TopSnps      <- TopSnps %>%tidyr::separate("Trait",c("Person","Scale","Subscale","Age"))

# TopSnps <- TopSnps %>% 
#   select(-Person)  %>% 
#   mutate(Age=str_replace_all(Age,"6mo","Q04")) %>% 
#   mutate(Age=str_replace_all(Age,"18mo","Q05")) %>% 
#   mutate(Age=str_replace_all(Age,"3yr","Q06")) %>% 
#   mutate(Age=str_replace_all(Age,"5yr","Q07")) %>% 
#   mutate(Age=str_replace_all(Age,"8yr","Q09")) %>% 
#   mutate(Age=str_replace_all(Age,"C_14yr","Q11C")) %>% 
#   mutate(Age=str_replace_all(Age,"M_14yr","Q11M")) %>% 
#   mutate(Age=str_replace_all(Age,"3yo","Q06")) %>% 
#   mutate(Age=str_replace_all(Age,"5yo","Q07")) %>% 
#   mutate(Age=str_replace_all(Age,"8yo","Q09")) %>% 
#   mutate(Age=str_replace_all(Age,"7yo","Q08")) %>% 
#   mutate(Age=str_replace_all(Age,"14yo","Q11C")) %>% 
#   mutate(Trait=paste0(Scale,Subscale,Age))

TopPhe  <- unique(TopSnps$Trait)
pDf     <- fread(paste0("Data/",PHE,".tsv")) %>% select(FID,IID,matches(TopPhe))
ID      <- which(GenObj2$fam$sample.ID %in%  pDf$IID)
cDf     <- fread("/cluster/projects/p805/crayner/data/Cov/20SnpPc3GenoBatchCentrePcN207283.tsv") %>% filter(IID %in%  pDf$IID)
cVar    <- names(cDf)[-c(1,2)]
pVar    <- names(pDf)[-c(1,2)]
pDf     <- pDf %>% inner_join(cDf,by=c("FID","IID")) 

# cDf     <- cDf %>% select(cVar)
# pDf     <- pDf %>% select(-cVar)
# pDf     <- pDf %>% mutate(across(pVar), function(x) residuals(lm(as.formula(paste0(x, "~", paste0(cVar,collapse="+"))),data=pDf)))

for(i in seq_along(pVar)){
  pDf[[paste0(pVar[i])]] <- residuals(lm(as.formula(paste0(pVar[i], "~", paste0(cVar,collapse="+"))),data=pDf, na.action=na.exclude))
}

pDf     <- pDf %>% select(-all_of(cVar))

GenObj3 <- snp_subset(GenObj2,ind.row=ID,backingfile=paste0("/cluster/projects/p805/crayner/vQTL/Data/",PHE,"_TopSnps.GenObj_2"))
GenObj3 <- snp_attach(paste0("/cluster/projects/p805/crayner/vQTL/Data/",PHE,"_TopSnps.GenObj_2.rds"))

gDf     <- 
  as_tibble(GenObj3$genotypes[,1:ncol(GenObj3$genotypes)]) %>% 
  purrr::set_names(GenObj3$map$marker.ID) %>% 
  mutate(across(everything(), function(x) round(x, 0)))

gDf2    <- 
  GenObj3$fam %>% select(FID=family.ID,IID=sample.ID,sex) %>% 
  left_join(pDf, by=c("FID","IID")) %>% cbind(gDf)

# names(gDf2) <- str_remove_all(names(gDf2),"Ch")

################################################################################
SnpPlot <- function(data, outcome, snp){
  pals <- sample(1:length(hcl.pals(type = "Qualitative")),1)
  pal_name <- c(hcl.pals()[c(pals)])
  pal_cols <-  hcl.colors(2,palette=paste0(pal_name))
  #
  test <- data %>% select(outcome=all_of(outcome),snp=all_of(snp),sex) %>%
    # mutate(snp = factor(snp, levels = c(0,1,2))) %>%
    mutate(sex = factor(sex, levels = c(1,2), labels = c("males","females"))) %>%
    na.omit()
  #
  # table(test %>% select(snp))
  # B_m <- lm(outcome ~ snp, data = test)
  # B_v <- car::leveneTest(outcome ~ as.factor(snp), data = test)
  mP   <- SnpTab %>% filter(rsid==paste0(snp)&Trait==paste0(outcome)) %>% .$mP
  vP   <- SnpTab %>% filter(rsid==paste0(snp)&Trait==paste0(outcome)) %>% .$vP
  chr  <- SnpTab %>% filter(rsid==paste0(snp)&Trait==paste0(outcome)) %>% .$chr
  maf  <- SnpTab %>% filter(rsid==paste0(snp)&Trait==paste0(outcome)) %>% .$af
  n    <- SnpTab %>% filter(rsid==paste0(snp)&Trait==paste0(outcome)) %>% .$n
  mlab <- paste("P[M] == ",formatC(mP,format="e",digits=2))
  vlab <- paste("P[V] == ",formatC(vP,format="e",digits=2))
  
  capt <- c(
    paste("Chromosome ",as.character(chr)),
    paste("MAF =",as.character(maf)),
    paste("N =",as.character(n)),
    paste("atop(",mlab,")"),
    paste("atop(",vlab,")")
    )
  
  yMax <- max(test$outcome)
  SD   <- sd(test$outcome,na.rm=T)
  
  plot <-
    ggplot(data=test,aes(x=as.factor(snp),y=outcome,colour=sex)) +
    geom_quasirandom(shape = 21, size = 2, alpha = 1, show.legend = F)+
    geom_boxplot(outlier.size = -1, alpha = 0, show.legend = F, colour = "#39ff14", width = 0.25) +
    scale_color_manual(values = c(pal_cols)) +
    ylab(paste0(outcome))+
    xlab(paste0(snp," (minor allele count)")) + # labs(title = capt)+
    theme_bw() +
    theme(
      legend.position="none" #,panel.border=element_blank(),
      #panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank()
      )+
    annotate(geom="text",size=2,hjust=0,
             label=paste0("Chromosome ",as.character(chr),"; MAF =",as.character(round(maf,2)),"; N =",as.character(n)),
             parse=F,x=0.5,y=yMax+(.5*SD))+
    annotate(geom="text",size=2,hjust=0,label=paste("atop(",mlab,")"),parse=T,x=0.5,y=yMax-(.5*SD))+
    annotate(geom="text",size=2,hjust=0,label=paste("atop(",vlab,")"),parse=T,x=0.5,y=yMax-(1*SD))
    
  return(plot)
}
################################################################################

SnpList   <- TopSnps$rsid
TraitList <- TopSnps$Trait

foreach(i=seq_along(SnpList)) %do% {
  plot <- SnpPlot(data=gDf2, outcome=TraitList[[i]], snp=SnpList[[i]])
  ggsave(
    plot,
    file = paste0("/cluster/projects/p805/crayner/vQTL/Plots_Snps/",PHE,"/",TraitList[[i]],"_",SnpList[[i]],"_Fsc_X.png"),
    device = "png",
    dpi = 320,
    width = 1200,
    height = 1200,
    units = "px"
  )
}

