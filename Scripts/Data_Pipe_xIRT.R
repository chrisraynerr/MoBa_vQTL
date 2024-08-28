#!/usr/bin/Rscript
cat("\n### Preparing workspace \n")
remove(list = ls()); start.time <- Sys.time()

setwd("N:/durable/projects/crayner/PipePsy/")
source("SetupLocalRprofile_PipePsy.R")
source("Functions_PhenoData.R")
revcode <- function(x) {len <- length(na.omit(unique(x)))+1;return((x*-1) + len)}

############################################################
cat("\n ... loading packages \n")

packages <- c(
  "optparse", "R.utils","data.table", "magrittr", "dplyr","tidyr",
  "stringr","rlang","purrr",  "mirt","mice","ggmirt"
)

installed_packages <-  packages %in% rownames(installed.packages())
if(any(installed_packages==F)){install.packages(packages[!installed_packages],dependencies=T)}

suppressPackageStartupMessages(invisible(lapply(packages,library,character.only=T)))

option_list = list(make_option(
    "--scale",action="store",default=NA,type="character",help="Scale label"
  ))

opt   =  parse_args(OptionParser(option_list=option_list))
SCALE <- opt$scale

# SCALE <- "ICQ"
# SCALE <- "SMFQ"
# SCALE <- "SWLS"
# SCALE <- "SCL"
# SCALE <- "RSS"
SCALE <- "NHPIC"
SCALE <- "SWLS"
###############################################################################
cat("\n### Loading data dictionary \n")


datd <- 
  fread("DataDictionary_ForPipe.csv") %>%  
  filter(measure==paste0(SCALE)) %>% 
  filter(item_code!="EE867")

codeD <- 
  datd %>%
  select(item_code,scale_item,wave,person) %>%
  distinct() %>%
  na.omit()

codeV  <- codeD$item_code
itemV  <- codeD$scale_item
waveV  <- codeD$wave

persV  <- unique(datd$person)
subsV  <- unique(datd$subscale)
UwaveV <- sort(unique(waveV))
if(SCALE=="ASQ"){subsV  <- str_remove_all(subsV, "Full")}
subsV  <- str_remove_all(subsV, "[[:punct:]]| ")
subsV  <- subsV[nzchar(subsV)]
respV  <- unique(datd$responses)[1]
rawMax <- str_count(respV, ";")+1

cat(" *** Analysing: ", SCALE, "\n")
cat(" *** There are ", length(codeV), " items\n")
cat(" *** There are ", length(UwaveV), " waves (time-points):\n")
cat("\t", UwaveV,"\n")
cat(" *** There are ", length(unique(persV)), " subjects:\n")
cat("\t", unique(persV),"\n")
cat(" *** There are ", length(unique(subsV)), " subscales:\n")
cat("\t", unique(subsV),"\n")

################################################################################
cat("\n### Loading data \n")

personCode <- str_extract(persV,"^.{2}")
scaleCode <- str_to_sentence(SCALE)

#P<-1;P<-2

for(P in c(1:length(personCode))){
  
  cat(personCode[P], "\n")
  if(!file.exists(paste0("Reports/",personCode[P],scaleCode,"_xIRT.html"))){
  if(!file.exists(paste0("Models_xIRT/",personCode[P],scaleCode,"_xIRTModelFitBest.tsv"))){
           
    cat("\n ")
    cat(" ... loading data\n")
    
    datd2  <- datd %>% filter(grepl(x=person, pattern=paste0(personCode[P]), ignore.case = T))
    codeD  <- datd2 %>% select(item_code,scale_item,wave,person) %>% distinct() %>% na.omit()
    codeV  <- codeD$item_code
    itemV  <- codeD$scale_item
    waveV  <- codeD$wave
    UwaveV <- sort(unique(waveV))
    subsV  <- unique(datd2$subscale)
    subsV  <- str_remove_all(subsV, "[[:punct:]]| ")
    if(SCALE=="ASQ"){subsV  <- str_remove_all(subsV, "Full")}
    subsV  <- subsV[nzchar(subsV)]
    subsV  <- sort(subsV)
    subsL  <- list()
    
    for(i in c(1:length(subsV))){
      subsL[[i]] <-
        datd2 %>%
        filter(grepl(pattern=paste0(subsV[i]), x=subscale)) %>%
        select(item_code,scale_item)
    }
    subsVrm <- NULL
    for(i in c(1:length(subsV))){
      if(nrow(subsL[[i]])<3){
        cat(" ... less than three items on the subscale \n")
        cat(" ... will not run IRT on this subscale \n")
        subsVrm <- c(subsVrm,i)
      }}
    
    if(!is.null(subsVrm)){
      subsV <- subsV[-subsVrm]
    }
  
    dat <-
      fread(
        paste0("Data/",personCode[P],scaleCode[1],"_ScaleItems.csv")
      ) %>%
      filter(if_any(matches("Q"), ~ !is.na(.))) %>%
      distinct(PrID, .keep_all = T) %>% 
      mutate(across(is.integer,as.factor))
    
    cat("\n ... impute missing at random \n")
    impDf0   <- mice(dat, maxit=0)
    impMeth  <- impDf0$method
    impPred  <- impDf0$predictorMatrix
    # impPred[,"MoID"] <- -2
    
    cat("\n ... running single imputation \n")
    impDf <- mice(dat,meth=impMeth,pred=impPred,m=1,seed=2019,printFlag=F)
    dat   <- complete(impDf)
    dat   <- dat %>% mutate(across(is.factor,as.integer))
    cat(" ... separate questionnaire waves \n")
  
  datQ <- list()
  
  for(W in c(1:length(UwaveV))){
    datQ[[W]] <- dat %>% select(matches("ID",ignore.case=F),matches(paste0("_",UwaveV[W])))
  }
  cat(" ... subsetting items based on sub-scales\n")
  
  subs <- list()
  
  for(i in c(1:length(subsV))){
    subs[[i]] <- list()
    for(W in c(1:length(UwaveV))){
      subs[[i]][[W]] <-  
        datQ[[W]]  %>% 
        select(matches("ID"),matches(paste0(subsL[[i]]$scale_item))) %>% 
        filter(if_any(matches(paste0(subsL[[i]]$scale_item)), ~ !is.na(.))) %>% 
        as.data.frame()
    
      cat(" ... checking for reverse coded items \n\n")
      
      cordat <-
        subs[[i]][[W]] %>% select(matches(paste0(subsL[[i]]$scale_item))) %>% 
        as.data.frame() %>% mutate(across(everything(),as.integer))
      
      if(ncol(cordat)>2){
      cors      <- cormat(cordat)
      neg_table <- cortab(cors) %>% filter(r < -0.1)
      
      if(nrow(neg_table)>0){
        neg_vars  <- data.frame("V1"=c(as.character(neg_table$Var1),as.character(neg_table$Var2))) %>% 
          group_by(V1) %>% mutate(n=n()) %>% ungroup() %>% distinct() %>% filter(n>=max(n))
        
        subs[[i]][[W]] <- subs[[i]][[W]] %>% mutate(across(matches(neg_vars$V1), revcode)) %>% as.data.frame()
        
        cordat <-
          subs[[i]][[W]] %>% select(matches(paste0(subsL[[i]]$scale_item))) %>% 
          as.data.frame() %>% mutate(across(everything(),as.integer))
        
        cors      <- cormat(cordat)
        neg_table <- cortab(cors) %>% filter(r < -0.1)
        
        if(nrow(neg_table)>0){
          neg_vars  <- data.frame("V1"=c(as.character(neg_table$Var1),as.character(neg_table$Var2))) %>% 
            group_by(V1) %>% mutate(n=n()) %>% ungroup() %>% distinct() %>% filter(n>=max(n))
          
          subs[[i]][[W]] <- subs[[i]][[W]] %>% mutate(across(matches(neg_vars$V1), revcode)) %>% as.data.frame()
          
          cordat <-
            subs[[i]][[W]] %>% select(matches(paste0(subsL[[i]]$scale_item))) %>% 
            as.data.frame() %>% mutate(across(everything(),as.integer))
          
          cors      <- cormat(cordat)
          
        }
      }
    } else {
      cat(" ... less than three items on the subscale \n")
      cat(" ... will not run IRT on this subscale \n")
      subs[[i]][[W]] <- NA
      }  
    }
  }
  # for(i in c(1:length(subs))){
  #   subs[[i]] <- Filter(Negate(is.null), subs[[i]])
  # }
  
cortab(cors)
################################################################################
cat("\n### IRT models \n\n")

IRTitems<- list()
LEVELS  <- list()

for(i in c(1:length(subsV))){
  IRTitems[[i]]  <- list()
  LEVELS[[i]] <- list()
  for(W in c(1:length(UwaveV))){
    if(is.data.frame(subs[[i]][[W]])){
    IRTitems[[i]][[W]] <- subs[[i]][[W]] %>% select(matches(SCALE))
    LEVELS[[i]][[W]]   <- nrow(lapply(IRTitems[[i]][[W]], function(x)unique(x[!is.na(x)])) %>% do.call("bind_cols",.))
    } else {
      IRTitems[[i]][[W]] <- NA 
      LEVELS[[i]][[W]]   <- NA
    }
  }
}

model1 <- list()
model2 <- list() 

TestBest  <- list()
ModelFit1 <- list()
ModelFit2 <- list()

for(i in c(1:length(subsV))){
  
  model1[[i]] <- list() # graded Rasch
  model2[[i]] <- list() # nominal `3`
  
  TestBest[[i]]  <- list()
  ModelFit1[[i]] <- list()
  ModelFit2[[i]] <- list()

  for(W in c(1:length(UwaveV))){
    
    # if(!file.exists(paste0("Models_xIRT/",personCode[P],scaleCode,subsV[i],UwaveV[W],"_xIRTModel.rds"))){
    
      cat("\n ")
      cat(" ... estimating models \n")
    
    if(is.data.frame(subs[[i]][[W]]) && LEVELS[[i]][[W]]>2){
    
    model1[[i]][[W]] <- 
      (mirt(
        data=IRTitems[[i]][[W]],
        model=1,
        verbose=F,
        itemtype="graded",
        SE=T,
        SE.type="sandwich",
        method = "EM",
        technical = list(NCYCLES = 10000)
      ))
    
    model2[[i]][[W]] <- 
      (mirt(
        data=IRTitems[[i]][[W]],
        model=1,
        verbose=F,
        itemtype="nominal",
        SE=T,
        SE.type="sandwich",
        method = "EM",
        technical = list(NCYCLES = 10000)
      ))
    
    } else if(is.data.frame(subs[[i]][[W]]) && LEVELS[[i]][[W]]<3) {
        
        model1[[i]][[W]] <- 
          (mirt(
            data=IRTitems[[i]][[W]],
            model=1,
            verbose=F,
            itemtype="Rasch",
            SE=T,
            SE.type="sandwich",
            technical = list(NCYCLES = 10000)
          ))
        
        model2[[i]][[W]] <- 
          (mirt(
            data=IRTitems[[i]][[W]],
            model=1,
            verbose=F,
            itemtype="3PL",
            SE=T,
            SE.type="sandwich",
            technical = list(NCYCLES = 10000)
          ))
        
      } else {
      model1[[i]][[W]] <- NA
      model2[[i]][[W]] <- NA
    }
      TestBest[[i]][[W]]  <- anova(model1[[i]][[W]], model2[[i]][[W]])
      # ModelFit1[[i]][[W]] <- M2(model1[[i]][[W]], type = "C2")
      # ModelFit2[[i]][[W]] <- M2(model2[[i]][[W]], type = "C2")
  }
}

# M2 similar to chi-square this is the chi-square statistic we obtain from the maximum likelihood statistic
# CFI is the comparative fit index – values can range between 0 and 1 (values greater than 0.90, conservatively 0.95 indicate good fit)
# TLI Tucker Lewis Index which also ranges between 0 and 1 with values greater than 0.90 indicating good fit.

cat(" ... evaluating models \n")

irt_models <- list()
fit_stats  <- list()
best_model <- list()

fitStats1  <- list()
fitStats2  <- list()
  
  for(i in c(1:length(subsV))){
    print(subsV[i])
    fitStats1[[i]] <- list()
    fitStats2[[i]] <- list()
    for(W in c(1:length(UwaveV))){
      if(is.data.frame(subs[[i]][[W]])){
      print(UwaveV[W])
        
        fitStats1[[i]][[W]] <- list()
        fitStats2[[i]][[W]] <- list()
        
        Subscale = subsV[[i]]
        Wave     = UwaveV[W]
        Model    = model1[[i]][[W]]@Model$itemtype[1]
        LL       = model1[[i]][[W]]@Fit$logLik
        AIC      = model1[[i]][[W]]@Fit$AIC
        BIC      = model1[[i]][[W]]@Fit$BIC
        
        fitStats1[[i]][[W]] <- c(Subscale, Wave, Model, round(LL,2), round(AIC,2), round(BIC,2) )
  
        Subscale = subsV[[i]]
        Wave     = UwaveV[W]
        Model    = model2[[i]][[W]]@Model$itemtype[1]
        LL       = model2[[i]][[W]]@Fit$logLik
        AIC      = model2[[i]][[W]]@Fit$AIC
        BIC      = model2[[i]][[W]]@Fit$BIC
        
        fitStats2[[i]][[W]] <- c(Subscale, Wave, Model, round(LL,2), round(AIC,2), round(BIC,2))
      
      } else {
        fitStats1[[i]][[W]] <- NA
        fitStats2[[i]][[W]] <- NA
    }
    }
  }
    
for(i in c(1:length(subs))){
  fitStats1[[i]] <- Filter(Negate(is.na), fitStats1[[i]])
  fitStats2[[i]] <- Filter(Negate(is.na), fitStats2[[i]])
}

fitStats1_2<-list()  
fitStats2_2<-list()  
for(i in c(1:length(subs))){
    fitStats1_2[[i]] <- 
      do.call("rbind", fitStats1[[i]]) %>% 
      as_tibble() %>%
      set_names("Subscale","Wave", "Model", "LL", "AIC", "BIC")
    
    fitStats2_2[[i]] <- 
      do.call("rbind", fitStats2[[i]]) %>% 
      as_tibble() %>%
      set_names("Subscale","Wave", "Model", "LL", "AIC", "BIC")
  
  }
  
  fitStats1_3 <- do.call("rbind", fitStats1_2) %>% as_tibble()
  fitStats2_3 <- do.call("rbind", fitStats2_2) %>% as_tibble()
      
  fit_stats <- 
    bind_rows(fitStats1_3, fitStats2_3) %>% 
    mutate(Y = paste0(Wave,Subscale)) %>% 
    arrange(Y,AIC)
    
  fwrite(fit_stats,
    paste0("Models_xIRT/",personCode[P],scaleCode,"_xIRTModelFit.tsv"
    ))
      
  best_model <- fit_stats %>% distinct(Y, .keep_all = T) %>%  select(-Y)
    
  fwrite(best_model,
    paste0("Models_xIRT/",personCode[P],scaleCode,"_xIRTModelFitBest.tsv"
    ))

  cat("\n ")
  cat(" ... loading best model \n")
  BestModel <- list()
  for(i in c(1:length(subsV))){
    BestModel[[i]] <- list()
    for(W in c(1:length(UwaveV))){
      if(is.data.frame(subs[[i]][[W]])){
        BEST <- best_model %>% filter(Subscale==paste0(subsV[i]) & Wave==paste0(UwaveV[W])) %>% .$Model
        if(BEST=="graded"|BEST=="Rasch"){BEST_MODEL<-"model1"}
        if(BEST=="nominal"|BEST=="3PL"){BEST_MODEL<-"model2"}
        BestModel[[i]][[W]] <- get(BEST_MODEL)[[i]][[W]]
      }  else {
        BestModel[[i]][[W]] <- NA
      } 
    saveRDS(BestModel[[i]][[W]],paste0("Models_xIRT/",personCode[P],scaleCode,subsV[i],UwaveV[W],"_xIRTModel.rds"))
    }
    # BestModel[[i]] <- Filter(is.data.frame, BestModel[[i]])
  }
  }
#}
  
  if(!file.exists(paste0("Data/",personCode[P],scaleCode,subsV,"_xIRTFsc.tsv"))){
    
  cat("\n ")
  cat(" ... loading best model \n")
  
  model1<-model2<-NULL
  BestModel <- list()
  for(i in c(1:length(subsV))){
    BestModel[[i]] <- list()
    for(W in c(1:length(UwaveV))){
    BestModel[[i]][[W]] <- readRDS(paste0("Models_xIRT/",personCode[P],scaleCode,subsV[i],UwaveV[W],"_xIRTModel.rds"))
    }
  }

  cat("\n### Computing factor scores (from the best model) \n")
  
  fScLs <- list()
  fScDf <-list()
  
  for(i in c(1:length(subsV))){
    if(!file.exists(paste0("Data/",personCode[P],scaleCode,subsV[i],"_xIRTFsc.tsv"))){
    fScLs[[i]] <- list()
    for(W in c(1:length(UwaveV))){
      if(is.data.frame(subs[[i]][[W]])){
          
      fScLs[[i]][[W]] <-
        fscores(
          BestModel[[i]][[W]],
          method="EAP",
          full.scores=T,
          full.scores.SE=T,
          MI=1
        )
      
      colnames(fScLs[[i]][[W]]) <- c(paste0(UwaveV[W]),paste0(UwaveV[W],"_SE"))
      ID <- subs[[i]][[W]] %>% select(matches("ID"))# %>% mutate(Wave =paste0(UwaveV[W]) )
      fScLs[[i]][[W]] <- cbind(ID, fScLs[[i]][[W]]) %>% na.omit()
    
      } else {
        fScLs[[i]][[W]] <- NA
      } 
    }
    fScLs[[i]] <- Filter(is.data.frame, fScLs[[i]])
    
    fScDf <- plyr::join_all(fScLs[[i]], by=c("FID","ChIID","MoIID","FaIID","PrID"))
    
      cat("\n ... writing file",i,"of",length(subsV),"\n\n")
      
      fwrite(
        fScDf,
        paste0("Data/",personCode[P],scaleCode,subsV[i],"_xIRTFsc.tsv"),
        sep = "\t",
        col.names = T,
        row.names = F
      )
    }
  }
  }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # rm(Subscale);rm(Wave)
  
    if(!file.exists(paste0("Reports/",personCode[P],scaleCode,"_xIRT.html"))){
      #################################################
      # WRITE FUNCTION TO RENDER REPORT
      render_report = function(SCALE, PERSON) {
        rmarkdown::render(
          "Report_xIRT.Rmd",
          params = list(
            SCALE=SCALE,
            PERSON=PERSON
          ),
          output_file =
            paste0(
              paste0("Reports/",PERSON,SCALE,"_xIRT.html"))
        )
      }

      #################################################
      # RENDER REPORT
      TITLE <- unique(datd2$fullname)[1]

      render_report(
        SCALE = scaleCode,
        PERSON = personCode[P]
      )
    }
  cat("\n ... Done \n")
  }
}


cat("\n*** END OF SCRIPT ***\n")

q(save="no")

