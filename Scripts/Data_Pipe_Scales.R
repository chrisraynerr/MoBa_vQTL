#!/usr/bin/Rscript
cat("\n### Preparing workspace \n")
remove(list = ls()); start.time <- Sys.time()
setwd("N:/durable/projects/crayner/PipePsy/")
source("Setup_LocalRprofile.R")

############################################################
cat("\n ... loading packages \n")

packages <- c(
  "optparse", "data.table", "dplyr","stringr","tidyr","haven","rlang",
  "phenotools", "magrittr", "purrr")

installed_packages <-  packages %in% rownames(installed.packages())

if(any(installed_packages==F)){install.packages(packages[!installed_packages],
    dependencies = T)}

suppressMessages( suppressWarnings( invisible(
      lapply(packages, library, character.only=T) )))

option_list = list(make_option(
    "--scale", action = "store", default = NA, type = "character", help = "Scale label"
  ))

opt = parse_args(OptionParser(option_list=option_list))

SCALE <- opt$scale

# SCALE <-"CBCL"

################################################################################
cat("\n### Loading data dictionary \n")

cat(" *** The scale is: ", SCALE, "\n")

datd <- 
  fread("DataDictionary_ForPipe.csv") %>%  
  filter(measure==paste0(SCALE)) %>% 
  filter(item_code!="EE867")

persV <- unique(datd$person)

cat(" *** There are", length(persV), "subjects:\n")
cat(persV)

personCode <- str_extract(persV,"^.{2}")
scaleCode  <- str_to_sentence(SCALE)

# P<-2
for(P in c(1:length(personCode))){
  
  if(!file.exists(paste0("Data/",personCode[P],scaleCode,"_ScaleItems.csv"))){
    
  cat("\n### Extracting items for",persV[P],"\n")
  
  codeD <- 
    datd %>%
    select(item_code,scale_item,wave,person,subscale,responses) %>%
    filter(grepl(x=person, pattern=personCode[P]))  %>%
    distinct(item_code, .keep_all = T) %>%
    distinct() %>%
    na.omit()

  testD <- codeD %>% select(item_code,scale_item,wave)%>%distinct()
  codeV <- testD$item_code
  itemV <- testD$scale_item
  waveV <- testD$wave
  
  subsV <- unique(codeD$subscale)
  subsV <- str_remove_all(subsV, "[[:punct:]]| ")
  respV <- unique(codeD$responses)[1]
  
  rawMax  <- str_count(respV, ";")+1
  
  cat(" *** There are ", length(codeV), " items\n\n")
  cat(" *** There are ", length(unique(waveV)), " waves (time-points):\n\n")
  cat(unique(waveV,"\n\n"))
  cat(" *** There are ", length(unique(subsV)), " subscales:\n\n")
  cat(unique(subsV),"\n\n")
  
  subsL <- list()
  for(i in c(1:length(subsV))){
    subsL[[i]] <-
      datd %>%
      filter(grepl(paste0(subsV[i]), subscale)) %>%
      select(item_code,scale_item) %>%
      distinct()
  }
  
  ################################################################################
  cat(" ... getting variables\n\n")
  # ERROR CAUSED WHEN THERE ARE DUPLICATED ITEMS IN THE LIST!!!!!
  # raw <- 
  #   curate_dataset2023(
  #     variables_required=c(codeV),
  #     moba_data_root_dir=paste0(QDIR),
  #     PDB="2601",
  #     moba_data_version=12,
  #     completion_threshold=0.5,
  #     return_items=T,
  #     consistent_items=F,
  #     out_format="merged_df"
  #     )
  
  NAMES <- fread(paste0(QDIR,"PDB2601_colnames.csv"),header=F)
  COLS  <- NAMES[,which(NAMES$V1 %in% c("preg_id","BARN_NR",codeV))]
  raw   <- fread(paste0(QDIR,"PDB2601_allQs.csv"),header=T,data.table=F)[,c(COLS)]
  ################################################################################
  cat(" ... quick QC\n\n")
  
  codeV2 <- codeV[codeV %in% names(raw)]
  itemV2 <- itemV[codeV %in% names(raw)]
  waveV2 <- waveV[codeV %in% names(raw)]
  
  items  <- make.unique(paste0(itemV2,"_",waveV2),sep="_")
  
  wide <- 
    raw %>% 
    zap_formats() %>%
    zap_labels() %>%
    zap_label() %>%
    mutate(PrID=paste0(preg_id,'_',BARN_NR)) %>% 
    select(PrID,all_of(codeV2)) %>%
    # select(PrID,matches("_raw")) %>%
    # rename_all(list(~str_replace_all(., '_raw',''))) %>% 
    rename_at(vars(c(all_of(codeV2))), ~ c(items)) %>%
    filter(if_any(matches(items), ~ !is.na(.))) %>%
    mutate(across(matches(items), ~ as.integer(.))) %>%
    mutate(across(matches(items), ~ na_if(.x, .x > rawMax) ))  %>%
    filter_at(vars(items), all_vars(. <= rawMax | is.na(.))) %>% 
    setDT()
  
  ################################################################################
  cat(" ... linking to geno data\n\n")
  
  Link <- 
    fread("../Data/LinkageFiles/MoBaTriosGeLinkage.txt") %>%
    mutate(PrID=paste0(PrID,"_",BaN)) %>%
    select(PrID,FID,ChIID,MoIID,FaIID)
  
  dat <- 
    wide %>%
    full_join(Link,"PrID")%>%
    select(FID,ChIID,MoIID,FaIID,PrID,matches("Q")) %>%
    #filter(!is.na(FID)) %>%
    filter(if_any(matches("Q"), ~ !is.na(.)))
  
  ################################################################################
  cat(" ... writing file\n\n")
  
  fwrite(
    dat,
    paste0("Data/",personCode[P],scaleCode,"_ScaleItems.csv"),
    col.names = T,
    row.names = F
    )

  }
  
  REPORT <- FALSE
  
  if(REPORT==TRUE && !file.exists(paste0("Reports_Items/",personCode[P],scaleCode,"_ScaleItems.html"))){
  #################################################
  # WRITE FUNCTION TO RENDER REPORT
  render_report = function(SCALE, PERSON) {
    rmarkdown::render(
      "Report_ScaleItems.Rmd",
      params = list(
        SCALE=SCALE,
        PERSON=PERSON
      ),
      output_file =
        paste0(
          paste0("Reports_Items/",PERSON,SCALE,"_ScaleItems.html"))
    )
  }
  
  #################################################
  # RENDER REPORT
  TITLE <- unique(datd$fullname)[1]
  
  render_report(
    SCALE = scaleCode,
    PERSON = personCode[P]
  )
  
  }
}

cat("*** END OF SCRIPT ***\n\n")

q(save="no")

