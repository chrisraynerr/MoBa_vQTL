#!/usr/bin/Rscript

cat("\n### Preparing workspace \n")

remove(list = ls())
start.time <- Sys.time()
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

#.libPaths(c("/cluster/projects/p805/crayner/software/R/4.0",
#            "/cluster/projects/p805/software/R/4.0"))

.libPaths("C:/Users/p805-christopherdr/Documents/R/win-library/4.0")


cat("\n... loading packages \n")
packages <- c("data.table","ggplot2","stringr","magrittr","colorspace","dplyr","tidyr","ggfittext","heatmaply")

installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages==F)){install.packages(packages[!installed_packages],dependencies = T)}
suppressPackageStartupMessages(invisible(lapply(packages,library,character.only=T)))

options(bitmapType="cairo")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plot Multivariable LDSC
setwd("N:/durable/projects/crayner/vQTL/GenSem")

h2Tab    <- 
  fread("mvLDSC_ChPsyFsxAsr_h2FullTable.txt",header=F,sep=" ") %>% 
  select(Trait=V1,h2=V2,se=V3,Int=V4,IntSE=V5,Ratio=V6,RatioSE=V7,Z=V8) %>% 
  #filter(Z > 2) %>% 
  select(Trait, h2, se, Z)  %>% 
  # filter(!(grepl(x=Trait,"Mot_")))%>% 
  distinct(Trait,.keep_all = T) 

mean(h2Tab$h2)
max(h2Tab$h2)
min(h2Tab$h2)

h2Tab_Qrs <- 
  h2Tab %>% 
  filter(grepl(x=Trait,"_Qrs")) %>% 
  mutate(across(matches("Trait"),str_remove_all, "Qrs_")) %>% 
  separate(Trait,c("Sex","Trait","Time"),"_") %>% 
  mutate(Z = ifelse(Z<0,0,Z)) %>% 
  mutate(alpha = ifelse(h2<0 | h2-se<0, 0.4, 1))

h2Tab_Qrs %>% 
  filter(Sex=="M" & Z>=2.5)

h2Tab_Qrs %>% filter(h2>0) %>% summarise(mean(h2))
mean(h2Tab_Qrs$h2)

scaleFUN <- function(x) sprintf("%.2f", x)

h2Plot_Qrs<-
ggplot(h2Tab_Qrs %>% mutate(h2=round(h2,2)), 
       aes(x=Sex, y = h2, colour = Sex, alpha = alpha )) +
  geom_point(aes(size = Z)) +
  geom_errorbar(aes(ymin = h2 - se, ymax = h2 + se), width = 0.2) +
  facet_grid(Trait~Time,scales = "free_y",) +
  theme_minimal()+
  theme(
    legend.position="bottom", 
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    strip.text.y=element_text(angle=0)
    )+ 
  scale_colour_viridis_d(option="plasma")+
  guides(size="none",alpha="none")+
  geom_hline(yintercept=0,colour="limegreen",linetype="dashed")+
  scale_y_continuous(
    labels=scaleFUN
    # labels = function(x){
    #   ifelse(round(min(x),2)==0.00, c(0,round(median(x),2),round(max(x),2)), c(round(min(x),2),0,round(max(x),2)))}, 
    # breaks = function(x){
    #   ifelse(round(min(x),2)==0.00, c(0,round(median(x),2),round(max(x),2)), c(round(min(x),2),0,round(max(x),2)))}, 
    # guide=guide_axis(n.dodge=2,)
    )

ggsave(
  h2Plot_Qrs,
  file = "h2Plot.png",
  device = "png",
  dpi = 320,
  width = 2400,
  height = 10000,
  units = "px"
)

h2TabF_Qrs <- 
  h2Tab %>% 
  filter(grepl(x=Trait1,"F_Qrs")) %>% 
  mutate(across(matches("Trait"),str_remove_all, "F_Qrs_"))
  
h2TabM_Tsc <- 
  h2Tab %>% 
  filter(grepl(x=Trait1,"M_Tsc")) %>% 
  mutate(across(matches("Trait"),str_remove_all, "M_Tsc_"))

h2TabF_Tsc <- 
  h2Tab %>% 
  filter(grepl(x=Trait1,"F_Tsc")) %>% 
  mutate(across(matches("Trait"),str_remove_all, "F_Tsc_"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rgTab  <- 
  fread("ldsc_rg.txt",header=F,sep=" ") %>% 
  select(Trait1=V1,Trait2=V2,rg=V3,se=V4) %>% 
  mutate(across(everything(),str_remove_all, "Ch")) %>% 
  mutate(across(everything(),str_replace_all, "FPsy", "F_")) %>% 
  mutate(across(everything(),str_replace_all, "MPsy", "M_")) %>% 
  mutate(across(everything(),str_replace_all, "TscQ", "_Q")) %>%
  mutate(across(everything(),str_remove_all, "\\(|\\)")) %>% 
  mutate(across(c(rg,se),as.numeric))%>% 
  mutate(rg=ifelse(rg>1,1,rg)) %>% 
  mutate(rg=ifelse(rg < -1,-1,rg)) %>% 
  mutate(se=ifelse(se>1,1,se))%>% 
  mutate(se=ifelse(se < -1,-1,se)) %>% 
  mutate(z=rg/se) %>% 
  filter(Trait1!=Trait2) %>% 
  filter(!(grepl(x=Trait1,"Mot_")|grepl(x=Trait2,"Mot_")))

rgTab2 <- rgTab %>% bind_rows(rgTab %>% select(Trait1=Trait2,Trait2=Trait1,rg,se,z))

rgTabM_Qrs <- 
  rgTab2 %>% 
  filter(grepl(x=Trait1,"M_Qrs")&grepl(x=Trait2,"M_Qrs")) %>% 
  mutate(across(matches("Trait"),str_remove_all, "M_Qrs_")) #%>% 
  # bind_rows(h2TabM_Qrs)

rgTabF_Qrs <- 
  rgTab2 %>% 
  filter(grepl(x=Trait1,"F_Qrs")&grepl(x=Trait2,"F_Qrs")) %>% 
  mutate(across(matches("Trait"),str_remove_all, "F_Qrs_")) #%>% 
  # bind_rows(h2TabF_Qrs)


rgTabM_Qrs2 <- 
  rgTabM_Qrs %>% 
  # mutate(Trait2=ifelse(Trait1==Trait2, "h2", Trait2)) %>% 
  select(Trait1,Trait2,rg) %>% 
  pivot_wider(id_cols = Trait1, names_from=Trait2,values_from = rg) %>% 
  tibble::column_to_rownames("Trait1")

rgTabM_QrsZ <- 
  rgTabM_Qrs %>% 
  # mutate(Trait2=ifelse(Trait1==Trait2, "h2", Trait2)) %>% 
  # mutate(z=rg/se) %>% 
  select(Trait1,Trait2,z) %>% 
  # mutate(Trait2=ifelse(Trait1==Trait2, "h2", Trait2)) %>% 
  pivot_wider(id_cols = Trait1, names_from=Trait2,values_from = z) %>% 
  tibble::column_to_rownames("Trait1")


rgTabF_Qrs2 <- 
  rgTabF_Qrs %>% 
  select(Trait1,Trait2,rg) %>% 
  pivot_wider(id_cols = Trait1, names_from=Trait2,values_from = rg) %>% 
  tibble::column_to_rownames("Trait1")

rgTabF_QrsZ <- 
  rgTabF_Qrs %>% 
  select(Trait1,Trait2,z) %>% 
  pivot_wider(id_cols = Trait1, names_from=Trait2,values_from = z) %>% 
  tibble::column_to_rownames("Trait1")

library(heatmaply)

rgPlotM <- 
heatmaply_cor(
  rgTabM_Qrs2, # %>% select(!h2),
  node_type="scatter",
  point_size_mat=abs(rgTabM_QrsZ),
  label_names=c("x", "y", "rg"),
  point_size_name = "z",
  fontsize_col=8,
  fontsize_row=8,
  dendrogram ="column",
  file=paste0("MQrsRg_heatmaply.html")
)

rgPlotF <- 
  heatmaply_cor(
    rgTabF_Qrs2,
    node_type="scatter",
    point_size_mat=abs(rgTabF_QrsZ),
    label_names=c("x", "y", "rg"),
    point_size_name = "z",
    fontsize_col=8,
    fontsize_row=8,
    dendrogram ="column",
    file=paste0("FQrsRg_heatmaply.html")
  )


h2Tab <- rgTab2 %>% select(h2)
h2TabE <- seTab %>% select(h2)

h2Plot <- 
  heatmaply_cor(
    h2Tab,
    node_type="scatter",
    point_size_mat=abs(-1*h2TabE),
    label_names=c("x", "y", "h2"),
    point_size_name = "-e",
    fontsize_col=8,
    fontsize_row=8,
    Rowv = NA,
    Colv = NA,
    file=paste0("Parth2_heatmaply.html")
    # file=paste0("N:/durable/projects/crayner/data/Plots/participation_pred_heatmaply_plot.html")
  )
         
fig <- subplot(rgPlot, h2Plot, widths = c(0.9, 0.1), shareY = TRUE)

library(knitr)
library(kableExtra)
x <-kable(h2Tab, format = "html")
class(x) <- c("kableExtra", class(x))
print(x)

datatable(
  h2Tab, 
  rownames = FALSE, 
  options = list(
    columnDefs = list(list(className = 'dt-right', targets = '_all')),
    columnDefs = list(list(className = 'dt-head-right', targets = '_all'))
    )
  )
