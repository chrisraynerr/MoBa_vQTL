#!/usr/bin/Rscript
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Preparing workspace \n")
remove(list = ls());start.time <- Sys.time()

cat("\n... loading packages \n")
packages <- c(
  "optparse","data.table","ggplot2","stringr","magrittr","colorspace","dplyr",
  "directlabels","extrafont","ggrepel" #,"DescTools"
)
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages==F)){install.packages(packages[!installed_packages],dependencies=T)}
suppressPackageStartupMessages(invisible(lapply(packages,library,character.only=T)))

# extrafont::font_import()
# windowsFonts(sans = windowsFont("Helvetica"))

# Cairo::CairoFonts(
#   regular="Helvetica:style=Regular",
#   bold="Helvetica:style=Bold",
#   italic="Helvetica:style=Italic",
#   bolditalic="Helvetica:style=Bold Italic,BoldItalic",
#   symbol="Symbol", usePUA=TRUE
# )

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n... loading options  \n")

options(bitmapType="cairo")

option_list = list(
  make_option("--PheDir",action="store",default=NA,type="character"),
  make_option("--Phe",action="store",default=NA,type="character")
  # make_option("--GwaSums",action="store",default=NA,type="character"),
  # make_option("--inDir",action="store",default=NA,type="character"),
  # make_option("--outDir",action="store",default=NA,type="character")
)
opt = parse_args(OptionParser(option_list=option_list))

PheName <- opt$Phe
PheDir  <- opt$PheDir
outDir  <- paste0("Plots_QQ/",PheDir)
dir.create(paste0(outDir))
tabDir  <- paste0("Tables_Snps/",PheDir)
dir.create(paste0(tabDir))

vDir    <- paste0("GwaSumStats/",PheDir,"Qrs/")
mDir    <- paste0("GwaSumStats/",PheDir,"/")

# for file in ./*_RegenieMQt.txt.gz; do mv "${file}" "${file//Xirt/}"; done
# vDir    <- "GwaSumStats/ChPsyGenQrsAsr/"
# mDir    <- "GwaSumStats/ChPsyGenFscAsr/"
# vDir    <- "./"
# mDir    <- "./"
# PheName <- "ChCbclShortQ06"
# outDir  <- "Plots_QQ"
# vDir    <- "GwaSumStats/ChAnthQrs/"
# mDir    <- "GwaSumStats/ChAnth/"

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RoundTo <- function(x, multiple = 1, FUN = round) {
  if(is.function(FUN)) {
    # if FUN is a function, then save it under new name and
    # overwrite function name in FUN, which has to be character
    fct <- FUN
    FUN <- "fct"
    FUN <- gettextf("%s", FUN)
  }
  # round will set digits to 0 by default, which is exactly what we need here
  return(eval(parse(text = gettextf("%s(x/multiple) * multiple", FUN))))
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n### Loading GwaSums: ", PheName, "\n")

vDf <- 
  fread(paste0(vDir,PheName,"_RegenieVQt.txt.gz"),data.table=F) %>% 
  mutate(chr=as.numeric(chr),vlogP=log10p,vP=10^-log10p,vB=beta) %>%
  # mutate(chr=as.numeric(chr),vP=10^-log10p) %>%
  select(rsid,chr,af,n,vP,vlogP,vB) %>% 
  na.omit()

mDf <- 
  fread(paste0(mDir,PheName,"_RegenieMQt.txt.gz"),data.table=F) %>% 
  mutate(chr=as.numeric(chr),mlogP=log10p,mP=10^-log10p,mB=beta) %>%
  # mutate(chr=as.numeric(chr),mP=10^-log10p) %>%
  select(rsid,chr,mP,mlogP,mB) %>% 
  na.omit()

Df  <- vDf %>% full_join(mDf,c("rsid","chr")) %>% na.omit()
vDf <- mDf <- NULL
N   <- nrow(Df) 

cat("\n### Make P-P plot \n")

Df$label <- ifelse((Df$vP<1e-5&Df$mP<1e-5)|(Df$vP<5e-8|Df$mP<5e-8),Df$rsid,NA)
Df2 <- Df %>% filter(!is.na(Df$label)) %>% mutate(Trait=paste0(PheName)) %>% select(Trait,rsid,chr,af,n,vP,vB,mP,mB)

fwrite(Df2,file=paste0(tabDir,"/",PheName,"_TopSnps.tsv"),sep="\t")

Df$label <- ifelse((Df$vP<1e-7&Df$mP<1e-7)|(Df$vP<5e-8|Df$mP<5e-8),paste0("Chr",Df$chr,":",Df$rsid),NA)

maxM <- max(RoundTo(max(Df$mlogP),multiple=2,"ceiling"),10)
maxV <- max(RoundTo(max(Df$vlogP),multiple=2,"ceiling"),10)

PpPlot   <-
  ggplot(Df,aes(x=mlogP,y=vlogP,color=chr))+ #,label=rsid
  #ggplot(Df %>% filter(mP<0.001|vP<0.001),aes(x=mlogP,y=vlogP,color=chr))+#,label=rsid))+
  geom_point(size=2,shape=1)+
  scale_colour_viridis_c(option="inferno",direction=-1,begin=0.2,end=0.8)+
  geom_abline(intercept=0,slope=1,alpha=0.5,color="limegreen") +
  geom_hline(yintercept=7.3,alpha=0.5,color="blue",linetype="dashed") +
  geom_vline(xintercept=7.3,alpha=0.5,color="blue",linetype="dashed") +
  scale_x_continuous(expand=c(0,0),limits=c(0,maxM),breaks=c(seq(0,maxM,by=2))) +
  scale_y_continuous(expand=c(0,0),limits=c(0,maxV),breaks=c(seq(0,maxV,by=2))) +
  coord_fixed() +
  geom_abline(intercept=0,slope=1,alpha=0.5,color="limegreen") +
  xlab(expression(paste("mGWA -log"[10], plain(p)))) +
  ylab(expression(paste("vGWA -log"[10], plain(p)))) +
  theme_classic() +
  guides(
    colour=guide_colourbar(ticks=F,nbin=22,limits=c(1,22),
                           breaks=c(1,22),title.position="top",barheight=0.3))+
  labs(colour="Chromosome")+  
  theme(
    text=element_text(family="Sans",size=6),
    panel.border=element_blank(),#plot.title=element_text(hjust=0.5),
    legend.position="top",legend.justification="left",
    legend.direction="horizontal",legend.margin=margin(0,0,0,0,unit="mm"),
    legend.title=element_text(size=6),legend.text=element_text(size=4),
    plot.margin=margin(.1,.1,.1,.1,unit="cm"),legend.box.spacing=unit(0,"pt")
    ) #+
  # geom_text(aes(x=mlogP,y=vlogP,label=label, colour=chr),size=3,vjust=2,hjust=2,family="serif")
  # directlabels::geom_dl(aes(label=label),method=list("smart.grid",cex=.5)) +
  # geom_text_repel(
  #   aes(label=label),
  #   nudge_x=11-subset(Df,!is.na(label))$vlogP,
  #   nudge_y=11-subset(Df,!is.na(label))$mlogP,
  #   # nudge_x=-0.25,direction="y",force=.5,
  #   # size=3,box.padding=.1,point.padding=.5,force=10,
  #   segment.size=0.2,segment.color="grey80") 

ggsave(
  PpPlot,file=paste0(outDir,"/",PheName,"_PpPlot.png"),
  device="png",dpi=320,width=1200,height=1200,units="px"
)
# saveRDS(PpPlot,file=paste0(outDir,"/",PheName,"_PpPlot.rds"))

##############################################################################
cat("\n### Make Q-Q plot \n")

mChisq  <- qchisq(1 - Df$mP, 1)
mLambda <- median(mChisq) / qchisq(0.5, 1)
vChisq  <- qchisq(1 - Df$vP, 1)
vLambda <- median(vChisq) / qchisq(0.5, 1)

E  <- -log10(ppoints(N))
vP <- data.frame(O=-log10(Df$vP),C=Df$chr,GWA="V") %>% arrange(-O) %>% mutate(E=E,LCI=-log10(qbeta(p=(1-.95)/2,shape1=1:N,shape2=N:1)),UCI=-log10(qbeta(p=(1+.95)/2,shape1=1:N,shape2=N:1)))
mP <- data.frame(O=-log10(Df$mP),C=Df$chr,GWA="M") %>% arrange(-O) %>% mutate(E=E,LCI=-log10(qbeta(p=(1-.95)/2,shape1=1:N,shape2=N:1)),UCI=-log10(qbeta(p=(1+.95)/2,shape1=1:N,shape2=N:1)))
P  <- vP %>% bind_rows(mP) %>% mutate(C2=as.factor(ifelse(GWA=="V",C,NA)),GWA=as.factor(GWA),C=as.factor(C)) 
vP <- mP <- NULL

# Make QQ
mlab <- paste0("lambda[M] == ",round(mLambda,2))
vlab <- paste0("lambda[V] == ",round(vLambda,2))
# LABS <- paste0(mlab,"\n",vlab)
LABS <- c(mlab,vlab)

maxM <- max(P[which(P$GWA=="M"),]$O)
maxV <- max(P[which(P$GWA=="V"),]$O)

P$label <- ifelse(P$O>7.28,P$C,ifelse(P$O==maxM,P$C,ifelse(P$O==maxV,P$C,NA)))
# P$label <- ifelse(P$O>7.28,P$C,NA)

maxLab  <- max(RoundTo(max(maxM,maxV),multiple=2,"ceiling"),10)

QqPlot <-
  ggplot(P,aes(x=E,y=O,group=GWA,shape=GWA,color=C)) +
  geom_point(size=2)+
  scale_shape_manual(values=c(3,4))+
  scale_colour_viridis_d(option="inferno",direction=-1,begin=0.2,end=0.8)+
  scale_x_continuous(expand=c(0,0),limits=c(0,8),breaks=c(0,2,4,6,8)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,maxLab),breaks=c(seq(0,maxLab, by=2))) +
  geom_abline(intercept=7.3,slope=0,alpha=0.5,color="blue",linetype="dashed") +
  geom_abline(intercept=5,slope=0,alpha=0.5,color="royalblue",linetype="dotted") +
  geom_abline(intercept=0,slope=1,alpha=0.5,color="limegreen") +
  coord_fixed() +
  xlab(expression(paste("Expected -log"[10], plain(p)))) +
  ylab(expression(paste("Observed -log"[10], plain(p)))) +
  # labs(title=bquote(paste0(mlab,"; ",vlab)), shape="GWA:") +
  # labs(title=bquote(paste0(mlab,"; ",vlab)), shape="GWA:") +
  # ggtitle(expression(paste("Title with ", italic("subscript"))))
  guides(color="none") +
  theme_classic() +
  theme(
    plot.title=element_blank(),# plot.title=element_text(color="grey",size=10,face="italic"),
    text=element_text(family="Sans",size=6),
    panel.border=element_blank(), # plot.title=element_text(hjust=0.5),
    legend.position="top",legend.justification="left",
    legend.direction="horizontal",legend.margin=margin(0,0,0,0,unit="mm"),
    legend.title=element_text(size=6),legend.text=element_text(size=6),
    plot.margin=margin(.1,.1,.1,.1,unit="cm"),legend.box.spacing=unit(0,"pt")
    ) +
  annotate(geom="text",size=2,x=0.5,y=c(maxLab-.5,maxLab-1.5),hjust=0,label=c(LABS),parse=T)

  # annotate(geom="text",size=3,x=c(.1,.1),y=c(9.2,8.8),hjust=0,label=c(mlab,vlab),parse=T) +
  # geom_text_repel(
  #   aes(label=label),
  #   nudge_x=11-subset(Df,!is.na(label))$E,
  #   nudge_y=11-subset(Df,!is.na(label))$O,
  #   size=3,box.padding=.1,point.padding=.5,force=10,
  #   segment.size=0.2,segment.color="grey80")

# saveRDS(QqPlot,file=paste0(outDir,"/",PheName,"_QqPlot.rds"))

ggsave(
  QqPlot,file=paste0(outDir,"/",PheName,"_QqPlot.png"),
  device="png",dpi=320,width=1000,height=1200,units="px"
  )

# FullPlot <- cowplot::plot_grid(QqPlot,PpPlot,align=c("v","h"),ncol=2,axis="b",rel_widths=c(.45,.54))
# 
# ggsave(
#   FullPlot,file=paste0(outDir,"/",PheName,"_FullPlot.png"),
#   device="png",type="cairo",dpi=320,width=2200,height=1200,units="px"
# )



