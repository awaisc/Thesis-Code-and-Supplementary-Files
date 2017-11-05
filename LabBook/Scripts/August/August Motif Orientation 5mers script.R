
library(magrittr)
library(GenomicRanges)
library(ggplot2)
library(magrittr)
library(tibble)
library(pander)
library(reshape2)
library(dplyr)
library(MotifDb)
library(BSgenome.Mmusculus.UCSC.mm9)
library(reshape2)
library(tidyr)


Genes<-import("~/DataFiles/Gene Tracks/Mouse/mm9.bed")
Enhancers<- import("~/DataFiles/Enhancer Tracks/Mouse/mouse_permissive_enhancers_phase_1_and_2.bed")
promoters<-promoters(Genes)
genome<-BSgenome.Mmusculus.UCSC.mm9
TAATTA

PlaindromicNoSpace<-matchPWM(cbind(round(PWM("TAATT")*5),round(PWM("TTAAT")*5)), genome, "100%")  
Plaindroic1Space<-matchPWM(cbind(round(PWM("TAATT")*5),0.25, round(PWM("TTAAT")*5)), genome, "100%")  
Plaindroic2Space<-matchPWM(cbind(round(PWM("TAATT")*5), 0.25, 0.25, round(PWM("TTAAT")*5)), genome, "100%") 


TandemNoSapce<-matchPWM(cbind(round(PWM("TAATT")*5),round(PWM("AATTA")*5)), genome, "100%")  
Tandem1Space<-matchPWM(cbind(round(PWM("TAATT")*5), 0.25, round(PWM("AATTA")*5)), genome, "100%")  
Tandem2Space<-matchPWM(cbind(round(PWM("TAATT")*5), 0.25, 0.25, round(PWM("AATTA")*5)), genome, "100%")  


#jolmaTFBS<-matchPWM(cbind(round(PWM("TAATT")*5), 0.25, round(PWM("AATTA")*5)), genome, "100%")

#Function Of ARX
genomicLocation<-function(x){
  dataFrame1<-rbind.data.frame(
    "Enhancers" =subsetByOverlaps(x, Enhancers)%>%length(),
    
    "Genes" = subsetByOverlaps(x, Genes)%>%length(),
    
    "Promoters"= subsetByOverlaps(x, promoters)%>%length(),
    
    "Other"= length(x)-sum(
      "Enhancers" =subsetByOverlaps(x, Enhancers)%>%length(),
      "Genes" = subsetByOverlaps(x, Genes)%>%length(),
      "Promoters"= subsetByOverlaps(x, promoters)%>%length()))
  
  colnames(dataFrame1) <-c("Number Of Motifs")
  return(dataFrame1)
}

ARX5merModelsMotifSpacing<-c( "Plaindromic No Space" =PlaindromicNoSpace,
                              "Plaindromic One Space" = Plaindroic1Space,
                              "Plaindromic Two Space"= Plaindroic2Space,
                              "Tandem No Space"= TandemNoSapce,
                              "Tandem 1 Space/ Jolma" = Tandem1Space,
                              "Tandem 2 Space" = Tandem2Space)
LocationsOf5mers<-lapply(ARX5merModelsMotifSpacing, genomicLocation)

DataFrameLocationOf5Mers<-do.call(rbind, LocationsOf5mers)%>%rownames_to_column()


dataFrameGgplot<-separate(DataFrameLocationOf5Mers,rowname, into=c("Motif Model", "Genomic Location"), sep= "\\.")

ggplot(dataFrameGgplot, aes(x=`Motif Model`, fill= `Genomic Location`, y= `Number Of Motifs`))+
  geom_bar(stat="identity")+
  theme_bw()+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=90, vjust=0.3))




startSites<-import("~/DataFiles/Gene Tracks/Mouse/FullMm9genome.GTF")%>%subset(. , type=="start_codon")

DistanceToTSS<-lapply(ARX5merModelsMotifSpacing,function(x){
  object1<-distanceToNearest(startSites, x)%>%
  as.data.frame()
  object2<-object1[,3]
  return(object2)
})

DistanceOf5MersDataFrame<-do.call(cbind, DistanceToTSS)
reshapedDistanceOf5mers<-melt(DistanceOf5MersDataFrame)

ggplot(reshapedDistanceOf5mers, aes(color=Var2, x= value))+
  geom_freqpoly(bins=1000)+
  theme_bw()+
  coord_cartesian(xlim = c(0, 200000))+
  xlab(label = "Distance To Transcription Start Site")+
  ylab(label = "Number Of Motifs")+
  guides(color=guide_legend(title= "Motif Model"))



