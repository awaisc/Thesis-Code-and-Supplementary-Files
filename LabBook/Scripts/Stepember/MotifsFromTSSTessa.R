


library(tidyr)
library(tibble)
arx6merTFBS<-readRDS("~/DataFiles/ChIPseq/Mouse/ARX6mermm9Sites")
ARXTandem2SpacedTFBS<-readRDS("~/DataFiles/ChIPseq/Mouse/ARXTande2SpacedSites")
Plaindromic4SpacedTFBS<-readRDS("~/DataFiles/ChIPseq/Mouse/Plaindromic4SpacedTFBS")
JolmaTFBS<-readRDS("~/DataFiles/ChIPseq/Mouse/JolmaTFBS")
###Inputs!
enhancerGrange <-
  import(con = "~/DataFiles/Enhancer Tracks/Mouse/mouse_permissive_enhancers_phase_1_and_2.bed")
UCSCgenes <- import("~/DataFiles/Gene Tracks/Mouse/mm9.bed")
startSites<-subset(import("~/DataFiles/Gene Tracks/Mouse/FullMm9genome.GTF"), type== "start_codon")
promoters <- promoters(UCSCgenes,upstream = 5000)
genome<-BSgenome.Mmusculus.UCSC.mm9

ARXMotifModelList<-c("6 Mer"= arx6merTFBS,
                     "Tandem 2 Spaced" = ARXTandem2SpacedTFBS,
                     "Palindromic 4 Spaced" = Plaindromic4SpacedTFBS,
                     "Jolma Model" = JolmaTFBS)

MotifsInPromoters<-lapply(ARXMotifModelList,function(x){subsetByOverlaps(x,promoters)})


ArxDistanceTOTSSInPromoters<-lapply(MotifsInPromoters,function(x){distanceToNearest(x, startSites)%>%as.data.frame()})


allMotifModelsSingleDataFrame<-do.call(rbind.data.frame,ArxDistanceTOTSSInPromoters)
allMotifModelsSingleDataFrame<-rownames_to_column(allMotifModelsSingleDataFrame, var= "Motif Model")

allMotifModelsSingleDataFrameWithModels<-separate(allMotifModelsSingleDataFrame, `Motif Model`,into = c("Motif Model", "Numbers"), sep = '\\.')
allMotifModelsSingleDataFrameWithModels<-subset(allMotifModelsSingleDataFrameWithModels, distance<=2000)
ggplot(allMotifModelsSingleDataFrameWithModels, aes(x=distance, color= `Motif Model`))+
  geom_freqpoly(bins=50)+
  facet_wrap(~`Motif Model`, scales = "free")+
  theme_bw()
