## Tandem repeats
library(Biostrings)
library(MotifDb)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm9)
library(ggplot2)
UCSCgenes<- import("~/Scripts/March/mm9.bed")
promoters<- promoters(UCSCgenes)
Arx6Mer<- MotifDb::query(MotifDb, "Arx")[[6]]
Arx6Mer<- round(Arx6Mer*100)
Arx6MerPerfect<-matchPWM(Arx6Mer, BSgenome.Mmusculus.UCSC.mm9, "100%")
Arx6MerTFBS<-matchPWM(Arx6Mer, BSgenome.Mmusculus.UCSC.mm9, "90%")
subset(Arx6mer, strand=="+")

distancetoNextArxMotif<-distanceToNearest(Arx6MerTFBS)
interger<-subset(distancetoNextArxMotif, distance<=00)%>%countLnodeHits()
motifsWithin200bpEachOther<-subset(Arx6MerTFBS, interger)
clusterDistanceToPromtoer<-distanceToNearest(motifsWithin200bpEachOther, promoters )%>%as.data.frame()
distancetoNextArxMotif%>%as.data.frame()%>%
ggplot(aes(x=`distance`))+
  geom_histogram(binwidth = 10000)+
  xlim(2000000)+
  theme_bw()

dim(subset(clusterDistanceToPromtoer, distance<10000))


distanceTogenes<- distanceToNearest(motifsWithin200bpEachOther, UCSCgenes)
numberOfTandemsInsideGenes<-subset(distanceToGenePromoter, distance<=0)


##Average Distance
Numeric2<-apply(distancetoNextArxMotif, 2, as.numeric)
Avg<-sum(Numeric2[,3])/dim(Numeric2)[1]
