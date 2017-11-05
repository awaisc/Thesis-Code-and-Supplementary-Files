

##New code stuff
library(seqLogo)
library(magrittr)
library(GenomicRanges)
library(ggplot2)
library(magrittr)
library(tibble)
library(pander)
library(reshape2)
library(plyr)
library(MotifDb)
library(BSgenome.Mmusculus.UCSC.mm9)

arx6Mer <-
  rbind(
    A = c(0, 1, 1, 0, 0, 1),
    C = c(0, 0, 0, 0, 0, 0),
    G = c(0, 0, 0, 0, 0, 0) ,
    T = c(1, 0, 0, 1, 1, 0)
  )

arx6MerPlaindrome<-
  rbind(
    A = c(1, 0, 0, 1, 1, 0),
    C = c(0, 0, 0, 0, 0, 0),
    G = c(0, 0, 0, 0, 0, 0) ,
    T = c(0, 1, 1, 0, 0, 1)
  )

grange6mer<-matchPWM(arx6Mer, BSgenome.Mmusculus.UCSC.mm9, "100%")
grange6merplaindrome<-matchPWM(arx6MerPlaindrome, BSgenome.Mmusculus.UCSC.mm9, "100%")


grange6merPlus<-subset(grange6mer, strand=="+")
grange6merMinus<-subset(grange6mer, strand=="-")
grange6merPlaindromePlus<-subset(grange6merplaindrome, strand =="+")
grange6merPlainDromeMinus<-subset(grange6merplaindrome, strand == "-")

TransArxMotif<-distanceToNearest(grange6merPlus, grange6merPlainDromeMinus, ignore.strand=TRUE)
intergerforTrans<-subset(TransArxMotif, distance>0 & distance<= 200 )%>%countRnodeHits()
nonOverlappingMinusTransClusters<-subset(grange6merPlainDromeMinus, intergerforTrans)


cisArxMotif<-distanceToNearest(grange6merPlus)
cisARxMotifs<-subset(grange6mer, seqnames=="chrX")
TransArxMotifs<-subset(grange6mer, !seqnames=="chrX")

intergerforCis<-subset(cisArxMotif, distance>0 & distance<= 200 )%>%countRnodeHits()
nonOverlappingMinusCisMotifs<-subset(grange6merPlainDromeMinus, intergerforCis)



##Data tables of wheere these motifs are!
UCSCgenes <- import("~/Scripts/March/mm9.bed")
promoters <- promoters(UCSCgenes)
enhancerGrange <- import(con = "~/DataFiles/Enhancer Tracks/Mouse/Enhanceresmm9.bed")


CisTransDataTable <- rbind(
  cbind(
    length(cisARxMotifs),
    Arx6mer <- sum(countOverlaps(cisARxMotifs, UCSCgenes)),
    sum(countOverlaps(cisARxMotifs, promoters)),
    sum(countOverlaps(cisARxMotifs, enhancerGrange))
),
cbind(
    length(TransArxMotifs),
    sum(countOverlaps(TransArxMotifs, UCSCgenes)),
    sum(countOverlaps(TransArxMotifs, promoters)),
    sum(countOverlaps(TransArxMotifs, enhancerGrange))))
CisTransDataTable<-CisTransDataTable%>%as.data.frame()

rownames(CisTransDataTable)<- c("Cis Motifs", "Trans Motifs")
colnames(CisTransDataTable)<- c("total",
                                "motifs in genes",
                                "motifs in Promoters",
                                "motifs in enhancers")
ggplotCisTransDataTable<-rownames_to_column(CisTransDataTable)

ggplot(ggplotCisTransDataTable, mapping= aes(x=rowname, y=total, fill=rowname))+
  geom_bar(stat="identity")+
  xlab(label= "Cis Or Trans Orientation")+
  ylab(label= " Number Of Arx Motifs")+
  guides(fill=guide_legend(title=NULL))+
  theme_bw()


##generating plots for distances

ggplotCisPromoters<-distanceToNearest(cisARxMotifs, promoters)%>%as.data.frame()
ggplotTransPromoters<-distanceToNearest(TransArxMotifs, promoters)%>%as.data.frame()

ggplotCisTrans<-merge(ggplotCisPromoters[3], ggplotTransPromoters[3], by=0, all= TRUE)
ggplotCisTransReshaped<-reshape(ggplotCisTrans,
          varying = c("distance.x", "distance.y"),
          v.names = "Distance",
          timevar = "Orientation",
          times = c("Cis", "Trans"),
          direction = "long")

ggplot(ggplotCisTransReshaped, aes(x=Distance, group=Orientation, fill=Orientation))+
  geom_histogram(bins = 1000)+
  xlab(label = "Distance To The Closest Promoter(Base Pairs)")+
  ylab(label= "Number of Motifs")+
  theme(text = element_text(size=40))+
  scale_x_continuous(limits = c(0, 2000000))+
  guides(fill=guide_legend(title=NULL))+
  theme_bw()


averageDistanceFromPromoterCis<-sum(apply(ggplotCisPromoters[3], 2, as.numeric))/dim(ggplotCisPromoters)[1]
averageDistanceFromPromoterTrans<-sum(apply(ggplotTransPromoters[3], 2, as.numeric))/dim(ggplotTransPromoters)[1]

cisNonPromoters<-subset(ggplotCisPromoters, distance>0)
transNonPromoters<-subset(ggplotTransPromoters, distance>0)
averageDistanceFromPromoterCisOutsidePromoter<-sum(apply(cisNonPromoters[3], 2, as.numeric))/dim(cisNonPromoters)[1]
averageDistanceFromPromoterTransOutsidePromoter<-sum(apply(transNonPromoters[3], 2, as.numeric))/dim(transNonPromoters)[1]

rbind(averageDistanceFromPromoterCis,averageDistanceFromPromoterCisOutsidePromoter, averageDistanceFromPromoterTrans,averageDistanceFromPromoterTransOutsidePromoter
      )
