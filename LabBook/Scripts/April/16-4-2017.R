
## Generation for PWM for Planidromic sequences
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
library(magrittr)

enhancerGrange <-
  import(con = "~/DataFiles/Enhancer Tracks/Mouse/Enhanceresmm9.bed")
UCSCgenes <- import("~/Scripts/March/mm9.bed")
promoters <- promoters(UCSCgenes)

arx6MerPWMnospace <- MotifDb::query(MotifDb, "arx")[[6]]
arx6MerPWM1space <-
  rbind(
    A = c(0, 1, 1, 0, 0, 1, 0.25, 1, 0, 0, 1, 1, 0),
    C = c(0, 0, 0, 0, 0, 0, 0.25, 0),
    G = c(0, 0, 0, 0, 0, 0, 0.25, 0) ,
    T = c(1, 0, 0, 1, 1, 0, 0.25, 0, 1, 1, 0, 0, 1)
  )

arx6MerPWM2space <-
  rbind(
    A = c(0, 1, 1, 0, 0, 1, 0.25, 0.25, 1, 0, 0, 1, 1, 0),
    C = c(0, 0, 0, 0, 0, 0, 0.25, 0.25),
    G = c(0, 0, 0, 0, 0, 0, 0.25, 0.25, 0) ,
    T = c(1, 0, 0, 1, 1, 0, 0.25, 0.25, 0, 1, 1, 0, 0, 1)
  )

arx6MerPWM3space <-
  rbind(
    A = c(0, 1, 1, 0, 0, 1, 0.25, 0.25, 0.25, 1, 0, 0, 1, 1, 0),
    C = c(0, 0, 0, 0, 0, 0, 0.25, 0.25, 0.25),
    G = c(0, 0, 0, 0, 0, 0, 0.25, 0.25, 0.25, 0) ,
    T = c(1, 0, 0, 1, 1, 0, 0.25, 0.25, 0.25, 0, 1, 1, 0, 0, 1)
  )

arx6MerPWM4space <-
  rbind(
    A = c(0, 1, 1, 0, 0, 1, 0.25, 0.25, 0.25, 0.25, 1, 0, 0, 1, 1, 0),
    C = c(0, 0, 0, 0, 0, 0, 0.25, 0.25, 0.25, 0.25),
    G = c(0, 0, 0, 0, 0, 0, 0.25, 0.25, 0.25, 0.25) ,
    T = c(1, 0, 0, 1, 1, 0, 0.25, 0.25, 0.25, 0.25, 0, 1, 1, 0, 0, 1)
  )

grangeplaindromic1space <-
  matchPWM(arx6MerPWM1space, BSgenome.Mmusculus.UCSC.mm9, "90%")
grangeplaindromic2space <-
  matchPWM(arx6MerPWM2space, BSgenome.Mmusculus.UCSC.mm9, "90%")
grangeplaindromic3space <-
  matchPWM(arx6MerPWM3space, BSgenome.Mmusculus.UCSC.mm9, "90%")
grangeplaindromic4space <-
  matchPWM(arx6MerPWM4space, BSgenome.Mmusculus.UCSC.mm9, "90%")

##Databale results

planindromicDataTable <- rbind(
  cbind(
    length(grangeplaindromic1space),
    Arx6mer <- sum(countOverlaps(grangeplaindromic1space, UCSCgenes)),
    sum(countOverlaps(grangeplaindromic1space, promoters)),
    sum(countOverlaps(grangeplaindromic1space, enhancerGrange))
  ),
  
  cbind(
    length(grangeplaindromic2space),
    sum(countOverlaps(grangeplaindromic2space, UCSCgenes)),
    sum(countOverlaps(grangeplaindromic2space, promoters)),
    sum(countOverlaps(grangeplaindromic2space, enhancerGrange))
  )
  ,
  cbind(
    numberOfArxSitesPlaindromic3Space <- length(grangeplaindromic3space),
    sum(countOverlaps(grangeplaindromic3space, UCSCgenes)),
    sum(countOverlaps(grangeplaindromic3space, promoters)),
    sum(countOverlaps(grangeplaindromic4space, enhancerGrange))
  ),
  
  cbind(
    numberOfArxSitesPlaindromic4Space <- length(grangeplaindromic4space),
    sum(countOverlaps(grangeplaindromic4space, UCSCgenes)),
    sum(countOverlaps(grangeplaindromic4space, promoters)),
    sum(countOverlaps(grangeplaindromic4space, enhancerGrange))
  )
) %>% as.data.frame()

colnames(planindromicDataTable) <- c("Total",
                                     "Motifs in genes",
                                     "Motifs in Promoters",
                                     "Motifs in Enhancers")
rownames(planindromicDataTable) <- c("1 Space",
                                     "2 Space",
                                     "3 Space",
                                     "4 Space")
planindromicDataTable %>% pander()
planindromicDataTable<-rownames_to_column(planindromicDataTable)
ggplot(planindromicDataTable, aes(x = rowname, y = Total, fill= rowname)) +
  geom_bar(stat = "identity") +
  xlab(label="Number of Nucleotides Between Motifs")+
  ylab(label= "Number of Arx Motifs")+
  guides(fill=guide_legend(title=NULL))+
  theme(axis.title.x = element_text(face="bold", colour="#990000", size=20),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=16))+
  theme_bw()


### Tandeom Sites

arxTandem1Space <-
  rbind(
    A = c(0, 1, 1, 0, 0, 1, 0.25, 0, 1, 1, 0, 0, 1),
    C = c(0, 0, 0, 0, 0, 0, 0.25, 0),
    G = c(0, 0, 0, 0, 0, 0, 0.25, 0) ,
    T = c(1, 0, 0, 1, 1, 0, 0.25, 1, 0, 0, 1, 1, 0)
  )
arxTandem2Space <-
  rbind(
    A = c(0, 1, 1, 0, 0, 1, 0.25, 0.25, 0, 1, 1, 0, 0, 1),
    C = c(0, 0, 0, 0, 0, 0, 0.25, 0.25, 0),
    G = c(0, 0, 0, 0, 0, 0, 0.25, 0.25, 0) ,
    T = c(1, 0, 0, 1, 1, 0, 0.25, 0.25, 1, 0, 0, 1, 1, 0)
  )

arxTandem3Space <-
  rbind(
    A = c(0, 1, 1, 0, 0, 1, 0.25, 0.25, 0.25, 0, 1, 1, 0, 0, 1),
    C = c(0, 0, 0, 0, 0, 0, 0.25, 0.25, 0.25, 0),
    G = c(0, 0, 0, 0, 0, 0, 0.25, 0.25, 0.25, 0) ,
    T = c(1, 0, 0, 1, 1, 0, 0.25, 0.25, 0.25, 1, 0, 0, 1, 1, 0)
  )

arxTandem4Space <-
  rbind(
    A = c(0, 1, 1, 0, 0, 1, 0.25, 0.25, 0.25, 0.25, 0, 1, 1, 0, 0, 1),
    C = c(0, 0, 0, 0, 0, 0, 0.25, 0.25, 0.25, 0.25, 0),
    G = c(0, 0, 0, 0, 0, 0, 0.25, 0.25, 0.25, 0.25, 0) ,
    T = c(1, 0, 0, 1, 1, 0, 0.25, 0.25, 0.25, 0.25, 1, 0, 0, 1, 1, 0)
  )

grangeTandem1space <-
  matchPWM(arxTandem1Space, BSgenome.Mmusculus.UCSC.mm9, "90%")
grangeTandem2space <-
  matchPWM(arxTandem2Space, BSgenome.Mmusculus.UCSC.mm9, "90%")
grangeTandem3space <-
  matchPWM(arxTandem3Space, BSgenome.Mmusculus.UCSC.mm9, "90%")
grangeTandem4space <-
  matchPWM(arxTandem4Space, BSgenome.Mmusculus.UCSC.mm9, "90%")

##Tandem DataTable
tandemDataTable <- rbind(
  cbind(
    numberofTandem1spaceSites <- length(grangeTandem1space),
    dataTable1SpaceGenes <-
      sum(countOverlaps(grangeTandem1space, UCSCgenes)),
    dataTable1SpacePromoters <-
      sum(countOverlaps(grangeTandem1space, promoters)),
    dataTable1SpaceEnhancer <-
      sum(countOverlaps(grangeTandem1space, enhancerGrange))
  ),
  cbind(
    numberofTandem2spaceSites <- length(grangeTandem2space),
    dataTable2SpaceGenes <-
      sum(countOverlaps(grangeTandem2space, UCSCgenes)),
    dataTable2SpacePromoters <-
      sum(countOverlaps(grangeTandem2space, promoters)),
    dataTable2SpaceEnhancer <-
      sum(countOverlaps(grangeTandem2space, enhancerGrange))
  ),
  cbind(
    numberofTandem3spaceSites <- length(grangeTandem3space),
    dataTable3SpaceGenes <-
      sum(countOverlaps(grangeTandem3space, UCSCgenes)),
    dataTable3SpacePromoters <-
      sum(countOverlaps(grangeTandem3space, promoters)),
    dataTable3SpaceEnhancer <-
      sum(countOverlaps(grangeTandem3space, enhancerGrange))
  ),
  cbind(
    numberofTandem4spaceSites <- length(grangeTandem4space),
    dataTable4SpaceGenes <-
      sum(countOverlaps(grangeTandem4space, UCSCgenes)),
    dataTable4SpacePromoters <-
      sum(countOverlaps(grangeTandem4space, promoters)),
    dataTable4SpaceEnhancer <-
      sum(countOverlaps(grangeTandem4space, enhancerGrange))
  )
) %>% as.data.frame

colnames(tandemDataTable) <- c("Total",
                               "Motifs in genes",
                               "Motifs in promoters",
                               "Motifs in enhancers")
rownames(tandemDataTable) <- c("1 Space",
                               "2 Space",
                               "3 Space",
                               "4 Space")
tandemDataTable %>% pander()

tandemDataTable <- rownames_to_column(tandemDataTable)

ggplot(tandemDataTable, aes(x = rowname, y = Total, fill= rowname)) +
  geom_bar(stat="identity")+
  xlab(label="Number of Nucleotides Between Motifs")+
  ylab(label= "Number of Arx Motifs")+
  guides(fill=guide_legend(title=NULL))+
  theme_bw()


### making histograms of distance of Tandem the Arx Start sites

dataFrameDistance1SpacePromoter <-
  distanceToNearest(grangeTandem1space, promoters) %>% 
  as.data.frame()
dataFrameDistance2SpacePromoter <-
  distanceToNearest(grangeTandem2space, promoters) %>%
  as.data.frame()
dataFrameDistance3SpacePromoter <-
  distanceToNearest(grangeTandem3space, promoters) %>%
  as.data.frame()
dataFrameDistance4SpacePromoter <-
  distanceToNearest(grangeTandem4space, promoters) %>%
  as.data.frame()


dataFrameMerger<-function(z,x,c,v,b,n,m){
  
test<-merge(z[3],x[3],by=0, all=TRUE, row.names=NULL)
test2<-merge(test, c[3], by=0, all=TRUE, row.names=NULL)
test3<- merge(test2, v[3], by=0,all=TRUE, row.names=NULL)
test4<-merge(test3, b[3], by=0,all=TRUE, row.names=NULL)
test5<-merge(test4, n[3], by=0,all=TRUE, row.names=NULL)
test6<-merge(test5, m[3], by=0,all=TRUE, row.names=NULL)


return(test6)
}

dataFrameDistanceofTandemMotifsFromPromoter<-dataFrameMerger(dataFrameDistance1SpacePromoter,
                      dataFrameDistance2SpacePromoter, 
                      dataFrameDistance3SpacePromoter, 
                      dataFrameDistance4SpacePromoter)
dataFrameDistanceofTandemMotifsFromPromoter<- dataFrameDistanceofTandemMotifsFromPromoter[4:7]
colnames(dataFrameDistanceofTandemMotifsFromPromoter)<- c("1 Space",
                                                              "2 Space",
                                                              "3 Space",
                                                              "4 Space")
ggplotdataFrameDistanceofTandemicMotifsFromPromoter<-reshape(dataFrameDistanceofTandemMotifsFromPromoter,
                                                                varying = c("1 Space", "2 Space", "3 Space", "4 Space"),
                                                                v.names = "Distance",
                                                                timevar = "Space",
                                                                times = c("1 Nucleotide", "2 Nucleotide", "3 Nucleotide", "4 Nucleotide"),
                                                                direction = "long")

ggplot(ggplotdataFrameDistanceofTandemicMotifsFromPromoter, aes(x=Distance, group=Space, fill=Space))+
  geom_histogram(bins = 1000)+
  xlab(label = "Distance To The Closest Promoter(Base Pairs)")+
  ylab(label= "Number of Motifs")+
  theme(text = element_text(size=40))+
  scale_x_continuous(limits = c(0, 2000000))+
  theme_bw()



##histogram of distances of Plaindromic Motifs


dataFrameDistancePlandromic1SpacePromoter  <-
  distanceToNearest(grangeplaindromic1space, promoters) %>% as.data.frame
dataFrameDistancePlandromic2SpacePromoter <-
  distanceToNearest(grangeplaindromic2space, promoters) %>% as.data.frame
dataFrameDistancePlandromic3SpacePromoter <-
  distanceToNearest(grangeplaindromic3space, promoters) %>% as.data.frame
dataFrameDistancePlandromic4SpacePromoter <-
  distanceToNearest(grangeplaindromic4space, promoters) %>% as.data.frame
dataFrameDistancePlandromic5SpacePromoter <-
  distanceToNearest(grangeplaindromic5space, promoters) %>% as.data.frame
dataFrameDistancePlandromic6SpacePromoter <-
  distanceToNearest(grangeplaindromic6space, promoters) %>% as.data.frame
dataFrameDistancePlandromic7SpacePromoter <-
  distanceToNearest(grangeplaindromic7space, promoters) %>% as.data.frame

dataFrameDistanceofPlandromicMotifsFromPromoter <- dataFrameMerger(dataFrameDistancePlandromic1SpacePromoter,
                                                                   dataFrameDistancePlandromic2SpacePromoter,
                                                                   dataFrameDistancePlandromic3SpacePromoter,
                                                                   dataFrameDistancePlandromic4SpacePromoter,
                                                                   dataFrameDistancePlandromic5SpacePromoter,
                                                                   dataFrameDistancePlandromic6SpacePromoter,
                                                                   dataFrameDistancePlandromic7SpacePromoter)


dataFrameDistanceofPlandromicMotifsFromPromoter<-dataFrameDistanceofPlandromicMotifsFromPromoter[4:10]
colnames(dataFrameDistanceofPlandromicMotifsFromPromoter)<- c("1 Space",
                                                              "2 Space",
                                                              "3 Space",
                                                              "4 Space",
                                                              "5 Space",
                                                              "6 Space", 
                                                              "7 Space")
ggplotdataFrameDistanceofPlaindromicMotifsFromPromoter<-reshape(dataFrameDistanceofPlandromicMotifsFromPromoter,
                                                                varying = c("1 Space", "2 Space", "3 Space", "4 Space", "5 Space",
                                                                            "6 Space", 
                                                                            "7 Space"),
                                                                v.names = "Distance",
                                                                timevar = "Space",
                                                                times = c("1 Nucleotide", "2 Nucleotide", "3 Nucleotide", "4 Nucleotide", "5 Nucleotide", "6 Nucleotide", "7 Nucleotide"),
                                                                direction = "long")
ggplotdataFrameDistanceofPlaindromicMotifsFromPromoter$Distance<-ggplotdataFrameDistanceofPlaindromicMotifsFromPromoter$Distance%>%as.character%>%as.numeric()
ggplot(ggplotdataFrameDistanceofPlaindromicMotifsFromPromoter, aes(x=Distance, group=Space, fill=Space))+
  geom_histogram(bins = 1000)+
  xlab(label = "Distance To Closest Promoter(Base Pairs)")+
  ylab(label= "Number Of Motifs")+
  scale_x_continuous(limits = c(0, 2000000))+
  theme_bw()




## Average distances
Numeric<-apply(dataFrameDistanceofTandemMotifsFromPromoter, 2, as.numeric)
Numeric<-apply(dataFrameDistanceofPlandromicMotifsFromPromoter, 2, as.numeric)
Space1Av<-sum(na.omit(Numeric[,1]))/length(na.omit(Numeric[,1]))
Space2Av<-sum(na.omit(Numeric[,2]))/length(na.omit(Numeric[,2]))
Space3Av<-sum(na.omit(Numeric[,3]))/length(na.omit(Numeric[,3]))
Space4Av<-sum(na.omit(Numeric[,4]))/length(na.omit(Numeric[,4]))



