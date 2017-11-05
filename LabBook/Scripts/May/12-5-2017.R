##taking older code and redoing it so that it presents number of motifs from the GENE start site. 


gtfUCSCexonscoding<-import("~/Scripts/March/FullMm9genome.GTF")
startsites<-subset(gtfUCSCexonscoding, type=="start_codon")

dataFrameDistance1SpacePromoter <-
  distanceToNearest(grangeTandem1space, startsites) %>% 
  as.data.frame()
dataFrameDistance2SpacePromoter <-
  distanceToNearest(grangeTandem2space, startsites) %>%
  as.data.frame()
dataFrameDistance3SpacePromoter <-
  distanceToNearest(grangeTandem3space, startsites) %>%
  as.data.frame()
dataFrameDistance6SpacePromoter <-
  distanceToNearest(grangeTandem6space, startsites) %>%
  as.data.frame()


dataFrameMerger<-function(z,x,c,v){
  
  test<-merge(z[3],x[3],by=0, all=TRUE, row.names=NULL)
  test2<-merge(test, c[3], by=0, all=TRUE, row.names=NULL)
  test3<- merge(test2, v[3], by=0,all=TRUE, row.names=NULL)
  return(test3)
}

dataFrameDistanceofTandemMotifsFromPromoter<-dataFrameMerger(dataFrameDistance1SpacePromoter,
                                                             dataFrameDistance2SpacePromoter, 
                                                             dataFrameDistance3SpacePromoter, 
                                                             dataFrameDistance6SpacePromoter)
dataFrameDistanceofTandemMotifsFromPromoter<- dataFrameDistanceofTandemMotifsFromPromoter[4:7]
colnames(dataFrameDistanceofTandemMotifsFromPromoter)<- c("1 Space",
                                                          "2 Space",
                                                          "3 Space",
                                                          "6 Space")
ggplotdataFrameDistanceofTandemicMotifsFromPromoter<-reshape(dataFrameDistanceofTandemMotifsFromPromoter,
                                                             varying = c("1 Space", "2 Space", "3 Space", "6 Space"),
                                                             v.names = "Distance",
                                                             timevar = "Space",
                                                             times = c("1 Nucleotide", "2 Nucleotide", "3 Nucleotide", "6 Nucleotide"),
                                                             direction = "long")

ggplot(ggplotdataFrameDistanceofTandemicMotifsFromPromoter, aes(x=Distance, group=Space, fill=Space))+
  geom_freqpoly(bins = 500, aes(colour=Space))+
  theme_bw()+
  xlab(label = "Distance To The Closest Transcription Start Site(Base Pairs)")+
  ylab(label= "Number of Motifs")+
  theme(text = element_text(size=12))+
  scale_x_continuous(limits = c(0, 200000))+
  scale_y_continuous(limits = c(0, 100))



##histogram of distances of Plaindromic Motifs


dataFrameDistancePlandromic1SpacePromoter  <-
  distanceToNearest(grangeplaindromic1space, startsites) %>% as.data.frame
dataFrameDistancePlandromic2SpacePromoter <-
  distanceToNearest(grangeplaindromic2space, startsites) %>% as.data.frame
dataFrameDistancePlandromic3SpacePromoter <-
  distanceToNearest(grangeplaindromic3space, startsites) %>% as.data.frame
dataFrameDistancePlandromic4SpacePromoter <-
  distanceToNearest(grangeplaindromic4space, startsites) %>% as.data.frame
head(dataFrameDistancePlandromic4SpacePromoter)
dataFrameDistanceofPlandromicMotifsFromPromoter <- dataFrameMerger(dataFrameDistancePlandromic1SpacePromoter,
                                                                   dataFrameDistancePlandromic2SpacePromoter,
                                                                   dataFrameDistancePlandromic3SpacePromoter,
                                                                   dataFrameDistancePlandromic4SpacePromoter)


dataFrameDistanceofPlandromicMotifsFromPromoter<-dataFrameDistanceofPlandromicMotifsFromPromoter[4:7]
colnames(dataFrameDistanceofPlandromicMotifsFromPromoter)<- c("1 Space",
                                                              "2 Space",
                                                              "3 Space",
                                                              "4 Space")
ggplotdataFrameDistanceofPlaindromicMotifsFromPromoter<-reshape(dataFrameDistanceofPlandromicMotifsFromPromoter,
                                                                varying = c("1 Space", "2 Space", "3 Space", "4 Space"),
                                                                v.names = "Distance",
                                                                timevar = "Space",
                                                                times = c("1 Nucleotide", "2 Nucleotide", "3 Nucleotide", "4 Nucleotide"),
                                                                direction = "long")
ggplot(ggplotdataFrameDistanceofPlaindromicMotifsFromPromoter, aes(x=Distance, group=Space))+
  geom_freqpoly(bins = 500, aes(colour=Space))+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(label = "Distance To Closest Transcription Start Site(Base Pairs)")+
  ylab(label= "Number Of Motifs")+
  scale_x_continuous(limits = c(0, 200000))+
  scale_y_continuous(limits = c(0, 100))




## Average distances
NumericTandem<-apply(dataFrameDistanceofTandemMotifsFromPromoter, 2, as.numeric)
NumericPlandrimoc<-apply(dataFrameDistanceofPlandromicMotifsFromPromoter, 2, as.numeric)
Space1Av<-sum(na.omit(NumericTandem[,1]))/length(na.omit(NumericTandem[,1]))
Space2Av<-sum(na.omit(NumericTandem[,2]))/length(na.omit(NumericTandem[,2]))
Space3Av<-sum(na.omit(NumericTandem[,3]))/length(na.omit(NumericTandem[,3]))
Space4Av<-sum(na.omit(NumericTandem[,4]))/length(na.omit(NumericTandem[,4]))


Space1AvPlandromic<-sum(na.omit(NumericPlandrimoc[,1]))/length(na.omit(NumericPlandrimoc[,1]))
Space2AvPlandromic<-sum(na.omit(NumericPlandrimoc[,2]))/length(na.omit(NumericPlandrimoc[,2]))
Space3AvPlandromic<-sum(na.omit(NumericPlandrimoc[,3]))/length(na.omit(NumericPlandrimoc[,3]))
Space4AvPlamdromic<-sum(na.omit(NumericPlandrimoc[,4]))/length(na.omit(NumericPlandrimoc[,4]))






