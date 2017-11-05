## code i wrote to examine density of ARX motifs
gtfUCSCexonscoding<-import("~/Scripts/March/FullMm9genome.GTF")
startsites<-subset(gtfUCSCexonscoding, type=="start_codon")


arxtandemMinus1<-rbind(A=c(0,1,1,0,0,1,1,0,0,1),
                       C=c(0,0,0,0,0,0,0,0,0,0),
                       G=c(0,0,0,0,0,0,0,0,0,0),
                       T=c(1,0,0,1,1,0,0,1,1,0))
arxJolma<-rbind( A=c(0,1,1,0,0,0.25,1,1,0,0,1), 
                 C=c(0,0,0,0,0,0.25,0,0,0,0,0),
                 G=c(0,0,0,0,0,0.25,0,0,0,0,0),
                 T=c(1,0,0,1,1,0.25,0,0,1,1,0))
ArxPlaindrmicMinus1<-rbind( A=c(0,1,1,0,0,0,0,1,1,0), 
                            C=c(0,0,0,0,0,0,0),
                            G=c(0,0,0,0,0,0,0),
                            T=c(1,0,0,1,1,1,1,0,0,1))







grangeJolmaMinus<-
  matchPWM(arxJolma, BSgenome.Mmusculus.UCSC.mm9, "100%")
grangeplaindromicMinus1 <-
  matchPWM(ArxPlaindrmicMinus1, BSgenome.Mmusculus.UCSC.mm9, "100%")
grangeTandemMinusOne <-
  matchPWM(arxtandemMinus1, BSgenome.Mmusculus.UCSC.mm9, "100%")



dataFrameDistanceJolma <-
  distanceToNearest(grangeJolmaMinus, startsites) %>% 
  as.data.frame()
dataFrameDistanceTandemMinusOne <-
  distanceToNearest(grangeTandemMinusOne, startsites) %>%
  as.data.frame()
dataFrameDistancePlaindromicMinusOne <-
  distanceToNearest(grangeplaindromicMinus1, startsites) %>%
  as.data.frame()

merge1<-merge(dataFrameDistanceJolma[3], dataFrameDistanceTandemMinusOne[3],by=0, all=TRUE)
dataFrameMinsOneJolmaDistanceFromPromoter<-merge(merge1, dataFrameDistancePlaindromicMinusOne[3],by=0, all=TRUE)

dataFrameMinsOneJolmaDistanceFromPromoter<-dataFrameMinsOneJolmaDistanceFromPromoter[3:5]
colnames(dataFrameMinsOneJolmaDistanceFromPromoter)<- c( "Jolma's",
                                                              "Tandem Minus 1",
                                                             "Plaindromic Minus 1")

ggplotdataFrameMinsOneJolmaDistanceFromPromoter<-reshape(dataFrameMinsOneJolmaDistanceFromPromoter,
                                                                varying = c("Plaindromic Minus 1",  "Tandem Minus 1", "Jolma's"),
                                                                v.names = "Distance",
                                                                timevar = "Model",
                                                                times = c("Plaindromic Minus 1",  "Tandem Minus 1", "Jolma's"),
                                                                direction = "long")
ggplot(ggplotdataFrameMinsOneJolmaDistanceFromPromoter, aes(x=Distance, group=Model))+
  geom_freqpoly(bins = 500, aes(colour=Model))+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(label = "Distance To Closest Transcription Start Site(Base Pairs)")+
  ylab(label= "Number Of Motifs")+
  scale_x_continuous(limits = c(0, 200000))+
  scale_y_continuous(limits = c(0, 250))
