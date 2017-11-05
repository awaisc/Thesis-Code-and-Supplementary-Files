## Checking to see if the numbers are robust

library(magrittr)
library(GenomicRanges)
library(ggplot2)
library(magrittr)
library(tibble)
library(pander)
library(reshape2)
library(plyr)
library(MotifDb)
library(magrittr)
library(reshape2)
library(BSgenome.Hsapiens.UCSC.hg19)

enhancerGrange <-
  import(con = "~/DataFiles/Enhancer Tracks/Mouse/mouse_permissive_enhancers_phase_1_and_2.bed")
UCSCgenes <- import("~/DataFiles/Gene Tracks/Mouse/mm9.bed")
startSites<-subset(import("~/DataFiles/Gene Tracks/Mouse/FullMm9genome.GTF"), type== "start_codon")
promoters <- promoters(UCSCgenes)
genome<-BSgenome.Mmusculus.UCSC.mm9


ArxPlaindrmicMinus1<-rbind( A=c(0,1,1,0,0,0,0,1,1,0), 
                            C=c(0,0,0,0,0,0,0),
                            G=c(0,0,0,0,0,0,0),
                            T=c(1,0,0,1,1,1,1,0,0,1))

arx6MerPWMNospace<-rbind( A=c(0,1,1,0,0,1,0,1,1,0,0,1), 
                              C=c(0,0,0,0,0,0,0),
                              G=c(0,0,0,0,0,0,0) ,
                              T=c(1,0,0,1,1,0,1,0,0,1,1,0))

arx6MerPWM1space<-rbind( A=c(0,1,1,0,0,1,0.25,1,0,0,1,1,0), 
                         C=c(0,0,0,0,0,0,0.25,0),
                         G=c(0,0,0,0,0,0,0.25,0),
                         T=c(1,0,0,1,1,0,0.25,0,1,1,0,0,1))

arx6MerPWM2space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,1,0,0,1,1,0), 
                         C=c(0,0,0,0,0,0,0.25,0.25),
                         G=c(0,0,0,0,0,0,0.25,0.25,0) ,
                         T=c(1,0,0,1,1,0,0.25,0.25,0,1,1,0,0,1))

arx6MerPWM3space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0.25,1,0,0,1,1,0),
                         C=c(0,0,0,0,0,0,0.25,0.25,0.25,0),
                         G=c(0,0,0,0,0,0,0.25,0.25,0.25,0),
                         T=c(1,0,0,1,1,0,0.25,0.25,0.25,0,1,1,0,0,1))

arx6MerPWM4space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0.25,0.25,1,0,0,1,1,0),
                         C=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25),
                         G=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25) ,
                         T=c(1,0,0,1,1,0,0.25,0.25,0.25,0.25,0,1,1,0,0,1))

arx6MerPWM5space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0.25,0.25,0.25,1,0,0,1,1,0),
                         C=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0),
                         G=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0),
                         T=c(1,0,0,1,1,0,0.25,0.25,0.25,0.25,0.25,0,1,1,0,0,1))
                         
arx6MerPWM6space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0.25,0.25,0.25,0.25,1,0,0,1,1,0),
                         C=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0.25,0),
                         G=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0.25,0),
                         T=c(1,0,0,1,1,0,0.25,0.25,0.25,0.25,0.25,0.25,0,1,1,0,0,1))

arx6MerPWM7space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0.25,0.25,0.25,0.25,0.25,1,0,0,1,1,0),
                         C=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0),
                         G=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0) ,
                         T=c(1,0,0,1,1,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,1,1,0,0,1))

### Tandeom Sites
arxtandemMinus1<-rbind(A=c(0,1,1,0,0,1,1,0,0,1),
                       C=c(0,0,0,0,0,0,0,0,0,0),
                       G=c(0,0,0,0,0,0,0,0,0,0),
                       T=c(1,0,0,1,1,0,0,1,1,0))
arxJolma<-rbind( A=c(0,1,1,0,0,0.25,1,1,0,0,1), 
                 C=c(0,0,0,0,0,0.25,0,0,0,0,0),
                 G=c(0,0,0,0,0,0.25,0,0,0,0,0),
                 T=c(1,0,0,1,1,0.25,0,0,1,1,0))
arxTandemNoSpace<-rbind( A=c(0,1,1,0,0,1,0,1,1,0,0,1),
                         C=c(0,0,0,0,0,0,0,0,0,0,0,0),
                         G=c(0,0,0,0,0,0,0,0,0,0,0,0) ,
                         T=c(1,0,0,1,1,0,1,0,0,1,1,0))

arxTandem1Space<-rbind( A=c(0,1,1,0,0,1,0.25,0,1,1,0,0,1),
                        C=c(0,0,0,0,0,0,0.25,0),
                        G=c(0,0,0,0,0,0,0.25,0) ,
                        T=c(1,0,0,1,1,0,0.25,1,0,0,1,1,0))

arxTandem2Space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0,1,1,0,0,1),
                        C=c(0,0,0,0,0,0,0.25,0.25,0),
                        G=c(0,0,0,0,0,0,0.25,0.25,0) ,
                        T=c(1,0,0,1,1,0,0.25,0.25,1,0,0,1,1,0))

arxTandem3Space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0.25,0,1,1,0,0,1),
                        C=c(0,0,0,0,0,0,0.25,0.25,0.25,0),
                        G=c(0,0,0,0,0,0,0.25,0.25,0.25,0) ,
                        T=c(1,0,0,1,1,0,0.25,0.25,0.25,1,0,0,1,1,0))
arxTandem4Space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0.25,0.25,0,1,1,0,0,1),
                        C=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0),
                        G=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0) ,
                        T=c(1,0,0,1,1,0,0.25,0.25,0.25,0.25,1,0,0,1,1,0))

arxTandem5Space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0.25,0.25,0.25,0,1,1,0,0,1), 
                        C=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0),
                        G=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0) ,
                        T=c(1,0,0,1,1,0,0.25,0.25,0.25,0.25,0.25,1,0,0,1,1,0))
arxTandem6Space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0.25,0.25,0.25,0.25,0,1,1,0,0,1), 
                        C=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0.25,0),
                        G=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0.25,0) ,
                        T=c(1,0,0,1,1,0,0.25,0.25,0.25,0.25,0.25,0.25,1,0,0,1,1,0))
arxTandem7Space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,1,1,0,0,1), 
                        C=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0),
                        G=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0) ,
                        T=c(1,0,0,1,1,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,1,0,0,1,1,0))
arxTandem8Space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0.25,0.25,0,1,1,0,0,1), 
                        C=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0),
                        G=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0) ,
                        T=c(1,0,0,1,1,0,0.25,0.25,0.25,0.25,1,0,0,1,1,0))
##requires code from the 16-4-2017 to run
grangeJolmaMinus<-
  matchPWM(arxJolma, genome, "100%")
grangeplaindromicMinus1 <-
  matchPWM(ArxPlaindrmicMinus1, genome, "100%")
grangeplaindromicNospace <-
  matchPWM(arx6MerPWMNospace, genome, "100%")
grangeplaindromic1space <-
  matchPWM(arx6MerPWM1space, genome, "100%")
grangeplaindromic2space <-
  matchPWM(arx6MerPWM2space, genome, "100%")
grangeplaindromic3space <-
  matchPWM(arx6MerPWM3space, genome, "100%")
grangeplaindromic4space <-
  matchPWM(arx6MerPWM4space, genome, "100%")
grangeplaindromic5space <-
  matchPWM(arx6MerPWM5space, genome, "100%")
grangeplaindromic6space <-
  matchPWM(arx6MerPWM6space, genome, "100%")
grangeplaindromic7space <-
  matchPWM(arx6MerPWM7space, genome, "100%")

grangeTandemMinusOne <-
  matchPWM(arxtandemMinus1, genome, "100%")
grangeTandemNoSpace<-
  matchPWM(arxTandemNoSpace, genome, "100%")
grangeTandem1space <-
  matchPWM(arxTandem1Space, genome, "100%")
grangeTandem2space <-
  matchPWM(arxTandem2Space, genome, "100%")
grangeTandem3space <-
  matchPWM(arxTandem3Space, genome, "100%")
grangeTandem4space <-
  matchPWM(arxTandem4Space, genome, "100%")
grangeTandem5space <-
  matchPWM(arxTandem5Space, genome, "100%")
grangeTandem6space <-
  matchPWM(arxTandem6Space, genome, "100%")
grangeTandem7space <-
  matchPWM(arxTandem7Space, genome, "100%")

#grangeplaindromic1space<-matchPWM(arx6MerPWM1space, genome, "90%")
#grangeplaindromic2space<-matchPWM(arx6MerPWM2space, genome, "90%")
#grangeplaindromic3space<-matchPWM(arx6MerPWM3space, genome, "90%")
#grangeplaindromic4space<-matchPWM(arx6MerPWM4space, genome, "90%")
#grangeplaindromic5space<-matchPWM(arx6MerPWM4space, genome, "90%")
#grangeplaindromic6space<-matchPWM(arx6MerPWM4space, genome, "90%")
#grangeplaindromic7space<-matchPWM(arx6MerPWM4space, genome, "90%")
#grangeTandem1space<-matchPWM(arxTandem1Space, genome, "90%")
#grangeTandem2space<-matchPWM(arxTandem2Space, genome, "90%")
#grangeTandem3space<-matchPWM(arxTandem3Space, genome, "90%")
#grangeTandem4space<-matchPWM(arxTandem4Space, genome, "90%")
#grangeTandem5space<-matchPWM(arxTandem5Space, genome, "90%")
#grangeTandem6space<-matchPWM(arxTandem6Space, genome, "90%")
#grangeTandem7space<-matchPWM(arxTandem7Space, genome, "90%")
tandemDataTable <- rbind(
  cbind(
  length(grangeJolmaMinus),
    sum(countOverlaps(grangeJolmaMinus, UCSCgenes)),
    sum(countOverlaps(grangeJolmaMinus, promoters)),
    sum(countOverlaps(grangeJolmaMinus, enhancerGrange)),
  (length(grangeJolmaMinus)-sum(countOverlaps(grangeJolmaMinus, enhancerGrange))-
     sum(countOverlaps(grangeJolmaMinus, promoters))-  sum(countOverlaps(grangeJolmaMinus, UCSCgenes)))
),
  cbind(
    numberofTandem <- length(grangeTandemMinusOne),
    dataTableNoGenesminus1 <-
      sum(countOverlaps(grangeTandemMinusOne, UCSCgenes)),
    dataTableMinus1 <-
      sum(countOverlaps(grangeTandemMinusOne, promoters)),
    dataTableMinus1r <-
      sum(countOverlaps(grangeTandemMinusOne, enhancerGrange)),
    (length(grangeTandemMinusOne)-sum(countOverlaps(grangeTandemMinusOne, enhancerGrange))-
       sum(countOverlaps(grangeTandemMinusOne, promoters))-  sum(countOverlaps(grangeTandemMinusOne, UCSCgenes)))
  ),
  cbind(
    numberofTandemNoSpaceSites <- length(grangeTandemNoSpace),
    dataTableNoGenes <-
      sum(countOverlaps(grangeTandemNoSpace, UCSCgenes)),
    dataTableNoSpacePromoters <-
      sum(countOverlaps(grangeTandemNoSpace, promoters)),
    dataTableNoSpaceEnhancer <-
      sum(countOverlaps(grangeTandemNoSpace, enhancerGrange)),
    (length(grangeTandemNoSpace)-sum(countOverlaps(grangeTandemNoSpace, enhancerGrange))-
      sum(countOverlaps(grangeTandemNoSpace, promoters))-  sum(countOverlaps(grangeTandemNoSpace, UCSCgenes)))
 ),
 cbind(
    numberofTandem1spaceSites <- length(grangeTandem1space),
    dataTable1SpaceGenes <-
      sum(countOverlaps(grangeTandem1space, UCSCgenes)),
    dataTable1SpacePromoters <-
      sum(countOverlaps(grangeTandem1space, promoters)),
    dataTable1SpaceEnhancer <-
      sum(countOverlaps(grangeTandem1space, enhancerGrange)),
    (length(grangeTandem1space)-sum(countOverlaps(grangeTandem1space, enhancerGrange))-
       sum(countOverlaps(grangeTandem1space, promoters))-  sum(countOverlaps(grangeTandem1space, UCSCgenes)))
  ),
  cbind(
    numberofTandem2spaceSites <- length(grangeTandem2space),
    dataTable2SpaceGenes <-
      sum(countOverlaps(grangeTandem2space, UCSCgenes)),
    dataTable2SpacePromoters <-
      sum(countOverlaps(grangeTandem2space, promoters)),
    dataTable2SpaceEnhancer <-
      sum(countOverlaps(grangeTandem2space, enhancerGrange)),
    (length(grangeTandem2space)-sum(countOverlaps(grangeTandem2space, enhancerGrange))-
       sum(countOverlaps(grangeTandem2space, promoters))-  sum(countOverlaps(grangeTandem2space, UCSCgenes)))
  ),
  cbind(
    numberofTandem3spaceSites <- length(grangeTandem3space),
    dataTable3SpaceGenes <-
      sum(countOverlaps(grangeTandem3space, UCSCgenes)),
    dataTable3SpacePromoters <-
      sum(countOverlaps(grangeTandem3space, promoters)),
    dataTable3SpaceEnhancer <-
      sum(countOverlaps(grangeTandem3space, enhancerGrange)),
    (length(grangeTandem3space)-sum(countOverlaps(grangeTandem3space, enhancerGrange))-
       sum(countOverlaps(grangeTandem3space, promoters))-  sum(countOverlaps(grangeTandem3space, UCSCgenes)))
  ),
  cbind(
    numberofTandem4spaceSites <- length(grangeTandem4space),
    dataTable4SpaceGenes <-
      sum(countOverlaps(grangeTandem4space, UCSCgenes)),
    dataTable4SpacePromoters <-
      sum(countOverlaps(grangeTandem4space, promoters)),
    dataTable4SpaceEnhancer <-
      sum(countOverlaps(grangeTandem4space, enhancerGrange)),
    (length(grangeTandem4space)-sum(countOverlaps(grangeTandem4space, enhancerGrange))-
       sum(countOverlaps(grangeTandem4space, promoters))-  sum(countOverlaps(grangeTandem4space, UCSCgenes)))
  ), 
cbind(
    numberofTandem5spaceSites <- length(grangeTandem5space),
    dataTable5SpaceGenes <-
      sum(countOverlaps(grangeTandem5space, UCSCgenes)),
    dataTable5SpacePromoters <-
      sum(countOverlaps(grangeTandem5space, promoters)),
    dataTable5SpaceEnhancer <-
      sum(countOverlaps(grangeTandem5space, enhancerGrange)),
    (length(grangeTandem5space)-sum(countOverlaps(grangeTandem5space, enhancerGrange))-
       sum(countOverlaps(grangeTandem5space, promoters))-  sum(countOverlaps(grangeTandem5space, UCSCgenes)))
  ),
  cbind(
    numberofTandem6spaceSites <- length(grangeTandem6space),
    dataTable6SpaceGenes <-
      sum(countOverlaps(grangeTandem6space, UCSCgenes)),
    dataTable6SpacePromoters <-
      sum(countOverlaps(grangeTandem6space, promoters)),
    dataTable6SpaceEnhancer <-
      sum(countOverlaps(grangeTandem6space, enhancerGrange)),
    (length(grangeTandem6space)-sum(countOverlaps(grangeTandem6space, enhancerGrange))-
       sum(countOverlaps(grangeTandem6space, promoters))-  sum(countOverlaps(grangeTandem6space, UCSCgenes)))
  ),
  cbind(
    numberofTandem7spaceSites <- length(grangeTandem7space),
    dataTable7SpaceGenes <-
      sum(countOverlaps(grangeTandem7space, UCSCgenes)),
    dataTable7SpacePromoters <-
      sum(countOverlaps(grangeTandem7space, promoters)),
    dataTable7SpaceEnhancer <-
      sum(countOverlaps(grangeTandem7space, enhancerGrange)),
    (length(grangeTandem7space)-sum(countOverlaps(grangeTandem7space, enhancerGrange))-
       sum(countOverlaps(grangeTandem7space, promoters))-  sum(countOverlaps(grangeTandem7space, UCSCgenes)))
  )
) %>% as.data.frame

colnames(tandemDataTable) <- c("Total",
                               "Motifs in genes",
                               "Motifs in promoters",
                               "Motifs in enhancers",
                               "Non Coding")


rownames(tandemDataTable) <- c("Arx Jolma",
                               "Minus one",
                               "No Space",
                               "1 Space",
                               "2 Space",
                               "3 Space",
                               "4 Space",
                               "5 Space",
                               "6 Space",
                               "7 Space")
tandemDataTable %>% pander()

tandemDataTable <- rownames_to_column(tandemDataTable)
reshapedTandemDataTable<-reshape(tandemDataTable,
                                varying = c( "Motifs in promoters", "Motifs in enhancers", "Non Coding", "Motifs in genes"),
                                v.names = "Numbers of Motif",
                                timevar = "Location",
                                times = c( "Promoters", "Enhancers", "Non coding","Genes" ),
                                direction = "long")
ggplot(reshapedTandemDataTable, aes(x = rowname, y = `Numbers of Motif`, fill = `Location`)) +
  geom_bar(stat = "identity") +
  xlab(label= "Number of Nucleotides Between Motifs")+
  ylab(label= "Number of ARX Motifs")+
  guides(fill=guide_legend(title="Genomic Location"))+
  theme_bw()+
  theme(axis.title = element_text(size=30, face = "bold"),
        axis.text =element_text(size=20),
        legend.text = element_text(size=20), 
        legend.title = element_text(size=30))+
  scale_color_manual(values=c(`Enhancer`="#999999", `Genes`="#E69F00", `Non-coding`="#56B4E9", `Promoters`= "#56B4E9"))





planindromicDataTable <- rbind(
  cbind(
    length(grangeJolmaMinus),
    sum(countOverlaps(grangeJolmaMinus, UCSCgenes)),
    sum(countOverlaps(grangeJolmaMinus, promoters)),
    sum(countOverlaps(grangeJolmaMinus, enhancerGrange)),
    (length(grangeJolmaMinus)-sum(countOverlaps(grangeJolmaMinus, enhancerGrange))-
       sum(countOverlaps(grangeJolmaMinus, promoters))-  sum(countOverlaps(grangeJolmaMinus, UCSCgenes)))
  ),
  cbind(
    length(grangeplaindromicMinus1),
    sum(countOverlaps(grangeplaindromicMinus1, UCSCgenes)),
    sum(countOverlaps(grangeplaindromicMinus1, promoters)),
    sum(countOverlaps(grangeplaindromicMinus1, enhancerGrange)),
    (length(grangeplaindromicMinus1)-sum(countOverlaps(grangeplaindromicMinus1, enhancerGrange))-
       sum(countOverlaps(grangeplaindromicMinus1, promoters))-  sum(countOverlaps(grangeplaindromicMinus1, UCSCgenes)))
  ),
  cbind(
    length(grangeplaindromicNospace),
    sum(countOverlaps(grangeplaindromicNospace, UCSCgenes)),
    sum(countOverlaps(grangeplaindromicNospace, promoters)),
    sum(countOverlaps(grangeplaindromicNospace, enhancerGrange)),
    (length(grangeplaindromicNospace)-sum(countOverlaps(grangeplaindromicNospace, enhancerGrange))-
       sum(countOverlaps(grangeplaindromicNospace, promoters))-  sum(countOverlaps(grangeplaindromicNospace, UCSCgenes)))
  ),
  cbind(
    length(grangeplaindromic1space),
    Arx6mer <- sum(countOverlaps(grangeplaindromic1space, UCSCgenes)),
    sum(countOverlaps(grangeplaindromic1space, promoters)),
    sum(countOverlaps(grangeplaindromic1space, enhancerGrange)),
    (length(grangeplaindromic1space)-sum(countOverlaps(grangeplaindromic1space, enhancerGrange))-
       sum(countOverlaps(grangeplaindromic1space, promoters))-  sum(countOverlaps(grangeplaindromic1space, UCSCgenes)))
  ),
  cbind(
    length(grangeplaindromic2space),
    sum(countOverlaps(grangeplaindromic2space, UCSCgenes)),
    sum(countOverlaps(grangeplaindromic2space, promoters)),
    sum(countOverlaps(grangeplaindromic2space, enhancerGrange)),
    (length(grangeplaindromic2space)-sum(countOverlaps(grangeplaindromic2space, enhancerGrange))-
       sum(countOverlaps(grangeplaindromic2space, promoters))-  sum(countOverlaps(grangeplaindromic2space, UCSCgenes)))
  )
  ,
  cbind(
    numberOfArxSitesPlaindromic3Space <- length(grangeplaindromic3space),
    sum(countOverlaps(grangeplaindromic3space, UCSCgenes)),
    sum(countOverlaps(grangeplaindromic3space, promoters)),
    sum(countOverlaps(grangeplaindromic4space, enhancerGrange)),
    (length(grangeplaindromic3space)-sum(countOverlaps(grangeplaindromic3space, enhancerGrange))-
       sum(countOverlaps(grangeplaindromic3space, promoters))-  sum(countOverlaps(grangeplaindromic3space, UCSCgenes)))
  ),
  cbind(
    numberOfArxSitesPlaindromic4Space <- length(grangeplaindromic4space),
    sum(countOverlaps(grangeplaindromic4space, UCSCgenes)),
    sum(countOverlaps(grangeplaindromic4space, promoters)),
    sum(countOverlaps(grangeplaindromic4space, enhancerGrange)),
    (length(grangeplaindromic4space)-sum(countOverlaps(grangeplaindromic4space, enhancerGrange))-
       sum(countOverlaps(grangeplaindromic4space, promoters))-  sum(countOverlaps(grangeplaindromic4space, UCSCgenes)))
  ),
  cbind(
    numberOfArxSitesPlaindromic5Space <- length(grangeplaindromic5space),
    sum(countOverlaps(grangeplaindromic5space, UCSCgenes)),
    sum(countOverlaps(grangeplaindromic5space, promoters)),
    sum(countOverlaps(grangeplaindromic5space, enhancerGrange)),
    (length(grangeplaindromic5space)-sum(countOverlaps(grangeplaindromic5space, enhancerGrange))-
       sum(countOverlaps(grangeplaindromic5space, promoters))-  sum(countOverlaps(grangeplaindromic5space, UCSCgenes)))
  ),
  cbind(
    numberOfArxSitesPlaindromic6Space <- length(grangeplaindromic6space),
    sum(countOverlaps(grangeplaindromic6space, UCSCgenes)),
    sum(countOverlaps(grangeplaindromic6space, promoters)),
    sum(countOverlaps(grangeplaindromic6space, enhancerGrange)),
    (length(grangeplaindromic6space)-sum(countOverlaps(grangeplaindromic6space, enhancerGrange))-
       sum(countOverlaps(grangeplaindromic6space, promoters))-  sum(countOverlaps(grangeplaindromic6space, UCSCgenes)))
  ),
  cbind(
    numberOfArxSitesPlaindromic7Space <- length(grangeplaindromic7space),
    sum(countOverlaps(grangeplaindromic7space, UCSCgenes)),
    sum(countOverlaps(grangeplaindromic7space, promoters)),
    sum(countOverlaps(grangeplaindromic7space, enhancerGrange)),
    (length(grangeplaindromic7space)-sum(countOverlaps(grangeplaindromic7space, enhancerGrange))-
       sum(countOverlaps(grangeplaindromic7space, promoters))-  sum(countOverlaps(grangeplaindromic7space, UCSCgenes)))
  )
) %>% as.data.frame()
colnames(planindromicDataTable) <- c("Total",
                                     "Motifs in genes",
                                     "Motifs in Promoters",
                                     "Motifs in Enhancers",
                                     "Non Coding")
rownames(planindromicDataTable) <-c("Arx Jolma",
                                    "Minus one",
                                    "No Space",
                                    "1 Space",
                                    "2 Space",
                                    "3 Space",
                                    "4 Space",
                                    "5 Space",
                                    "6 Space",
                                    "7 Space")


planindromicDataTable %>% pander()
planindromicDataTable<- rownames_to_column(planindromicDataTable)

reshapedPlaindromicDataTable<-reshape(planindromicDataTable,
                                 varying = c( "Motifs in Promoters", "Motifs in Enhancers", "Non Coding", "Motifs in genes"),
                                 v.names = "Numbers of Motif",
                                 timevar = "Location",
                                 times = c( "Promoters", "Enhancers", "Non coding","Genes" ),
                                 direction = "long")
ggplot(reshapedPlaindromicDataTable, aes(x = rowname, y = `Numbers of Motif`, fill = `Location`)) +
  geom_bar(stat = "identity") +
  xlab(label= "Number of Nucleotides Between Motifs")+
  ylab(label= "Number of ARX Motifs")+
  guides(fill=guide_legend(title="Genomic Location"))+
  theme_bw()+
  theme(axis.title = element_text(size=30, face = "bold"),
        axis.text =element_text(size=20) )+
  scale_color_manual(values=c(`Enhancer`="#999999", `Genes`="#E69F00", `Non-coding`="#56B4E9", `Promoters`= "#56B4E9"))



### making histograms of distance of Tandem the Arx Start sites

dataFrameDistance1SpacePromoter <-
  distanceToNearest(grangeTandem1space, startSites) %>% 
  as.data.frame()
dataFrameDistance2SpacePromoter <-
  distanceToNearest(grangeTandem2space, startSites) %>%
  as.data.frame()
dataFrameDistance3SpacePromoter <-
  distanceToNearest(grangeTandem3space, startSites) %>%
  as.data.frame()
dataFrameDistance6SpacePromoter <-
  distanceToNearest(grangeTandem6space, startSites) %>%
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

ggplot(ggplotdataFrameDistanceofTandemicMotifsFromPromoter, aes(x=Distance, group=Space, color=Space))+
  geom_freqpoly(bins = 100)+
  xlab(label = "Distance To The Closest TSS(Base Pairs)")+
  ylab(label= "Number of Motifs")+
  theme_bw()+
  theme(axis.title = element_text(size=30, face = "bold"),
        axis.text =element_text(size=20) )+
  scale_x_continuous(limits = c(0, 100000))

ggplot(ggplotdataFrameDistanceofTandemicMotifsFromPromoter, aes(x=Distance, group=Space, color=Space))+
  geom_density(bins = 100)+
  xlab(label = "Distance To The Closest TSS(Base Pairs)")+
  ylab(label= "Number of Motifs")+
  theme_bw()+
  theme(axis.title = element_text(size=30, face = "bold"),
        axis.text =element_text(size=20) )+
  scale_x_continuous(limits = c(0, 100000))

##histogram of distances of Plaindromic Motifs


dataFrameDistancePlandromic1SpacePromoter  <-
  distanceToNearest(grangeplaindromic1space, startSites) %>% as.data.frame
dataFrameDistancePlandromic2SpacePromoter <-
  distanceToNearest(grangeplaindromic2space, startSites) %>% as.data.frame
dataFrameDistancePlandromic3SpacePromoter <-
  distanceToNearest(grangeplaindromic3space, startSites) %>% as.data.frame
dataFrameDistancePlandromic4SpacePromoter <-
  distanceToNearest(grangeplaindromic4space, startSites) %>% as.data.frame


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
ggplotdataFrameDistanceofPlaindromicMotifsFromPromoter$Distance<-ggplotdataFrameDistanceofPlaindromicMotifsFromPromoter$Distance%>%as.character%>%as.numeric()
ggplot(ggplotdataFrameDistanceofPlaindromicMotifsFromPromoter, aes(x=Distance, group=Space, color=Space))+
  geom_freqpoly(bins = 100)+
  xlab(label = "Distance To Closest TSS(Base Pairs)")+
  ylab(label= "Number Of Motifs")+
  scale_x_continuous(limits = c(0, 1000000))+
  theme_bw()+
  theme(axis.title = element_text(size=30, face = "bold"),
        axis.text =element_text(size=20))

ggplot(ggplotdataFrameDistanceofPlaindromicMotifsFromPromoter, aes(x=Distance, group=Space, color=Space))+
  geom_density(bins = 100)+
  xlab(label = "Distance To Closest TSS(Base Pairs)")+
  ylab(label= "Number Of Motifs")+
  scale_x_continuous(limits = c(0, 1000000))+
  theme_bw()+
  theme(axis.title = element_text(size=30, face = "bold"),
        axis.text =element_text(size=20))




## Average distances
Numeric<-apply(dataFrameDistanceofTandemMotifsFromPromoter, 2, as.numeric)
Numeric<-apply(dataFrameDistanceofPlandromicMotifsFromPromoter, 2, as.numeric)
Space1Av<-sum(na.omit(Numeric[,1]))/length(na.omit(Numeric[,1]))
Space2Av<-sum(na.omit(Numeric[,2]))/length(na.omit(Numeric[,2]))
Space3Av<-sum(na.omit(Numeric[,3]))/length(na.omit(Numeric[,3]))
Space4Av<-sum(na.omit(Numeric[,4]))/length(na.omit(Numeric[,4]))

