##Productive R DAY betches


library(Biostrings)
library(seqLogo)
library(GenomicRanges)
library(pander)
library(rtracklayer)
library(GenomicFeatures)
library(MotifDb)
library(BSgenome.Mmusculus.UCSC.mm9)
library(magrittr)
library(ggplot2)

gtfUCSCexonscoding<-import("~/Shiny App tutorials/Tut1/March/FullMm9genome.GTF")
gtfUCSCgenes<- import("~/Shiny App tutorials/Tut1/March/mm9.bed")
enhancerGrange <- import( con = "~/DataFiles/Enhancer Tracks/Mouse/Enhanceresmm9.bed")
methylationGrange<- import(con = "~/DataFiles/Methylation Tracks/Mouse/CpGIslands.bed")
promoterGrange<- promoters(gtfUCSC)

## Getting the other motifs
arx6Mer<-MotifDb::query(MotifDb, "Arx")[[6]]
arx6Mer <- round(Arx6mer*100)
arx6merTFBS<- matchPWM(arx6Mer, BSgenome.Mmusculus.UCSC.mm9, "90%", with.score = TRUE)

arx4mer<-rbind(A=c(0, 1, 1,0), C=c(0,0,0,0), G=c(0,0,0,0),
     T=c(1,0,0,1))
arx4merTFBS<-matchPWM(arx4mer, BSgenome.Mmusculus.UCSC.mm9, "100%")

arxJolma<-MotifDb::query(MotifDb, "Mmusculus-jolma2013-Arx")[[1]]
arxJolma <- round(arxJolma*100)
arxJolmaTFBS<- matchPWM(arxJolma, BSgenome.Mmusculus.UCSC.mm9, "90%", with.score = TRUE)




## make me some pie charts
arx6MerInsideExons<-length(findOverlaps(gtfUCSCexonscoding,arx6merTFBS))
arx6MerInsideGenes<-length(findOverlaps(gtfUCSCgenes, arx6merTFBS))
arx6MersInsideEnhancers<- length(findOverlaps(enhancerGrange, arx6merTFBS))
arx6MerInsideMethylation<- length(findOverlaps(methylationGrange, arx6merTFBS))
arx6merInsidePromoter<- length(findOverlaps(promoterGrange, arx6merTFBS))

arxJolmaInsideExon<- length(findOverlaps(gtfUCSCexonscoding, arxJolmaTFBS))
arxJolmaInsideGenes<- length(findOverlaps(gtfUCSCgenes, arxJolmaTFBS))
arxJolmaInsideEnhancer<- length(findOverlaps(enhancerGrange, arxJolmaTFBS))
arxJolmaInsideMethylation<- length(findOverlaps(methylationGrange, arxJolmaTFBS))
arxJolmaInsidePromoter<- length(findOverlaps(promoterGrange, arxJolmaTFBS))

arx4merInsideExon<- length(findOverlaps(gtfUCSCexonscoding, arx4merTFBS))
arx4merInsideGenes<- length(findOverlaps(gtfUCSCgenes, arx4merTFBS))
arx4merInsideEnhancer<- length(findOverlaps(enhancerGrange, arx4merTFBS))
arx4merInsideMethylation<- length(findOverlaps(methylationGrange, arx4merTFBS))
arx4merInsidePromoter<- length(findOverlaps(promoterGrange, arx4merTFBS))

##CHART TIME betches
table6mer<- cbind("6 mer Motifs in Total mm9 Genome" = length(arx6merTFBS),
              "6 mer Motifs Inside Genes"= arx6MerInsideGenes,
              "6 mer Motifs Inside Exons"=arx6MerInsideExons,
              "6 mer Motifs Inside Enhancers"=arx6MersInsideEnhancers,
              "6 mer Motifs Inside CpG Islands"= arx6MerInsideMethylation,
              "6 mer Motif Inside Promoter Site" =arx6merInsidePromoter)%>%pander()

tableJolma<- cbind("Jolma Motifs In Total mm9"= length(arxJolmaTFBS),
                   "Jolma Motifs Inside Genes"=arxJolmaInsideGenes ,
                   "Jolma Motifs Inside Exons"= arxJolmaInsideExon,
                   "Jolma Motifs Inside Enhancers" = arxJolmaInsideEnhancer,
                   "Jolma Motifs Inside CpG Islands" = arxJolmaInsideMethylation,
                   "Jolma Motif Inside Promoter Site" =arxJolmaInsidePromoter)%>%pander()

table<- cbind( cbind("4 mer Motifs in Total mm9 Genome" = length(arx6merTFBS),
                     "4 mer Motifs Inside Genes"= arx4merInsideGenes,
                     "4 mer Motifs Inside Exons"=arx4merInsideExon,
                     "4 mer Motifs Inside Enhancers"=arx4merInsideEnhancer,
                     "4 mer Motifs Inside CpG Islands"= arx4merInsideMethylation,
                     "4 mer Motif Inside Promoter Site" =arx4merInsidePromoter)%>%pander())
  
  

  
