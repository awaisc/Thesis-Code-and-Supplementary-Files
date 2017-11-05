library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm9)
library(pander)
library(GenomicFeatures)
library(GenomicRanges)
arx6Mer <-
  rbind(
    A = c(0, 1, 1, 0, 0, 1),
    C = c(0, 0, 0, 0, 0, 0),
    G = c(0, 0, 0, 0, 0, 0) ,
    T = c(1, 0, 0, 1, 1, 0)
  )


gtfUCSCexonscoding<-import("~/Scripts/March/FullMm9genome.GTF")
gtfUCSCgenes<- import("~/Scripts/March/mm9.bed")
enhancerGrange <- import( con = "~/DataFiles/Enhancer Tracks/Mouse/Enhanceresmm9.bed")
methylationGrange<- import(con = "~/DataFiles/Methylation Tracks/Mouse/CpGIslands.bed")
promoterGrange<- promoters(gtfUCSCgenes)

arx6MerTFBS<-matchPWM(arx6Mer, BSgenome.Mmusculus.UCSC.mm9, "100%")

arx6MerTable<-cbind( length(arx6MerTFBS),
arx6MerInsideExons<-length(findOverlaps(gtfUCSCexonscoding,arx6MerTFBS)),
arx6MerInsideGenes<-length(findOverlaps(gtfUCSCgenes, arx6MerTFBS)),
arx6MersInsideEnhancers<- length(findOverlaps(enhancerGrange, arx6MerTFBS)),
arx6MerInsideMethylation<- length(findOverlaps(methylationGrange, arx6MerTFBS)),
arx6merInsidePromoter<- length(findOverlaps(promoterGrange, arx6MerTFBS))
)
arx6MerTable[1]-sum(arx6MerTable[2:6])
