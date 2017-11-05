library(Gviz)
library(rtracklayer)
#UCSC Data Tracks

from <- 28000000
to   <- 29000000
chr <- "chrX"
##human Tracks

setwd("~/Research Proposal/")
#UCSC Track GCcontent
gcContent <- UcscTrack(genome = "hg19",
                         chromosome = chr,
                          track = "GC Percent",
                         table = "gc5Base",
                         trackType = "DataTrack",
                         start = "start",
                          end = "end", 
                         data = "score", 
                         type = "hist", 
                         window = -1,
                          windowSize = 1500,
                         fill.histogram = "black",
                          col.histogram = "black",
                         ylim = c(30, 70),
                         name = "GC Percent",
                       from = from,
                       to = to)

## gene Track

knownGenes <- UcscTrack(genome = "hg19", 
                        chromosome = chr,
                        track = "knownGene", 
                        trackType = "GeneRegionTrack",
                        rstarts = "exonStarts",
                        rends = "exonEnds",
                        gene = "name",
                        symbol = "name",
                        transcript = "name",
                        strand = "strand",
                        fill = "#8282d2",
                        name = "UCSC Genes",
                        stacking = "dense",
                        showID = TRUE)



ChIpconfirmedbed<- import.bed("5col.bed")
chipConfirmedSitesV3<-DataTrack(ChIpconfirmedbed,
                 genome = "hg19",
                 chromosome = chr,
                 to = to, 
                 from = from, 
                 name = "ChIPSeqV3Bed",
                 names = "name",
                 type = "histogram",
                 fill = "green",
                 color = "green")


##Ideogram Track
chromosomeBandTrack<- IdeogramTrack(genome = "hg19", 
              chromosome = chr)



##MYmade GC content Track
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(magrittr)
library(biomaRt)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(zoo)
library(magrittr)
library(parallel)
library(biomaRt)
library(MotifDb)
library(seqLogo)

mdb <- MotifDb
matrices.human <- MotifDb::query(mdb, 'hsapiens')

pfm.arx_jolma2013.jaspar <- MotifDb::query(mdb, 'Hsapiens-jolma2013-ARX')[[1]]
pcm.arx_jolma2013.jaspar <- round(100 * pfm.arx_jolma2013.jaspar)
seqLogo(pfm.arx_jolma2013.jaspar)
arxMotifs <- matchPWM(pcm.arx_jolma2013.jaspar,BSgenome.Hsapiens.UCSC.hg19[["chrX"]], "90%")
ranges<-cbind.data.frame(as.data.frame(ranges(arxMotifs)), as.data.frame(arxMotifs), seqnames = c(chr))

ARXMotifsGRangeForChrX<- makeGRangesFromDataFrame(ranges,
                                                  keep.extra.columns=TRUE,
                                                  ignore.strand=FALSE,
                                                  seqinfo=NULL,
                                                  seqnames.field="seqnames",
                                                  start.field="start",
                                                  end.field=c("end", "stop"),
                                                  strand.field="strand",
                                                  starts.in.df.are.0based=FALSE)

arxMotifsTrack <- AnnotationTrack(ARXMotifsGRangeForChrX,
                            genome= "hg19",
                            chromosome = chr,
                            name= "ARX binding sites",
                            stacking = "dense",
                            colour = "red")


### DNase Sec Track

dnaseSeq <- import("ENCFF752YMC.bigWig", which =
                     GRanges("chrX", IRanges(1, end = .Machine$integer.max - 1)))
dnaseSeqTrack <- DataTrack(dnaseSeq,
                           from = from,
                           to = to,
                           colour = "green",
                           name = "DnaseSeq105brainmale")


## Enhancer Tracks

enhancerTrack <- import("human_permissive_enhancers_phase_1_and_2.bed") %>% DataTrack( to = to, from = from, chromosome = chr, type = "boxplot")

plotTracks(enhancerTrack)


###plotTracks

plotTracks(list(chromosomeBandTrack, knownGenes, arxMotifsTrack, chipConfirmedSitesV3 ,dnaseSeqTrack, enhancerTrack), to= to, from = from )
