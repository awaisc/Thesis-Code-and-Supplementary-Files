## Mouse Tracks because they're likely necessary as i'll be ultising mouse models
library(rtracklayer)
library(Gviz)
fromM<- 1
toM<-   1000000
chrM<- "chrY"

setwd("~/Research Proposal/")
##Ideogram

ideogramTrackMouse <- IdeogramTrack(genome = "mm9",
                                    chrM)
##Genes

knownGenesMouse <- UcscTrack(genome = "mm9", 
                        chromosome = chrM,
                        track = "knownGene", 
                        trackType = "GeneRegionTrack",
                        rstarts = "exonStarts",
                        rends = "exonEnds",
                        gene = "name",
                        symbol = "name",
                        transcript = "name",
                        strand = "strand",
                        fill = "Green",
                        name = "UCSC Genes",
                        stacking = "dense",
                        featureAnnotation = "id")

##ChIP confirmed TFBS
#ineedto do this and get it working

##GC content


## Mouse Arx Motifs

library(BSgenome.Mmusculus.UCSC.mm9)
library(Gviz)
library(magrittr)
mdb <- MotifDb
mouseArx <- MotifDb::query(mdb, 'Mmusculus-jolma2013-Arx')[[1]]
mouseArxPWM <- round(mouseArx*100)
mouseArxDataFrame<- matchPWM(mouseArxPWM, BSgenome.Mmusculus.UCSC.mm9[[chrM]], "90%")%>%
  ranges()%>%
  as.data.frame()%>% cbind(seqnames=(chrM))
grangesMouseArx<- makeGRangesFromDataFrame(mouseArxDataFrame,
                                           keep.extra.columns=FALSE,
                                           ignore.strand=FALSE,
                                           seqinfo=NULL,
                                           seqnames.field="seqnames",
                                           start.field="start",
                                           end.field=c("end", "stop"),
                                           strand.field="strand",
                                           starts.in.df.are.0based=FALSE)
MouseARXmotifs<-AnnotationTrack(grangesMouseArx,
                                genome= "hg19",
                                chromosome = chrM,
                                name= "ARX binding sites",
                                stacking = "dense",
                                colour = "red")

##DNase seq track

dnaseSeqMouse <- import(con = "~/DataFiles/ChIPseq/Mouse/cerebrummousednaseseq.bigWig", which =
                          GRanges(chrM, IRanges(1, end = .Machine$integer.max - 1)))%>% DataTrack(name = "DNase Seq Track",
                                                                                                  type = "l",
                                                                                                  color = "orange",
                                                                                                  max = 80)

##enhancer Track
enhancerTrack <- import( con = "~/DataFiles/Enhancer Tracks/Mouse/Enhanceresmm9.bed",  which =
                           GRanges(chrM, IRanges(1, end = .Machine$integer.max - 1))) %>%AnnotationTrack(name = "Enhancers",
                                                                                                         stacking= "dense")


##Replace Chipconfirmed sites
plotTracks(list(ideogramTrackMouse, knownGenesMouse, enhancerTrack, MouseARXmotifs, dnaseSeqMouse),
           to = toM, 
           from= fromM)
