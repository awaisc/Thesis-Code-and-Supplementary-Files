### cluster identification

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
                                genome= "mm9",
                                chromosome = chrM,
                                name= "ARX binding sites",
                                stacking = "dense",
                                colour = "red")


#in order for this to work right i need 1. End of ARX motif is less than 200bp form the next hence start
q<- 2
w<-1
distances <- NA
for(i in 1:dim(mouseArxDataFrame)[1]){
distances<-cbind(distances ,mouseArxDataFrame[q,1]-mouseArxDataFrame[w,2])
print(w)
q<- q+1
w<- w+1
}

dataFrameDistances<-list(t(distances))%>%as.data.frame()
combinedDataFrameDistances<- cbind(mouseArxDataFrame, dataFrameDistances[1:dim(dataFrameDistances)[1]-1,], (stored[4]<=200)  )
names(combinedDataFrameDistances) <- c("start", "end", "width", "seqnames", "distancebetween", "clustered")

grangeDataFrame<-clustered[clustered$clustered==TRUE,]
grangeDataFrame2<- grangeDataFrame[grangeDataFrame$distance>=0,]
clustereDataTrack<-makeGRangesFromDataFrame(grangeDataFrame2[-1,],
                                                                   seqinfo=NULL,
                                                                   seqnames.field="seqnames",
                                                                   start.field="start",
                                                                   end.field=c("end", "stop"),
                                                                   strand.field="strand",
                                                                   starts.in.df.are.0based=FALSE) %>%AnnotationTrack(genome= "mm9",
                                                                                                                     chromosome = chrM,
                                                                                                                     name= "ARX binding sites",
                                                                                                                     stacking = "dense",
                                                                                                                     colour = "red")
### Conservation Tracks


library(rtracklayer)
library(GenomicRanges)


session <- browserSession()
genome(session) <- "mm9"
trackNames(session) ## list the track names
## choose the Conservation track for a portion of mm9 chr1
query <- ucscTableQuery(session, "Conservation",
                        GRangesForUCSCGenome("mm9", "chrX",
                                             IRanges(1, .Machine$integer.max - 1)))
## list the table names
tableNames(query)
## get the phastCons30way track
tableName(query) <- "phastConsElements30way"
## retrieve the track data
conservationTrack<-track(query) # a GRanges object
## get a data.frame summarizing the multiple alignment
tableName(query) <- "multiz30way"
object<-getTable(query)
dataFrameGrange<- as.data.frame(cbind(object[2:4], object[7]))
GrangeConservation<- makeGRangesFromDataFrame(dataFrameGrange,
                                              keep.extra.columns=TRUE,
                                              ignore.strand=FALSE,
                                              seqinfo=NULL,
                                              seqnames.field=c("seqnames", "seqname",
                                                               "chromosome", "chrom",
                                                               "chr", "chromosome_name",
                                                               "seqid"),
                                              start.field="start",
                                              end.field=c("end", "stop"),
                                              strand.field="strand",
                                              starts.in.df.are.0based=FALSE)%>% DataTrack(type = "p")





##gene region Tracks
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(GenomicFeatures)
mouse<-GeneRegionTrack(TxDb.Mmusculus.UCSC.mm9.knownGene, stacking = "dense", chromosome = chrM)
plotTracks(mouse, to = 5100000, from = 5400000)


##Biomart gene regions
biomTrack <- BiomartGeneRegionTrack(genome = "mm9",
                                     chromosome = chrM,
                                     name = "ENSEMBL Genes")
plotTracks(biomTrack, from = fromM, to = toM)


##filter for differentially expressed genes
#1 get a Grange object for Genes from UCSC,
genes <- ucscTableQuery(session, "UCSC Genes",
                        GRangesForUCSCGenome("mm9", "chrX",
                                             IRanges(1, .Machine$integer.max - 1)))
tableName(genes) <- "knownGene"
UCSC<-getTable(genes)
ucscGrange<- makeGRangesFromDataFrame(UCSC)

##

pointer<-import(con= "~/DataFiles/Tessa differential Expressed/PA16col.bed")%>%as.data.frame()

library(readr)
geneSymbols <- read_tsv("~/DataFiles/Tessa differential Expressed/2colmultiplenames.txt")
                                 


library(dplyr)
differentiallyExpressed <- left_join(pointer, geneSymbols, by = c("name"="#kgID"))

selectCols <- cbind(differentiallyExpressed[1:3],differentiallyExpressed[8])
toM<- 10000+selectCols[3,input$gene]






