## boring table Script


library(rtracklayer)
library(GenomicRanges)
library(MotifDb)
library(BSgenome.Mmusculus.UCSC.mm9)
library(readxl)
library(ggplot2)
library(magrittr)
arx6Mer<- query(MotifDb, "Arx")[[6]]
arx6Mer<- round(arx6Mer*100)
arx6MerLow<- matchPWM(arx6Mer, BSgenome.Mmusculus.UCSC.mm9, "70%" )

arxJolma<- query(MotifDb, "Arx")[[3]]
arxJolma<- round(100*arxJolma)
arxJolmaLow<- matchPWM(arxJolma, BSgenome.Mmusculus.UCSC.mm9, "70%")

##Importing the Data to see lengths
n2aChipSeq <- read_excel("~/DataFiles/ChIPseq/Mouse/ChIPseqDataQuille2011.xls", 
                         sheet = "only N2a")%>%as.data.frame
embyroChipSeq  <- read_excel("~/DataFiles/ChIPseq/Mouse/ChIPseqDataQuille2011.xls", sheet = "only emb brain")%>%as.data.frame
commonChipSeq <- read_excel("~/DataFiles/ChIPseq/Mouse/ChIPseqDataQuille2011.xls",  sheet = "common genes")%>%as.data.frame
##Convert to Grange and clean

chipSeqDataCleaner<-function(x){
  splitColoumnMinus<-data.frame(do.call('rbind', strsplit(as.character(x$location),'-',fixed=TRUE)))
  colnames(splitColoumnMinus)<- c("X1", "end")#re naming the coloumns
  splitColoumnSemiColon<-data.frame(do.call('rbind', strsplit(as.character(splitColoumnMinus$X1),':',fixed=TRUE)))
  colnames(splitColoumnSemiColon)<- c("chromosome", "start")#renaming those two
  
  geneSymbolMetaDataFromOriginalData<-x[2:3]
  
  dataFrameOfChipSeqData<- cbind(geneSymbolMetaDataFromOriginalData, splitColoumnSemiColon, splitColoumnMinus[2])%>%na.omit()
  removingTheNegatives<- cbind(dataFrameOfChipSeqData, (as.data.frame(as.numeric(as.character(dataFrameOfChipSeqData$end)))-as.data.frame(as.numeric(as.character(dataFrameOfChipSeqData$start)))))
  
  
  negativesRemoved<-subset(removingTheNegatives, removingTheNegatives$`as.numeric(as.character(dataFrameOfChipSeqData$end))`>0)
  grangeChipSeq<-makeGRangesFromDataFrame(negativesRemoved,
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
                                          starts.in.df.are.0based=FALSE)
}
grangeN2aChipSeq<-chipSeqDataCleaner(n2aChipSeq)
grangeBrainChipSeq<-chipSeqDataCleaner(embyroChipSeq)
grangeCommonChipSeq<-chipSeqDataCleaner(commonChipSeq)

##how many 6mers are in the cihp confirmed promoters

intergerN2a6mer<-findOverlaps(arx6MerLow, grangeN2aChipSeq)%>%countRnodeHits()
genesMotifArxN2a6mer<-length(subset(grangeN2aChipSeq, intergerN2a6mer))
totalNumberGenesN2a<-length(grangeN2aChipSeq)

intergerCommom6mer<-findOverlaps(arx6MerLow, grangeCommonChipSeq)%>%countRnodeHits()
genesMotifArxCommomN2a6mer<-length(subset(grangeCommonChipSeq, intergerCommom6mer))
totalNumberGenesCommon<-length(grangeCommonChipSeq)

intergerBrain6mer<-findOverlaps(arx6MerLow, grangeBrainChipSeq)%>%countRnodeHits()
genesMotifArxBrain6mer<-length(subset(grangeBrainChipSeq, intergerBrain6mer))
totalNumberGenesBrain<-length(grangeBrainChipSeq)
totalNumberArxMotifs6mer<-length(arx6MerLow)

#how many Jolmas are in the chipconfirmed promoters!
intergerBrainJolma<-findOverlaps(arxJolmaLow, grangeBrainChipSeq)%>%countRnodeHits()
genesMotifArxBrainJolma<-length(subset(grangeBrainChipSeq, intergerBrainJolma))
totalNumberGenesBrain<-length(grangeBrainChipSeq)

intergerCommonJolma<-findOverlaps(arxJolmaLow, grangeCommonChipSeq)%>%countRnodeHits()
genesMotifArxCommonJolma<-length(subset(grangeCommonChipSeq, intergerCommonJolma))
totalNumberGenesCommon<-length(grangeCommonChipSeq)

intergerN2aJolma<-findOverlaps(arxJolmaLow, grangeN2aChipSeq)%>%countRnodeHits()
genesMotifArxN2aJolma<-length(subset(grangeN2aChipSeq, intergerN2aJolma))
totalNumberGenesN2a<-length(grangeN2aChipSeq)
totalNumberArxMotifsJolma<-length(arxJolmaLow)









