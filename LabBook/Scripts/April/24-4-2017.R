##Images for Power Point

library(MotifDb)
library(seqLogo)

arxJolma<- query(motifDb, "Arx")[[1]]

arx4Mer<-rbind(A=c(0, 1, 1,0), C=c(0,0,0,0), G=c(0,0,0,0),
               T=c(1,0,0,1))
arx4MerTFBS<-matchPWM(arx4Mer, BSgenome.Mmusculus.UCSC.mm9, "100%")
arxJolma<-MotifDb::query(MotifDb, "Mmusculus-jolma2013-Arx")[[1]]
arxJolma <- round(arxJolma*100)
SeqLogo(ArxJolma)
seqLogo::seqLogo(arxJolma)

arx6Mer<-MotifDb::query(MotifDb, "Arx")[[6]]
seqLogo(arx6Mer)
matchPWM(arx6Mer, BSgenome.Mmusculus.UCSC.mm9, "90%")


##4mer locations
##importing all locations
library(rtracklayer)
library(magrittr)
library(BSgenome.Mmusculus.UCSC.mm9)
library(MotifDb)
library(readxl)
n2aChipSeq <- read_excel("~/DataFiles/ChIPseq/Mouse/ChIPseqDataQuille2011.xls", 
                         sheet = "only N2a")%>%as.data.frame
embyroChipSeq  <- read_excel("~/DataFiles/ChIPseq/Mouse/ChIPseqDataQuille2011.xls", sheet = "only emb brain")%>%as.data.frame
commonChipSeq <- read_excel("~/DataFiles/ChIPseq/Mouse/ChIPseqDataQuille2011.xls",  sheet = "common genes")%>%as.data.frame



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
                                          keep.extra.columns=FALSE,
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


# 
grangeN2aChipSeq<-chipSeqDataCleaner(n2aChipSeq)
grangeBrainChipSeq<-chipSeqDataCleaner(embyroChipSeq)
grangeCommonChipSeq<-chipSeqDataCleaner(commonChipSeq)

subet<-findOverlaps(grangeN2aChipSeq, arx4MerTFBS)%>%countLnodeHits()
length(subset(grangeN2aChipSeq, subet))


subet<-findOverlaps(grangeBrainChipSeq, arx4MerTFBS)%>%countLnodeHits()
length(subset(grangeBrainChipSeq, subet))

subet<-findOverlaps(grangeCommonChipSeq, arx4MerTFBS)%>%countLnodeHits()
length(subset(grangeCommonChipSeq, subet))



