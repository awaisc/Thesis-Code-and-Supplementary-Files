##plotting to see if there is an inverse relationship between arx6mer TFBS score and DNaseSeq Score

library(rtracklayer)
library(GenomicRanges)
library(MotifDb)
library(BSgenome.Mmusculus.UCSC.mm9)
library(readxl)
library(ggplot2)

arx6Mer<- query(MotifDb, "Arx")[[6]]
arx6Mer<- round(arx6Mer*100)
arx6MerLow<- matchPWM(arx6Mer, BSgenome.Mmusculus.UCSC.mm9, "70%" )
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
grangeN2aChipSeq<-chipSeqDataCleaner(n2aChipSeq)
grangeBrainChipSeq<-chipSeqDataCleaner(embyroChipSeq)
grangeCommonChipSeq<-chipSeqDataCleaner(commonChipSeq)

grangeMotifs<-subset(arx6MerLow, findOverlaps(grangeN2aChipSeq, arx6MerLow)%>%countRnodeHits())

importDnaseGrange<-resize(grangeMotifs, 1000)
dnaseSeqMouse <- import(con = "~/DataFiles/DNase/Mouse/ENCFF292LVM.bigWig", which = importDnaseGrange)

seqlengths(dnaseSeqMouse)= seqlengths(arx6MerLow)[names(seqlengths(dnaseSeqMouse))]
seqlevels(arx6MerLow, force=TRUE)= seqlevels(dnaseSeqMouse)
intergers<-nearest( grangeMotifs, dnaseSeqMouse)
ggplotDataFrame<-cbind(grangeMotifs, as.data.frame(values(dnaseSeqMouse[intergers,])))%>%makeGRangesFromDataFrame()


colnames(ggplotDataFrame)<- c("seqnames", "start", "end", "motif width",
                              "strand",  "motif score", "DNA string", "DNase" )
logGgplotDataFrame<-cbind(ggplotDataFrame, "logDnase"=log(as.data.frame(ggplotDataFrame$DNase)))
colnames(logGgplotDataFrame)<- c("seqnames", "start", "end", "motif width","strand",
                                 "motif score", "DNA string", "DNase", "logDnase" )















##Code for generating Graphs

ggplot(logGgplotDataFrame, mapping = aes(x=`motif score`, y= DNase))+
     geom_smooth( )+
  ggtitle("Arx motif score vs DNaseI Score")+
     theme_bw()

ggplot(logGgplotDataFrame[logGgplotDataFrame$`motif score`>900,], mapping = aes(x=`motif score`, y= DNase))+
  geom_line( )+
  ggtitle("Arx motif score vs DNaseI Score")+
  theme_bw()
ggplot(logGgplotDataFrame[logGgplotDataFrame$DNase>100,], mapping = aes(x=`motif score`, y= DNase))+
  geom_point()+
  geom_smooth()+
  ggtitle("Arx motif distribution in 15.5dpc brain")
theme_bw()

ggplot(logGgplotDataFrame[logGgplotDataFrame$DNase>100,], mapping = aes(x=`motif score`))+
  geom_freqpoly()+
  ggtitle("Arx motif Distribution in the genome")
  theme_bw()
  
  ggplot(logGgplotDataFrame[logGgplotDataFrame$DNase>100,], mapping = aes(x=DNase))+
    geom_freqpoly()+
    ggtitle("Arx motif distribution in 15.5dpc brain")
    theme_bw()

ggplot(logGgplotDataFrame, mapping = aes(x=`motif score`, y= DNase))+
  geom_point( )+
  ggtitle("Arx motif score vs DNaseI Score")+
  theme_bw()

