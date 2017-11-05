##importing all locations
library(rtracklayer)
library(magrittr)
library(BSgenome.Mmusculus.UCSC.mm9)
library(MotifDb)
library(readxl)
library(Biostrings)
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





## Our confirmed PWM

arx6MerPWM4space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0.25,0.25,1,0,0,1,1,0),
                         C=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25),
                         G=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25) ,
                         T=c(1,0,0,1,1,0,0.25,0.25,0.25,0.25,0,1,1,0,0,1))
arxTandem2Space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0,1,1,0,0,1),
                        C=c(0,0,0,0,0,0,0.25,0.25,0),
                        G=c(0,0,0,0,0,0,0.25,0.25,0) ,
                        T=c(1,0,0,1,1,0,0.25,0.25,1,0,0,1,1,0))
arxTandem6Space<-rbind( A=c(0,1,1,0,0,1,0.25,0.25,0.25,0.25,0.25,0.25,0,1,1,0,0,1), 
                        C=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0.25,0),
                        G=c(0,0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0.25,0) ,
                        T=c(1,0,0,1,1,0,0.25,0.25,0.25,0.25,0.25,0.25,1,0,0,1,1,0))


grangeplaindromic4space <-
  matchPWM(arx6MerPWM4space, BSgenome.Mmusculus.UCSC.mm9, "100%")
grangeTandem2space <-
  matchPWM(arxTandem2Space, BSgenome.Mmusculus.UCSC.mm9, "100%")
grangeTandem6space <-
  matchPWM(arxTandem6Space, BSgenome.Mmusculus.UCSC.mm9, "100%")

c(grangeplaindromic4space, grangeTandem2space, grangeTandem6space)
##Removing the weirdly massive Peak.
grangeN2aChipSeq<-subset(grangeN2aChipSeq, !names==512)
Datable<-rbind(cbind(
  length(findOverlaps(grangeN2aChipSeq, grangeplaindromic4space)),
  findOverlaps(grangeBrainChipSeq, grangeplaindromic4space)%>%length(),
  findOverlaps(grangeCommonChipSeq, grangeplaindromic4space)%>%length()
),cbind(
  length(findOverlaps(grangeN2aChipSeq, grangeTandem2space)),
  findOverlaps(grangeBrainChipSeq, grangeTandem2space)%>%length,
  findOverlaps(grangeCommonChipSeq, grangeTandem2space)%>%length()
),cbind(
  length(findOverlaps(grangeN2aChipSeq, grangeTandem6space)),
  findOverlaps(grangeBrainChipSeq, grangeTandem6space)%>%length(),
  findOverlaps(grangeCommonChipSeq, grangeTandem6space)%>%length()
))

## Seeing how many probes overlap with our motifs
gff<-gff%>%as.data.frame()
gff<-gff[,2:4]

colnames(gff)<- c( "chromosome", "start" ,"end")
gff$chromosome<-gff$chromosome%>%as.character()
gff$start<-gff$start%>%as.character()%>%as.numeric()
gff$end<-gff$end%>%as.character()%>%as.numeric()
gff<-na.omit(gff)
gff<-subset(gff, gff$end-gff$start>0)
grange<-makeGRangesFromDataFrame(gff)
seqlevelsStyle(grange)<-"UCSC"


##How many probes actually have our motif models in them

cbind(
  findOverlaps(grangeTandem6space, grange),
  findOverlaps(grangeTandem2space, grange),
  findOverlaps(grangeplaindromic4space, grange))






