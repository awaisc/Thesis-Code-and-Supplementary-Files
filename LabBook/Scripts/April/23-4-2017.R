##Sequence enrichment SCript
library(readxl)
library(magrittr)
library(GenomicRanges)
library(tibble)

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


arx6Mer <-
  rbind(
    A = c(0, 1, 1, 0, 0, 1),
    C = c(0, 0, 0, 0, 0, 0),
    G = c(0, 0, 0, 0, 0, 0) ,
    T = c(1, 0, 0, 1, 1, 0)
  )


####sequence Enrichement
library(PWMEnrich)
library(BSgenome.Mmusculus.UCSC.mm9)
library(JASPAR2016)
library(rGADEM)


grange6mer<-matchPWM(arx6Mer, BSgenome.Mmusculus.UCSC.mm9, "100%")

##Strings for Sequence Enrichement WITHOUT ARX motifs! 6mer specically

##Motifs in Common genes with no Arx motifs!
overlaps<- distanceToNearest(grangeCommonChipSeq, grange6mer)
nonOverLappingInterger<- subset(overlaps, !distance==0)%>%countLnodeHits()
nonArxContainingCommon<- subset(grangeCommonChipSeq, nonOverLappingInterger)
stringsCommon<-getSeq(BSgenome.Mmusculus.UCSC.mm9,nonArxContainingCommon)

##Motifs in Brain genes with no Arx motifs!
overlapsBrain<- distanceToNearest(grangeBrainChipSeq, grange6mer)
nonOverLappingIntergerBrain<- subset(overlapsBrain, !distance==0)%>%countLnodeHits()
nonArxContainingBrain<- subset(grangeBrainChipSeq, nonOverLappingIntergerBrain)
stringsBrain<-getSeq(BSgenome.Mmusculus.UCSC.mm9,nonArxContainingBrain)

##Motifs in N2a genes with no Arx motifs!
overlapsN2a<- distanceToNearest(grangeN2aChipSeq, grange6mer)
nonOverLappingInterger<- subset(overlapsN2a, !distance==0)%>%countLnodeHits()
nonArxContainingN2a<- subset(grangeN2aChipSeq, nonOverLappingInterger)
stringsN2a<-getSeq(BSgenome.Mmusculus.UCSC.mm9,nonArxContainingN2a)


# ## exporting
writeXStringSet(stringsN2a, filepath= "~/DataFiles/ChIPseq/N2aStrings", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
writeXStringSet(stringsBrain, filepath= "~/DataFiles/ChIPseq/stringsBrain", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
writeXStringSet(x=stringsCommon, filepath= "~/DataFiles/ChIPseq/CommonStrings", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

## AFter i run it in WEEDER2.0 Yeah yeah
motifsNotIncludingArx <- read_excel("~/DataFiles/ChIPseq/Motifs.xlsx",col_names = FALSE)%>%as.matrix
rownames(Motifs)<- Motifs[,1]
motifsNumeric<-apply(motifsNotIncludingArx,2, as.numeric)

rownames(motifsNumeric)<- motifsNotIncludingArx[,1]
brainMotif1<-seqLogo(motifsNumeric[2:5, 2:7])
brainMotif2<-seqLogo(motifsNumeric[6:9, 2:11])
brainMotif3<-seqLogo(motifsNumeric[10:13, 2:9])

n2aMotif1<-seqLogo(motifsNumeric[15:18, 2:7])
n2aMotif2<-seqLogo(motifsNumeric[19:22, 2:11])
n2aMotif3<-seqLogo(motifsNumeric[23:26, 2:9])

commonMotif1<-seqLogo(motifsNumeric[28:31, 2:7])
commonMotif2<-seqLogo(motifsNumeric[32:35, 2:11])
commonMotif3<-seqLogo(motifsNumeric[36:39, 2:9])






