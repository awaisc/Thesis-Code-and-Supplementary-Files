## trying to import the chipseq Data

library(Biobase)
library(GEOquery)
library(limma)
library(magrittr)

# load series and platform data from GEO
gset <- getGEO("GSE29985", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]
chipDataSet<-cbind(as.data.frame(gset@featureData@data$`Chromosome annotation`), exprs(gset))
colnames(chipDataSet)<- c("probe", "n2a1", "n2a2", "n2a3", "embyrobrain1", "embyrobrain2", "embyrobrain3")
chipFilterNa<-chipDataSet[!chipDataSet$probe=="",]
splitColoumnComma<-data.frame(do.call('rbind', strsplit(as.character(chipDataSet$probe),',')))
splitColoumnSlash<-data.frame(do.call('rbind', strsplit(as.character(splitColoumnComma$X1),'///')))


spiltColumnSpace<-data.frame(do.call('rbind', strsplit(as.character(splitColoumnComma$X2),' ',fixed=TRUE)))
spiltColumnDotDot<-data.frame(do.call('rbind', strsplit(as.character(spiltColumnSpace$X3),'..',fixed=TRUE)))

spiltColumnOpener<-data.frame(do.call('rbind', strsplit(as.character(spiltColumnDotDot$X1),'(',fixed=TRUE)))
spiltColumnCloser<-data.frame(do.call('rbind', strsplit(as.character(spiltColumnDotDot$X2),')',fixed=TRUE)))

probeAnnotationDataFrame<-cbind(splitColoumnSlash[1],spiltColumnOpener[2], spiltColumnCloser[1], spiltColumnSpace[2])
colnames(probeAnnotationDataFrame)<- c("chromosome","start", "end", "probeID")
dataFrameGrange<-cbind(probeAnnotationDataFrame, chipFilterNa[2:7])
dataFrameGrange$start <- as.character(dataFrameGrange$start)%>%as.numeric()
dataFrameGrange$end <- as.character(dataFrameGrange$end)%>%as.numeric()

##The CHIPSEQDATA in a dataframe
chiPSeqRawDataFrame<-makeGRangesFromDataFrame(na.omit(dataFrameGrange),
                         keep.extra.columns=TRUE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid",1),
                         start.field="start",
                         end.field=c("end", "stop"),
                         starts.in.df.are.0based=FALSE)


### Checking the variability between samples!
grangeRawN2a1<- chiPSeqRawDataFrame[,2]
grangeRawN2a2<-chiPSeqRawDataFrame[,3]
grangeRawN2a3<-chiPSeqRawDataFrame[,4]

grangeRawBrain1<-chiPSeqRawDataFrame[,5]
grangeRawBrain2<-chiPSeqRawDataFrame[,6]
grangeRawBrain3<-chiPSeqRawDataFrame[,7]

library(Gviz)
library(magrittr)
grangeRawN2a1Track<-grangeRawN2a1%>%DataTrack(type="l", name="N2a Sample 1")
grangeRawN2a2Track<-grangeRawN2a2%>%DataTrack(type="l", name="N2a Sample 2")
grangeRawN2a3Track<-grangeRawN2a3%>%DataTrack(type="l", name="N2a Sample 3")

grangeRawBrain1Track<-grangeRawBrain1%>%DataTrack(type="l", name="Embyronic Brain 1")
grangeRawBrain2Track<-grangeRawBrain2%>%DataTrack(type="l", name="Embyronic Brain 2")
grangeRawBrain3Track<-grangeRawBrain3%>%DataTrack(type="l", name="Embyronic Brain 3")



dataFrameN2aAverage<-cbind(dataFrameGrange[1:3], rowMeans(dataFrameGrange[5:7], na.rm = TRUE))
DataFrameBrainAverage<-cbind(dataFrameGrange[1:3], rowMeans(dataFrameGrange[8:10], na.rm = TRUE))

grangeAverageN2aRaw<-makeGRangesFromDataFrame(na.omit(dataFrameN2aAverage),
                                              keep.extra.columns=TRUE,
                                              ignore.strand=TRUE,
                                              seqinfo=NULL,
                                              seqnames.field=c("seqnames", "seqname",
                                                               "chromosome", "chrom",
                                                               "chr", "chromosome_name",
                                                               "seqid",1),
                                              start.field="start",
                                              end.field=c("end", "stop"),
                                              starts.in.df.are.0based=FALSE)


grangeAverageBrainRaw<-makeGRangesFromDataFrame(na.omit(DataFrameBrainAverage),
                                                keep.extra.columns=TRUE,
                                                ignore.strand=TRUE,
                                                seqinfo=NULL,
                                                seqnames.field=c("seqnames", "seqname",
                                                                 "chromosome", "chrom",
                                                                 "chr", "chromosome_name",
                                                                 "seqid",1),
                                                start.field="start",
                                                end.field=c("end", "stop"),
                                                starts.in.df.are.0based=FALSE)



##Making the Averages into Data Tracks
averageBrainRawTrack<-grangeAverageBrainRaw%>%DataTrack(type="histogram", name = "Average Brain")
averageN2aRawTrack<-grangeAverageN2aRaw%>%DataTrack(type="histogram", name="Average N2a")

plotTracks(c(averageN2aRawTrack, averageBrainRawTrack))

##Checking to see where the plot falls. 

## load packages
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)

## Manually changing the Seqlevels
i<-1
df2<-NULL
df3<-NULL
df<-as.data.frame(grangeAverageBrainRaw)
for(i in 1:19){
df2<-cbind(paste0("chrX"), df[df[1]==paste0("Chromosome X"),])
colnames(df2)<-c(1,2,3,4,5,6,7)
df3<-rbind(df3, df2)
i<- i+1
}
seqlevels(SeqlevelsGrangeRaw)
DataFrameSeqLevels<-cbind(df3[1], df3[3:7])

colnames(DataFrameSeqLevels)<- c("chromosome", "start", "end", "width", "strand", "ChipSeqData")
SeqlevelsGrangeRaw<-makeGRangesFromDataFrame(DataFrameSeqLevels, keep.extra.columns = TRUE)
averageInterger<-findOverlaps(SeqlevelsGrangeRaw, grangeN2aChipSeq)%>%countLnodeHits()
avgchipscores<-subset(SeqlevelsGrangeRaw, averageInterger)%>%as.data.frame()
sum(abs(avgchipscores$ChipSeqData))/length(avgchipscores$ChipSeqData)
