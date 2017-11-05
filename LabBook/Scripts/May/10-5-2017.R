

## ChIP-chip Data workflow thingo magico
library(Ringo)
library(magrittr)
library(limma)
library(mclust)
library(ggplot2)



arrayfiles <- list.files(path="/home/a1649239/DataFiles/ChIPseq/Mouse/",
                         pattern="GSM742106_US45103054_251471711563_S01_ChIP")
RG <- read.maimages(arrayfiles, 
                    source="agilent", 
                    path="/media/awais/NewDrivewho/Downloads/chip-chip/")

par(mar=c(0.01,0.01,0.01,0.01), bg="black")
image(RG, 1, channel="red", dim1="Col", dim2="Row",
      mycols=c("sienna","darkred","orangered"))

image(RG,arrayno,channel=c("red","green","logratio"),
      mycols=NULL, mybreaks=NULL, dim1="PositionX", dim2="PositionY",
      ppch=20, pcex=0.3, verbose=TRUE)

pA <- extractProbeAnno(RG, "agilent", genome="mouse",
                       microarray="Agilent Tiling N2a")
X <- preprocess(RG[RG$genes$ControlType==0,], method="nimblegen",
                idColumn="ProbeName")

probeDists <- diff(pA["Y.start"])
br <- c(0, 100, 200, 300, 500, 1000, 10000, max(probeDists))
table(cut(probeDists, br))

##Working with my data up to here

smoothX <- computeRunningMedians(X, modColumn="FileName",
                                 winHalfSize=500, min.probes=3, probeAnno=pA)
sampleNames(smoothX) <- paste(sampleNames(X),"smooth",sep=".")

combX <- combine(X, smoothX)

frameData<-as.data.frame(combX@featureData@data$SystematicName)

chromsomeSplit<-data.frame(do.call('rbind', 
                                   strsplit(as.character(frameData$`combX@featureData@data$SystematicName`),
                                            ':',fixed=TRUE))) 
startSplit<-data.frame(do.call('rbind', 
                               strsplit(as.character(chromsomeSplit$X2),
                                        '-',fixed=TRUE))) 
genes<-RG$genes
gff<-cbind("name",chromsomeSplit$X1, as.character(startSplit$X1)%>%as.numeric, as.character(startSplit$X2)%>%as.numeric,"*")%>%as.data.frame
colnames(gff)<- c("name", "chr", "start", "end", "strand")

ggplot()
plot(combX, pA, chr="X",
     gff=gff,
     maxInterDistance=450, paletteName="Paired")



y0 <- upperBoundNull(exprs(smoothX))
y0G <- twoGaussiansNull(exprs(smoothX), max.adj.p=0.01)
##baseplot
#  hist(exprs(smoothX), n=100000, main=NA,
#         xlab="GSM742106_US45103054_251471711563_S01_ChIP.smooth")
# abline(v=y0, col="red", lwd=2)
# abline(v=y0G, col="blue", lwd=2)
# legend(x="topright", lwd=2, col=c("red","blue"),
#        legend=c("Non-parametric symmetric Null", "Gaussian Null"))
#

##Ggplot
ggplot(as.data.frame(exprs(smoothX)), aes(x=GSM742106_US45103054_251471711563_S01_ChIP.smooth, fill=GSM742106_US45103054_251471711563_S01_ChIP.smooth))+
  geom_histogram(bins = 1000)+
  geom_vline(xintercept=y0G,show.legend = TRUE)+
  geom_vline(xintercept=y0,show.legend = TRUE)+
  theme_bw()
chersX <- findChersOnSmoothed(smoothX, probeAnno=pA, threshold=y0)
gff$start<- gff$start%>%as.character%>%as.numeric%>%as.data.frame()
gff$end<- gff$end%>%as.character%>%as.numeric%>%as.data.frame()
gff$chr<-gff$chr%>%as.character()%>%as.data.frame()
gff<-gff[gff$end-gff$start>0,]

gff<-na.omit(gff)


##converting it to a Grange
library(GenomicRanges)

arxChiPPeaks<-chersX%>%as.data.frame
arxChiPPeaks<-cbind(arxChiPPeaks$chr, arxChiPPeaks$start, arxChiPPeaks$end, arxChiPPeaks$maxLevel, arxChiPPeaks$score)%>%as.data.frame 
colnames(arxChiPPeaks)<-c("chromosome", "start", "end", "max level", "score")



peaks<-arxChiPPeaks%>%na.omit()%>%makeGRangesFromDataFrame(
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
seqlevelsStyle(peaks)<-"UCSC"
##Seeing the overlap of sites containing a 6mer vs not containing a 6mer
library(Biostrings)
library(BiocInstaller)
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

arx6MerTFBS<-matchPWM(arx6Mer, BSgenome.Mmusculus.UCSC.mm9, "100%")

length(subset(x = peaks,findOverlaps(arx6MerTFBS, peaks)%>%countRnodeHits()))

