###
library(rtracklayer)
library(magrittr)
library(BSgenome.Mmusculus.UCSC.mm9)
library(MotifDb)
library(realxl)

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
grangeCommonChipSeq<-chipSeqDataCleaner(embyroChipSeq)
chipcConfirmedSites<- grangeCommonChipSeq

##inputs for the for loop
arx6Mer<-MotifDb::query(MotifDb, "Arx")[[6]]
arx6Mer <- round(arx6Mer*100)
arx6MerTFBS<- matchPWM(arx6Mer, BSgenome.Mmusculus.UCSC.mm9, "90%", with.score = TRUE)
transcriptionFactorModel<- arx6MerTFBS





i<-1
storage<-NULL
for(i in 1:21){
  chromosome<-c("chr1",
                "chr2", 
                "chr3",
                "chr4",
                "chr5",
                "chr6",
                "chr7",
                "chr8",
                "chr9",
                "chr10",
                "chr11",
                "chr12",
                "chr13",
                "chr14",
                "chr15",
                "chr16",
                "chr17",
                "chr18",
                "chr19",
                "chrX",
                "chrY")[i]
  
  dnaseSeqMouse <- import(con = "~/DataFiles/DNase/Mouse/ENCFF292LVM.bigWig", which =  GRanges(chromosome, IRanges(1, end = .Machine$integer.max - 1)))
  arx6chrM<-subset(transcriptionFactorModel, seqnames==chromosome)
  
  seqlengths(dnaseSeqMouse)= seqlengths(arx6chrM)[names(seqlengths(dnaseSeqMouse))]
  seqlevels(arx6chrM, force=TRUE)= seqlevels(dnaseSeqMouse)
  intergers<-nearest(arx6chrM, dnaseSeqMouse)
  grangesArxMotifWithDnaseSeq<-dnaseSeqMouse[intergers,]
  tfbsWithScore<-cbind(arx6chrM, as.data.frame(grangesArxMotifWithDnaseSeq)[6])%>%na.omit()%>%makeGRangesFromDataFrame(
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
  
  
  average2aDnaseI<-import(con = "~/DataFiles/DNase/Mouse/ENCFF292LVM.bigWig", which =chipcConfirmedSites)%>%as.data.frame()
  averageN2a<-sum(average2aDnaseI[6])/dim(average2aDnaseI)[1]
  potentialMotifs<-subset(tfbsWithScore, score.1>averageN2a)
  storage<-rbind(storage, as.data.frame(potentialMotifs))
  i<-i+1
  print(storage)
}
grangeStorage<-makeGRangesFromDataFrame(storage, keep.extra.columns = TRUE)


#ACtive Enhancers
enhancerGrange <- import( con = "~/DataFiles/Enhancer Tracks/Mouse/Enhanceresmm9.bed")

enhancerWithPotentialMotifs<-findOverlaps(enhancerGrange, grangeStorage)%>%countLnodeHits()
ActiveArxEnhancers<-subset(enhancerGrange, enhancerWithPotentialMotifs)




##seeing if there is an DnaseSeq SCore between and motif score!
findOverlaps(grangeCommonChipSeq, grangeStorage)
