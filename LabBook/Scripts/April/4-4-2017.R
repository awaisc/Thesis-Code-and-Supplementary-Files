### Motif Identificaiton
library(MotifDb)
library(seqLogo)
library(magrittr)
library(BSgenome.Mmusculus.UCSC.mm9)
chrM<- "chrX"
mouseArx6mer<-query(MotifDb, "Arx")[[6]]
                      mouseArxPWM6mer <- round(mouseArx6mer*100)
                      mouseArxDataFrame6mer<- matchPWM(mouseArxPWM6mer, BSgenome.Mmusculus.UCSC.mm9[[chrM]], "90%", with.score = TRUE)%>%
                        ranges()%>%
                        as.data.frame()%>% cbind(seqnames=("chrX"))
                      grangesMouseArx6mer<- makeGRangesFromDataFrame(mouseArxDataFrame46er,
                                                                 keep.extra.columns=TRUE,
                                                                 ignore.strand=FALSE,
                                                                 seqinfo=NULL,
                                                                 seqnames.field="seqnames",
                                                                 start.field="start",
                                                                 end.field=c("end", "stop"),
                                                                 strand.field="strand",
                                                                 starts.in.df.are.0based=FALSE)
                      MouseARXmotifs6mer<-AnnotationTrack(grangesMouseArx6mer,
                                                      genome= "hg19",
                                                      chromosome = "chrX",
                                                      name= "ARX binding sites",
                                                      stacking = "dense",
                                                      colour = "red")
                  
  
  ## cluster Identification!
                 function(x, y){   
                   q<- 2
                      w<-1
                      distances <- NA
                      for(i in 1:dim(mouseArxDataFrame6mer)[1]){
                        distances<-cbind(distances ,mouseArxDataFrame6mer[q,1]-mouseArxDataFrame6mer[w,2])
                        q<- q+1
                        w<- w+1
                      }
                      
                      
                      dataFrameDistances<-list(t(distances))%>%as.data.frame()
                      combinedDataFrameDistances<- cbind(mouseArxDataFrame6mer, dataFrameDistances[1:dim(dataFrameDistances)[1]-1,])
                      combinedDataFrameDistances<- cbind(combinedDataFrameDistances, (combinedDataFrameDistances[5]<=200))
                      names(combinedDataFrameDistances) <- c("start", "end", "width", "seqnames", "distancebetween", "clustered")
                      
                      grangeDataFrame<-combinedDataFrameDistances[combinedDataFrameDistances$clustered==TRUE,]
                      grangeDataFrame2<- grangeDataFrame[grangeDataFrame$distance>=0,]
                      y<-makeGRangesFromDataFrame(grangeDataFrame2[-1,],
                                                                  seqinfo=NULL,
                                                                  seqnames.field="seqnames",
                                                                  start.field="start",
                                                                  end.field=c("end", "stop"),
                                                                  strand.field="strand",
                                                                  starts.in.df.are.0based=FALSE) %>%AnnotationTrack(genome= "mm9",
                                                                                                                    chromosome = chrM,
                                                                                                                    name= "6mer Cluster binding sites",
                                                                                                                    stacking = "dense",
                                                                                                                    colour = "red")
                      }
                      

## Script for Identification of locaiton of ARX Motifs
library(rtracklayer)                

gtfUCSCgeneTrack<-import("~/Shiny App tutorials/Tut1/Tutorial1/FullMm9genome.GTF", which =## below i have subsetted for a chrM! okay cool so we only have genes on X
                           GRanges(chrM, IRanges(1, end = .Machine$integer.max - 1)))%>%as.data.frame()## because data frames are easy to manipulate

inGene<-NULL#emptying the filters and intitalising stuff for the the for loop
addme<-NULL
i<- 1
for(i in 1:length(gtfUCSCgeneTrack$type=="start_codon")){##Start of for loop, the dim is an to change how many iterations based on the number of start codons in the chromosome
if(as.data.frame(gtfUCSCgeneTrack$strand[gtfUCSCgeneTrack$type=="start_codon"])[i,]=="+"){#on the plus strand the Arx motif start site if its bigger then the start codon and less than the stop its all good
inGene<-as.data.frame(gtfUCSCgeneTrack$start[gtfUCSCgeneTrack$type=="start_codon"])[1,]< mouseArxDataFrame6mer[1] &  mouseArxDataFrame6mer[1] < as.data.frame(gtfUCSCgeneTrack$start[gtfUCSCgeneTrack$type=="stop_codon"])[i,]
addme<-cbind(addme, inGene)
print("plus strand")
}else{
  if(as.data.frame(gtfUCSCgeneTrack$strand[gtfUCSCgeneTrack$type=="start_codon"])[i,]=="-")## if on the negative its reveresd it needs to be smaller than the stop as its down stream and bigger than the start codon as its upstream
    inGene<-as.data.frame(gtfUCSCgeneTrack$start[gtfUCSCgeneTrack$type=="start_codon"])[1,]> mouseArxDataFrame6mer[1] &  mouseArxDataFrame6mer[1] > as.data.frame(gtfUCSCgeneTrack$start[gtfUCSCgeneTrack$type=="stop_codon"])[i,]
  addme<-cbind(addme, inGene)#piling the data into a dataframme
  print("minus strand")#beucause why not?
}
    
i<- i+1
print(i)
}









i<-1
for(i in i:25) {
chrM<- c("chr1",
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
  "chr20",
  "chr21",
  "chr22",
  "chr23",
  "chrX",
  "chrY")[i]

gtfUCSCgeneTrack<-import("~/Shiny App tutorials/Tut1/Tutorial1/FullMm9genome.GTF", which =## below i have subsetted for a chrM! okay cool so we only have genes on X
                           GRanges(chrM, IRanges(1, end = .Machine$integer.max - 1)))%>%as.data.frame()## because data frames are easy to manipulate

##building a new dataframe
geneStopStartCodon<-NULL
ucscGeneId<- as.data.frame(unique(gtfUCSCgeneTrack$gene_id))
i<- 1
for(i in i:dim(allStartCodons)[1]){
allStartCodons<-gtfUCSCgeneTrack[gtfUCSCgeneTrack$type=="start_codon",]
allStopCodons<- gtfUCSCgeneTrack[gtfUCSCgeneTrack$type=="stop_codon", ]
minGene<-allStartCodons[allStartCodons$gene_id==as.character(ucscGeneId[i,]),]
maxGene<-allStopCodons[allStopCodons$gene_id==as.character(ucscGeneId[i,]), ]
geneStopStartCodon<-rbind(geneStopStartCodon, maxGene, minGene)
i<-1+i

}

##now to write a chunk of code that excutes the same thing as above but on this new dataframe
holder<-NULL
i<- 1
for(i in 1:length(geneStopStartCodon$start)){
  if(as.data.frame(geneStopStartCodon$strand[geneStopStartCodon$type=="start_codon"])[i,]=="+"){
    insideGene<-geneStopStartCodon$start[geneStopStartCodon$type=="start_codon"][i]<mouseArxDataFrame6mer[1] & mouseArxDataFrame6mer[1] <geneStopStartCodon$start[geneStopStartCodon$type=="stop_codon"][i]

    holder<- cbind(holder, insideGene)
  } else{
    insideGene<-geneStopStartCodon$start[geneStopStartCodon$type=="start_codon"][i] > mouseArxDataFrame6mer[1] & mouseArxDataFrame6mer[1] > geneStopStartCodon$start[geneStopStartCodon$type=="stop_codon"][i]

    holder<- cbind(holder, insideGene)
    print(holder)
  }
  i<-i+1
}
print(sum(holder))
i<- i+1
}
##IT FUCKING WORKS!!!!!!! now i need to figure out how to loop it efficetively



###