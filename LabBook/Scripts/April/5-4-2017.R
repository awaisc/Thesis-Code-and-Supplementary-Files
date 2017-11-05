### Awais writing a script in a methological method to isolate each iteration


function(x){
  library(rtracklayer)
  library(dplyr)
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
                 "chrX",
                 "chrY")[x]
  
  #data import for the chromosome
  gtfUCSC<-import("~/Shiny App tutorials/Tut1/Tutorial1/FullMm9genome.GTF", which =
                    GRanges(chrM, IRanges(1, end = .Machine$integer.max - 1)))%>%as.data.frame()
  ##Arx Motif Location Identification
  library(MotifDb)
  library(magrittr)
  library(BSgenome.Mmusculus.UCSC.mm9)
  Arx6mer<-MotifDb::query(MotifDb, "Arx")[[6]]
  ArxPWM6mer <- round(Arx6mer*100)
  arxMotif<- matchPWM(ArxPWM6mer, BSgenome.Mmusculus.UCSC.mm9[[chrM]], "90%", with.score = TRUE)%>%
    ranges()%>%
    as.data.frame()%>% cbind(seqnames=(chrM))
  
  
  
  ##first we need to isolate each gene
  geneStarts<-gtfUCSC[gtfUCSC$type=="start_codon",]
  geneStops<- gtfUCSC[gtfUCSC$type=="stop_codon", ]
  geneStartStop<-left_join(geneStarts, geneStops, by = c("gene_id" = "gene_id"))
  
  #Use the above to see if a motif is between the sotp and start codon!
  ##figureout to apply to dataframe properly
  function(arxMotif){
  i<-1
  counter<-NULL 
    if(geneStartStop$strand.x=="+"){  
    pastStart<-geneStartStop[geneStartStop$start.x < arxMotif$start[1],] #3first subset for if the motif is located past the start codon
    beforeStop<- pastStart[pastStart$end.y> arxMotif$start[i],] ## take these values and chek to see if they are before the end codon
    numberOfGenes<-rbind(numberOfGenes,beforeStop) ##put these into a dataFrame!!
    counter<- cbind(counter, sum(geneStartStop$end.y> arxMotif$start[1])
    } else{
      pastStart<-geneStartStop[geneStartStop$start.x > arxMotif$start[1],]
      beforeStop<- pastStart[pastStart$end.y< arxMotif$start[1],]
      numberOfGenes<-rbind(numberOfGenes,beforeStop)
      counter<- cbind(counter, sum(geneStartStop$end.y<arxMotif$start[1])
    }
  }
  
  
  
  ##isolating exon sequences
  exons<- gtfUCSC[gtfUCSC$type=="exon",]
  cds<- gtfUCSC[gtfUCSC$type=="CDS",]

  
  ##identifying if the motifs are located in the gene
  if(exons$strand[1]=="+"){
   pastStart<- exons[exons$start<arxMotif$start[1],]
   beforeStop<-pastStart$end>arxMotif$start[1]
   digitcounter<- cbind(digitcounter, sum(beforeStop))
  }else{
    pastStart<- exons$start>arxMotif$start[1]
  beforeStop<-exons$end<arxMotif$start[1]
  digitcounter<- cbind(digitcounter, sum(beforeStop))
  }

  ##identfying if the motiofs are located in the CDS regions
  digitcounter<-NULL
  if(cds$strand[1]=="+"){
    pastStart<- cds[cds$start<arxMotif$start[1],]
    beforeStop<-pastStart$end>arxMotif$start[1]
    digitcounter<- cbind(digitcounter, sum(beforeStop))
  }else{
    pastStart<- cds$start>cds$start[1]
    beforeStop<-exons$end<arxMotif$start[1]
    digitcounter<- cbind(digitcounter, sum(beforeStop))
  }
  
  
  return(sum(digitcounter))
}