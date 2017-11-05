#mm9 Script

library(rtracklayer)
library(magrittr)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(AnnotationHub)
library(annotation)

genome<-BSgenome.Hsapiens.UCSC.hg19
ConservationScoreCuttOff<-0.9
InteractionStrength<-5

enhancers<-import("~/DataFiles/Enhancer Tracks/Human/human_permissive_enhancers_phase_1_and_2.bed")
genes<-import("~/DataFiles/Gene Tracks/Human/hg.bed")
promoters<-promoters(genes)


#Converting ARX binding sites to the mm10 genome for getting thhe scoresobject
arx6mer<-matchPWM( round(PWM("TAATTA")*7), genome, "100%" )

#arx4mer<-matchPWM(round(PWM("AATT")*4), genome, "100%")

conservation<-import("~/DataFiles/Conservation/Human/hg19.100way.phastCons.bw", which= arx6mer)
conservationScore<-subset(conservationLifted, score>=ConservationScoreCuttOff)


MotifModel<-arx6mer
BindingSitespromoters<-subsetByOverlaps(MotifModel, promoters)
BindingSitesenhancers<-subsetByOverlaps(MotifModel, enhancers)

ConservedBindingSitesenhancers<-subsetByOverlaps(BindingSitesenhancers, conservationScore)
ConservedBindingSitespromoters<-subsetByOverlaps(BindingSitespromoters, conservationScore)

conservedBindingSitesEnhancersAndPromoters<-c(ConservedBindingSitesenhancers, ConservedBindingSitespromoters)%>%unlist()



## HiC Training  Data
library(GenomicInteractions)
HiCInteractions<-makeGenomicInteractionsFromFile(fn = "~/DataFiles/HiC/Mouse/Interaction5kbResolution.bedpe", type = "bedpe",
                                                 experiment_name = "Mouse Brain", description = " Experiment 2")

validInteractions<-HiCInteractions[HiCInteractions$counts>=InteractionStrength]

##Getting Promoter Enhancer Interactions
colnames(mcols(enhancers))<-c("id", "score", "itmeRgb", "thick", "blacks")
colnames(mcols(ConservedBindingSitesenhancers))<-c("id", "score")
colnames(mcols(promoters))<-c("id"   , "score",   "itemRgb", "thick",   "blocks" )
genomicFeaturesForHiC <- list(enhancer = ConservedBindingSitesenhancers, 
                              promoters = promoters%>%GRanges())

annotateInteractions(validInteractions, genomicFeaturesForHiC)

PromoterEnhancerInteractions<-validInteractions[isInteractionType(validInteractions, "promoters", "enhancer")]

InteractionToGenomicRanges<-function(x){
  Test<-x%>%as.data.frame()
  Test1<-cbind(Test$seqnames1%>%as.character(), Test$start1, Test$end1, Test$strand1%>%as.character)%>%as.data.frame()
  Test2<-cbind(Test$seqnames2%>%as.character(), Test$start2, Test$end2, Test$strand2%>%as.character)%>%as.data.frame()
  colnames(Test1)<- c("chromosome", 
                      "start",
                      "end",
                      "strand")
  colnames(Test2)<- c("chromosome", 
                      "start",
                      "end",
                      "strand")
  
  Grange1<-makeGRangesFromDataFrame(Test1)
  Grange2<-makeGRangesFromDataFrame(Test2)
  
  seqlevelsStyle(Grange1)<-"ucsc"
  seqlevelsStyle(Grange2)<-"ucsc"
  mcols(Grange2)<-Test$counts
  mcols(Grange1)<-Grange2
  return(Grange1)
}

validInteractionsGranges<-InteractionToGenomicRanges(PromoterEnhancerInteractions)
InteractingPromoters<-subsetByOverlaps(promoters,mcols(validInteractionsGranges)$X)%>%GRanges

#Creating the seqInfo Object
mm9Seqinfo<-seqinfo(TxDb.Mmusculus.UCSC.mm9.knownGene)
seqlevels(mm9Seqinfo)<-seqlevels(validInteractionsGranges)


seqinfo(validInteractionsGranges)<-mm9Seqinfo

seqlevels(InteractingPromoters, force=TRUE)<-seqlevels(mm9Seqinfo)
seqinfo(InteractingPromoters)<-mm9Seqinfo

seqlevels(enhancercircle, force=TRUE)<-seqlevels(mm9Seqinfo)
seqinfo(enhancercircle, force= TRUE)<-mm9Seqinfo


##Gene Symbols instead of UCSC Gene IDS
Mus.muculsus<-makeOrganismDbFromTxDb(TxDb.Mmusculus.UCSC.mm9.knownGene,orgdb = org.Mm.eg.db)
symbolGenes <- transcripts(Mus.muculsus, columns="SYMBOL")
promoterGenes<-subsetByOverlaps(promoters(symbolGenes), 
                                (mcols(validInteractionsGranges))$X )
metaData<-((mcols(promoterGenes))%>%as.matrix())
mcols(promoterGenes)<-metaData

seqlevels(promoterGenes)<-seqlevels(mm9Seqinfo)
seqinfo(promoterGenes)<-mm9Seqinfo

library(ggbio)

p<-ggbio()+circle(validInteractionsGranges, geom= "link", linked.to="X", aes(color=seqnames, y= `X.1`), radius =50)+
  circle(promoterGenes, geom= "text", aes(color=seqnames, label= `SYMBOL`), angle= 90, size = 3)+
  circle(promotersUnique, geom= "point", aes(color=`X`), radius= 49, color="red")+
  circle(enhancercircle, geom="point", aes(color=`X`), radius=50, color="blue")+
  circle(InteractingPromoters, geom= "ideogram", aes(fill=seqnames))+
  circle(InteractingPromoters, geom= "scale", aes(color=seqnames))






### Epigenetic
brainChromHMMMouse<-import("~/DataFiles/ChromHMM/mouse/brain_cStates_HMM.bed")

unique((mcols(brainChromHMMMouse))$name)


