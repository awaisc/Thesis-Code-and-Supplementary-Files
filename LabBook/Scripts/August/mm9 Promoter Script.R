#mm9 Script

library(rtracklayer)
library(magrittr)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm9)
library(AnnotationHub)

genome<-BSgenome.Mmusculus.UCSC.mm9
ConservationScoreCuttOff<-0.7
InteractionStrength<-3

enhancers<-import("~/DataFiles/Enhancer Tracks/Mouse/mouse_permissive_enhancers_phase_1_and_2.bed")
genes<-import("~/DataFiles/Gene Tracks/Mouse/mm9.bed")
promoters<-promoters(genes, upstream= 10000)

arx6mer<-matchPWM( round(PWM("TAATTA")*7), genome, "100%" )
#arx4mer<-matchPWM(round(PWM("AATT")*4), genome, "100%")

MotifModel<-arx6mer

#Converting ARX binding sites to the mm10 genome for getting thhe scores

ah<-AnnotationHub()
mm10Chain<-AnnotationHub::query(ah, "mm9ToMm10.over.chain.gz")
mm10ChainFile<-mm10Chain[["AH14596"]]
arx6merTFBSmm10<-liftOver(chain = mm10ChainFile, x= arx6mer)%>%unlist

conservation<-import("~/DataFiles/Conservation/Mouse/mm10.60way.phastCons.bw", which= arx6merTFBSmm10)

##Lifting Over the mm10 Conservation Track
mm9ChainFile<-AnnotationHub::query(ah, "mm10ToMm9.over.chain.gz")
mm9Chainfile<-mm9ChainFile[["AH14535"]]

#Conservation Of All Binidng sites Of ARX!
conservationLifted<-liftOver(chain = mm9Chainfile, x = conservation)%>%unlist()
conservationScore<-subset(conservationLifted, score>=ConservationScoreCuttOff)


BindingSitespromoters<-subsetByOverlaps(MotifModel, promoters)
BindingSitesenhancers<-subsetByOverlaps(MotifModel, enhancers)

ConservedBindingSitesenhancers<-BindingSitesenhancers[countOverlaps(BindingSitesenhancers, conservationScore)>=6]
ConservedBindingSitespromoters<-BindingSitespromoters[countOverlaps(BindingSitespromoters, conservationScore)>=6]

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
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(org.Mm.eg.db)
library(OrganismDbi)
mm9Seqinfo<-seqinfo(TxDb.Mmusculus.UCSC.mm9.knownGene)
seqlevels(mm9Seqinfo)<-seqlevels(validInteractionsGranges)

# Converting everything to the same Seqinfo level
seqinfo(validInteractionsGranges)<-mm9Seqinfo
seqlevels(mcols(validInteractionsGranges)$X)<-seqlevels(mm9Seqinfo)
seqinfo(mcols(validInteractionsGranges)$X)<-mm9Seqinfo
seqlevels(InteractingPromoters, force=TRUE)<-seqlevels(mm9Seqinfo)
seqinfo(InteractingPromoters)<-mm9Seqinfo

enhancercircle<-subsetByOverlaps(enhancers, conservedBindingSitesEnhancersAndPromoters)

seqlevels(enhancercircle, force=TRUE)<-seqlevels(mm9Seqinfo)
seqinfo(enhancercircle, force= TRUE)<-mm9Seqinfo


##Gene Symbols instead of UCSC Gene IDS
Mus.muculsus<-makeOrganismDbFromTxDb(TxDb.Mmusculus.UCSC.mm9.knownGene,orgdb = org.Mm.eg.db)
symbolGenes <- transcripts(Mus.muculsus, columns="SYMBOL")
promoterGenes<-subsetByOverlaps(promoters(symbolGenes), 
                                (mcols(validInteractionsGranges))$X )

metadata<-lapply((mcols(promoterGenes))$SYMBOL, as.data.frame)%>%do.call(rbind, .)
colnames(metadata)<-"SYMBOL"
mcols(promoterGenes)<-metadata
seqlevels(promoterGenes)<-seqlevels(mm9Seqinfo)
seqinfo(promoterGenes)<-mm9Seqinfo


library(ggbio)

p<-ggbio()+circle(validInteractionsGranges, geom= "link", linked.to="X", aes(color=seqnames, y= `X.1`), radius = 50)+
  circle(promoterGenes, gap.geom= "chevron", aes(color=seqnames, label= `SYMBOL`), size = 3, width= 5, angle = 90)+
  circle(promoterGenes, geom= "point", aes(color=`X`), radius= 49, color="red")+
  circle(enhancercircle, geom="point", aes(color=`X`), radius=50, color="blue")+
  circle(InteractingPromoters, geom= "ideogram", aes(fill=seqnames))+
  circle(InteractingPromoters, geom= "scale", aes(color=seqnames))
  
p
##Gene Symbols

  

### Epigenetic
ggbioCircularMaker<-function(x){
#allPossible<-unique((mcols(brainChromHMMMouse))$name)%>%as.data.frame()

##Active genes/promoters/Enhancers in the brain
activeRegionsList<-
  c("9_Strong_Enhancer",
  "8_Strong_Enhancer",
  "1_Txn_Elogantion",
  "5_Active_Promoter",
  "6_Strong_Enhancer",
  "7_Active_Promoter",
  "4_Poised_Enhancer",
  "10_Poised_Promoter")


activeRegionsGRange<-sort(subset(x, name %in% activeRegionsList))

BrainActiveEnhancer<-subsetByOverlaps(enhancercircle, activeRegionsGRange, minoverlap = 200)

HiCInteractonsBrainActive<-subsetByOverlaps(validInteractionsGranges, BrainActiveEnhancer)


ActivePromoters<-subsetByOverlaps(promoterGenes, activeRegionsGRange)


BrainActiveGenes<-subsetByOverlaps(ActivePromoters, (mcols(HiCInteractonsBrainActive))$X)

seqlevels(BrainActiveGenes, force=TRUE)<-seqlevels(mm9Seqinfo)
seqinfo(BrainActiveGenes)<-mm9Seqinfo
seqlevels(BrainActiveGenes)<-seqlevels(mm9Seqinfo)
seqinfo(BrainActiveGenes)<-mm9Seqinfo
seqinfo(InteractingPromoters)<-mm9Seqinfo
BrainActiveGenes<-sort(BrainActiveGenes)
BrainActiveEnhancer<-sort(BrainActiveEnhancer)
InteractingPromoters<-sort(InteractingPromoters)


#Mat
HiCInteractonsBrainActive<-sort(HiCInteractonsBrainActive)
BrainActiveGenes<-sort(BrainActiveGenes)
InteractingPromoters<-sort(InteractingPromoters)


cirularPlot<-ggbio()+
  circle(HiCInteractonsBrainActive, geom= "link", linked.to="X",radius = 50)+
  circle(BrainActiveGenes, gap.geom= "rect", which=BrainActiveGenes,  aes(color=seqnames), radius= 49, color="red")+
  circle(BrainActiveEnhancer, gap.geom= "rect", aes(color=seqnames), radius=50, color="blue")+
  circle(HiCInteractonsBrainActive, geom= "ideogram", aes(fill=seqnames))+
  circle(InteractingPromoters, geom= "scale", aes(color=seqnames))

return(cirularPlot)
}


colorDataFrame<-cbind.data.frame(paste0("chr", 1:21), print('=', quote=FALSE), print(ggcggplotColours(n=21)[1:21]))
colorVector<-c("chr1"="#F8766D" ,            
"chr2"= "#EB8335",
"chr3"="#DA8F00",
"chr4"= "#C49A00",
"chr5"="#A9A400" ,
"chr6"= "#86AC00",
"chr7"="#53B400",
"chr8"="#00BA38",
"chr9"= "#00BE6D",
"chr10"="#00C094",
"chr11"="#00C0B5",
"chr12"="#00BDD2",
"chr13"="#00B6EB",
"chr14"="#00ABFD",
"chr15"="#619CFF",
"chr16"="#A58AFF",
"chr17"="#D078FF" ,
"chr18"="#EC69EF",
"chr19"="#FB61D7",
"chrX" ="#FF63B9",
"chrY" ="#FF6B96"
)
## Testies
testChromHmmMouse<-import("~/DataFiles/ChromHMM/mouse/testes_cStates_HMM.bed")
mESCembyronicMouseHmm<-import("~/DataFiles/ChromHMM/mouse/mESC_cStates_HMM.bed")

ggbioCircularMaker(sort(testChromHmmMouse))
ggbioCircularMaker(brainChromHMMMouse)
ggbioCircularMaker(mESCembyronicMouseHmm)




## Data For Heat Map
InactiveRegionList<-c("13_Heterochrom",
                      "14_Heterochrom",
                      "15_Insulator",
                      "9_Txn_Transition",
                      "11_Repressed",
                      "12_Heterochrom",
                      "2_Weak_Txn")


inactiveRegionsGRange<-subset(brainChromHMMMouse, name %in% InactiveRegionList)
conservedBindingSitesEnhancersAndPromoters







