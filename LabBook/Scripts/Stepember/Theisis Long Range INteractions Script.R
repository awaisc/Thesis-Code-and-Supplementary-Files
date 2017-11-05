
## GEnerating tablle of enhancer - promoter (in this case gene Symbol) Table for thesis


library(GenomicInteractions)
library(magrittr)
library(readr)
library(dplyr)


## Read in data tables
GiTrack<-readRDS("~/DataFiles/HiC/Human/StasticallySignificantHg19BetweenPromoterEnhancer")
dataFrameGiTrack<-GiTrack%>%as.data.frame()


## Split Them up based on anchor location eg Anchor 1 is promoter anchor 2 is Enhancer or enhancer is anchor 1 and promoter is anchor 2. 
Spilt1<-subset(dataFrameGiTrack, !promoters.id1=="NA")
Spilt2<-subset(dataFrameGiTrack, !promoters.id2=="NA")


Spilt1<-cbind.data.frame("PromoterIds"=Spilt1$promoters.id1,
                         "EnhancerIds"= Spilt1$enhancers.id2,
                         "Counts" = Spilt1$counts, 
                         "q_value" = Spilt1$SignificantValues.q_value)    

Split2<-cbind.data.frame("PromoterIds"=Spilt2$promoters.id2,
                          "EnhancerIds"= Spilt2$enhancers.id1,
                          "Counts" = Spilt2$counts, 
                          "q_value" = Spilt2$SignificantValues.q_value)    


#Combine the two Datafmraes into 1
PromoterIds<-rbind(Spilt1,Split2)

#Convert from ASIS format to dataframe
PromoterIds$PromoterIds<-lapply(PromoterIds$PromoterIds, list)
PromoterIds$PromoterIds<-lapply(PromoterIds$PromoterIds,function(x){do.call(rbind, x)})




##Paste Enhancers for each promoter 
dataFrame2<-NULL
for(i in 1:dim(PromoterIds)[1]){
  test<-PromoterIds[i,]
  
  dataFrame1<-NULL
  for(t in 1:(test$PromoterIds%>%as.data.frame%>%dim)[2]){
    dataFrame1<-rbind(dataFrame1,test[2:4])
  }
  
  dataFrame2<-rbind.data.frame(dataFrame2,
                               cbind.data.frame("Promoter Ids" =t(test$PromoterIds%>%as.data.frame()),
                                                dataFrame1)
  )
}
dataFrame2$`Promoter Ids`<-dataFrame2$`Promoter Ids`%>%as.character()





#Convert from UCSC to gene symbol because biomart wasn't working/ was only returning maybe 5-10 genes out of the 50!

##Bedfile From ucsc with all the metadata columns ;-) 
UCSCConvter <- read_delim("~/DataFiles/Gene Tracks/Human/hg19WithNames.bed", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

#We can left join based on the transcript Ids from the promoter id's gene
GeneSymbolTHesisTable<-left_join(dataFrame2, UCSCConvter, by= c("Promoter Ids" = "hg19.kgXref.kgID" ))


#Now we can make a table
EnhancerPromoterInteractionScoreCounts<-cbind.data.frame("Gene Symbol"=GeneSymbolTHesisTable$hg19.kgXref.geneSymbol,
                                                         "UCSC Transcript Id"= GeneSymbolTHesisTable$`Promoter Ids`,
                                                         "Chromosome"= GeneSymbolTHesisTable$hg19.knownGene.chrom,
                                                         "Start" = GeneSymbolTHesisTable$hg19.knownGene.txStart,
                                                         "End" =GeneSymbolTHesisTable$hg19.knownGene.txEnd,
                                                         "Enhancer"=GeneSymbolTHesisTable$EnhancerIds, 
                                                         "Q_value"=GeneSymbolTHesisTable$q_value)

#Subset for unique Ids
UniqueInteractionsBetweenHumansAndGeneSYmbo<-EnhancerPromoterInteractionScoreCounts[isUnique(EnhancerPromoterInteractionScoreCounts$`Gene Symbol`),]
#TABLE

write.table(UniqueInteractionsBetweenHumansAndGeneSYmbo, file = "~/Thesis/uniquePromoterEnhancersHuman",
            append=FALSE,
            quote=FALSE, 
            col.names=TRUE,
            row.names = FALSE,
            sep = "\t")

#############################################################
#########Identifying Unique Mouse Interactions BAM BAM MAMAN.
###############################################################



library(GenomicInteractions)
library(magrittr)
GiTrack<-readRDS("~/DataFiles/HiC/Mouse/mm9StasticallySignificantInteractions")
dataFrameGiTrack<-GiTrack%>%as.data.frame()



Spilt1<-subset(dataFrameGiTrack, !enhancer.id2=="NA")
Spilt2<-subset(dataFrameGiTrack, !enhancer.id1=="NA")

Spilt1<-cbind.data.frame("PromoterIds"=Spilt1$promoters.id1,
                         "EnhancerIds"= Spilt1$enhancer.id2,
                         "Counts" = Spilt1$counts, 
                         "q_value" = Spilt1$q_value)    

Split2<-cbind.data.frame("PromoterIds"=Spilt2$promoters.id2,
                          "EnhancerIds"= Spilt2$enhancer.id1,
                          "Counts" = Spilt2$counts, 
                          "q_value" = Spilt2$q_value)    
PromoterIds<-rbind(Spilt1,Split2)


PromoterIds$PromoterIds<-lapply(PromoterIds$PromoterIds, list)


PromoterIds$PromoterIds<-lapply(PromoterIds$PromoterIds,function(x){do.call(rbind, x)})



dataFrame2<-NULL
for(i in 1:dim(PromoterIds)[1]){
  test<-PromoterIds[i,]
  
  dataFrame1<-NULL
  for(t in 1:(test$PromoterIds%>%as.data.frame%>%dim)[2]){
    dataFrame1<-rbind(dataFrame1,test[2:4])
  }
  
  dataFrame2<-rbind.data.frame(dataFrame2,
                               cbind.data.frame("Promoter Ids" =t(test$PromoterIds%>%as.data.frame()),
                                                dataFrame1)
  )
}
dataFrame2$`Promoter Ids`<-dataFrame2$`Promoter Ids`%>%as.character()




library(readr)
library(dplyr)
UCSCConvter <- read_delim("~/DataFiles/Gene Tracks/Mouse/mm9withnames.bed", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)


GeneSymbolTHesisTable<-left_join(dataFrame2, UCSCConvter, by= c("Promoter Ids" = "mm9.kgXref.kgID" ))

EnhancerPromoterInteractionScoreCountsMouse<-cbind.data.frame("Gene Symbol"=GeneSymbolTHesisTable$mm9.kgXref.geneSymbol,
                                                         "UCSC Transcript Id"= GeneSymbolTHesisTable$`Promoter Ids`,
                                                         "Chromosome"= GeneSymbolTHesisTable$mm9.knownGene.chrom,
                                                         "Start" = GeneSymbolTHesisTable$mm9.knownGene.txStart,
                                                         "End" =GeneSymbolTHesisTable$mm9.knownGene.txEnd,
                                                         "Enhancer"=GeneSymbolTHesisTable$EnhancerIds, 
                                                         "Q_value"=GeneSymbolTHesisTable$q_value)

UniqueInteractionsBetweenMouseAndGeneSYmbo<-EnhancerPromoterInteractionScoreCountsMouse[isUnique(EnhancerPromoterInteractionScoreCountsMouse$`Gene Symbol`),]


write.table(UniqueInteractionsBetweenMouseAndGeneSYmbo, file = "~/Thesis/uniquePromoterEnhancerMouse",
            append=FALSE,
            quote=FALSE, 
            col.names=TRUE,
            row.names = FALSE,
            sep = "\t")


