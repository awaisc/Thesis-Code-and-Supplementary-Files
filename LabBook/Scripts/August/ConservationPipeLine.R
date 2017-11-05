

genome<-BSgenome.Mmusculus.UCSC.mm9
arx6merTFBS<-matchPWM(round(PWM("TAATTA")*7), genome, "100%")




enhancers<-import("~/DataFiles/Enhancer Tracks/Mouse/mouse_permissive_enhancers_phase_1_and_2.bed")
genes<-import("~/DataFiles/Gene Tracks/Mouse/mm9.bed")
promoters<-promoters(genes)


MotifModel<-arx6merTFBS

#Converting ARX binding sites to the mm10 genome for getting thhe scores

ah<-AnnotationHub()
mm10Chain<-AnnotationHub::query(ah, "mm9ToMm10.over.chain.gz")
mm10ChainFile<-mm10Chain[["AH14596"]]

mm9ChainFile<-AnnotationHub::query(ah, "mm10ToMm9.over.chain.gz")
mm9Chainfile<-mm9ChainFile[["AH14535"]]

motifModelmm10<-liftOver(chain = mm10ChainFile, x= MotifModel)%>%unlist



conservationPhastCon<-import("~/DataFiles/Conservation/Mouse/mm10.60way.phastCons.bw", which= motifModelmm10)

conservationPhloPy<-import("~/DataFiles/Conservation/Mouse/mm10.60way.phyloP60way.bw",which =motifModelmm10)


##Lifting Over the mm10 Conservation Track
#Conservation Of All Binidng sites Of ARX!
PhastConMouseLifted<-liftOver(chain = mm9Chainfile, x = conservationPhastCon)%>%unlist()

PhlyoPMouseLifted<-liftOver(chain = mm9Chainfile, x = conservationPhloPy)%>%unlist()

##Subset for only the most conserved regions
#we do a 1.3 instead of 0.9 as phylop is scores between -14 and 3 hence, 90% conservation is 15.3
PhlyoPMouseLiftedConservedRegions<-subset(PhlyoPMouseLifted, score>=1.3)
PhastConMouseLiftedConservedRegions<-subset(PhastConMouseLifted, score >=0.9)
##See which Motifs fall into these regions
polyPConserved<-subset(arx6merTFBS,countOverlaps( arx6merTFBS,PhlyoPMouseLiftedConservedRegions)>=6)

##See which Motifs fall into these regions
phastConserved<-subset(arx6merTFBS,countOverlaps(arx6merTFBS, PhastConMouseLiftedConservedRegions)>=6)



genomicLocation<-function(x){
  dataFrame1<-rbind.data.frame(
    "Enhancers" =subsetByOverlaps(x, enhancers)%>%length(),
    
    "Genes" = subsetByOverlaps(x, genes)%>%length(),
    
    "Promoters"= subsetByOverlaps(x, promoters)%>%length(),
    
    "Other"= length(x)-sum(
      "Enhancers" =subsetByOverlaps(x, enhancers)%>%length(),
      "Genes" = subsetByOverlaps(x, genes)%>%length(),
      "Promoters"= subsetByOverlaps(x, promoters)%>%length()))
  
  colnames(dataFrame1) <-c("Number Of Motifs")
  return(dataFrame1)
}

ARXMotifs<-cbind.data.frame(genomicLocation(phastConserved), 
                            genomicLocation(polyPConserved),
                            genomicLocation(arx6merTFBS))%>%rownames_to_column()
colnames(ARXMotifs)<-c("Genomic Location",
                       "PhastCon",
                       "PhyloP", 
                       "All")

ARXMotifs%>%pander()

reshaped<-ARXMotifs%>%melt()

ggplot(data = reshaped, mapping = aes(x=`variable`,  y= `value`, fill= `Genomic Location`))+
  geom_bar(stat="identity")+
  theme_bw()+
  ylab(label = "Number of Motifs")+
  xlab(label= "Conservation Method")+
  guides(fill=guide_legend(title= "Genomic Location"))
