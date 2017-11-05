##Table developer

##First thing is that we need to convert all relevant data to dataframes for easy queries




library(readr)
library(dplyr)
#importing the relevant Data

geneSymbols <- read_tsv("~/DataFiles/Tessa differential Expressed/2colmultiplenames.txt")
pointer<-import(con= "~/DataFiles/Tessa differential Expressed/PA16col.bed")%>%as.data.frame()
enhancerDataFrame<- import( con = "~/DataFiles/Enhancer Tracks/Mouse/Enhanceresmm9.bed",  which =
                              GRanges(chrM, IRanges(1, end = .Machine$integer.max - 1)))%>%ranges()%>%as.data.frame()
dataFrameConversation<- dataFrameGrange
differentiallyExpressed <- left_join(pointer, geneSymbols, by = c("name"="#kgID"))
geneNames<- as.data.frame(unique(geneSymbols[2]))#subsetting for intergert
x<-1#intializing X

for(i in 1:dim(geneNames[1])){


selectedGene<-differentiallyExpressed[differentiallyExpressed$geneSymbol==as.character(geneNames[x,1]),]#this will select for each gene ID from the differentially expressed data
toM<-max(selectedGene[3])+10000##max to get the longest transcript distance
fromM<-min(selectedGene[2])-10000##to get the max utilising this as a promoter 
chrM<-selectedGene[1,1]##plotting the chromosome! this will be fucked need to revalute the chromosme


##closet ARX motif, done via subtraction of gene start site - ARX data frame and locating the closest TFBS

##distance to ARX motif, Select for motif starts and the start site of gene find the lowest value
tableDistanceToArxTFBS<-  min(abs(mouseArxDataFrame[2]-min(selectedGene[2])))## utlising min and abs to get the closest DataFrame distnace

#remove the NA from the Cluster Data frame ##i shoudl rename this ;') 
##find distance between the start site of gene and the lo
grangeDataFrame3<-grangeDataFrame2[-1,]
tableDistanceToArxCluster<-min(abs(grangeDataFrame3[1]-min(selectedGene[2])),
                               na.rm = FALSE)
##distance of motifs to enhancer?
tableDistanceToEnhancer<- min(abs(enhancerDataFrame-(min(selectedGene[2])+tableDistanceToArxTFBS)))
                                  
## okay so this is gonna be confusing but simply                             
tableConservationScoreOfTFBS<-dataFrameConversation$score[abs(min(selectedGene[2]+tableDistanceToArxTFBS)-dataFrameConversation[2])==min(abs(min(selectedGene[2]+tableDistanceToArxTFBS)-dataFrameConversation[2]))]
tableConservationScoreOfCluster<-dataFrameConversation$score[abs(min(selectedGene[2]+tableDistanceToArxCluster)-dataFrameConversation[2])==min(abs(min(selectedGene[2]+tableDistanceToArxCluster)-dataFrameConversation[2]))]

if(!exists("ArxTable3"))##if the ArxTable3 doesn't exist it will create it 
  {

ArxTable3<- cbind("Gene"=selectedGene[1,8],
             "Chromosome" = as.character(selectedGene[1,1]),
            "Start"= min(selectedGene[2]),
             "End"=max(selectedGene[3]), 
             "Distance to closest Arx TFBS"=tableDistanceToArxTFBS,
            "Conservation score of TFBS"= tableConservationScoreOfTFBS,
             "Distance to closest cluster"=tableDistanceToArxCluster,
             "Conservation score of ARX cluster"=tableConservationScoreOfCluster)
} else {## if it does exist it will excute and add theac line usinig R line
  ArxTable3<- rbind(ArxTable3,cbind(
                   "Gene"=selectedGene[1,8],
                   "Chromosome" = as.character(selectedGene[1,1]),
                   "Start"= min(selectedGene[2]),
                   "End"=max(selectedGene[3]), 
                   "Distance to closest Arx TFBS"=tableDistanceToArxTFBS,
                   "Conservation score of TFBS"= tableConservationScoreOfTFBS,
                   "Distance to closest cluster"=tableDistanceToArxCluster,
                   "Conservation score of ARX cluster"=tableConservationScoreOfCluster))

}
x<<- x+1
print(x) 
}

