## Jaspar 

library(JASPAR2016)
library(magrittr)
library(TFBSTools)
library(RSQLite)
mouse<-list()
mouse[["species"]]<-"Mus musculus"
test<-getMatrixSet(JASPAR2016,mouse)
matrixGetter<-function(x){x@profileMatrix}
test<-lapply(test, matrixGetter)
test<-lapply(test, t)
i<-1
for(i in 1:length(test)){
  colnames(test[i])<-NULL  
  i<-i+1
}


sink("JASPAR musMusculus matrices")
test
sink()
