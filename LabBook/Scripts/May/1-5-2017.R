library(BSgenome.Mmusculus.UCSC.mm9)
library(Biostrings)
library(magrittr)
library(parallel)

ArxPlaindrmicMinus1<-rbind( A=c(0,1,1,0,0,1,1,0,0,1), 
                           C=c(0,0,0,0,0,0,0,0,0,0),
                           G=c(0,0,0,0,0,0,0,0,0,0),
                           T=c(1,0,0,1,1,0,0,1,1,0))
ArxPlaindromicNoSpace<-rbind( A=c(0,1,1,0,0,1,0,1,1,0,0,1), 
                              C=c(0,0,0,0,0,0,0),
                              G=c(0,0,0,0,0,0,0) ,
                              T=c(1,0,0,1,1,0,1,0,0,1,1,0))

arxTandemNoSpace<-rbind( A=c(0,1,1,0,0,1,0,1,1,0,0,1),
                        C=c(0,0,0,0,0,0,0,0,0,0,0,0),
                        G=c(0,0,0,0,0,0,0,0,0,0,0,0) ,
                        T=c(1,0,0,1,1,0,1,0,0,1,1,0))

ArxTandemMinus1<-rbind( A=c(0,1,1,0,0,1,1,0,0,1),
                        C=c(0,0,0,0,0,0,0),
                        G=c(0,0,0,0,0,0,0) ,
                        T=c(1,0,0,1,1,0,0,1,1,0))

grangeArxMinus1<-matchPWM(ArxPlaindrmicMinus1, BSgenome.Mmusculus.UCSC.mm9, "100%")
grangeArxNoSpaceTandem<-matchPWM(arxTandemNoSpace, BSgenome.Mmusculus.UCSC.mm9, "100%")
grangeArxNoSpacePlaindrome<-matchPWM(ArxPlaindromicNoSpace, BSgenome.Mmusculus.UCSC.mm9, "100%")




every6Mer <- unique(DNAStringSet(
  sapply(
    sample(c(120), 120),
    function(size)
      paste(sample(DNA_BASES, 6, replace=TRUE), collapse="")
  )
))



all6merTandemMatrixs<-NULL
a<-1

for(a in 1:length(every6Mer)){
all6merTandemMatrixs<-rbind(all6merTandemMatrixs ,
                             cbind(PWM(every6Mer[a,]), 0.25, 0.25, PWM(every6Mer[a,])))
a<-a+1
}


numberOf6Mers<-NULL
k<-1
for(i in 1:dim(all6merTandemMatrixs)[1]){
  grange6MersTandem2Space<-matchPWM(all6merTandemMatrixs[k:(k+3),], BSgenome.Mmusculus.UCSC.mm9, "100%")
  numberOf6Mers<-rbind(numberOf6Mers, length(grange6MersTandem2Space))
  k<-k+4
}


function(x){
matchPWM(all6merTandemMatrixs[i:i+3,], BSgenome.Mmusculus.UCSC.mm9, "100%")
}

