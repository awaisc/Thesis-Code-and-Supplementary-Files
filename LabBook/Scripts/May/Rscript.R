##

library(BSgenome.Mmusculus.UCSC.mm9)
library(Biostrings)
library(magrittr)
library(parallel)


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

t.test(numberOf6Mers)


all6merPlaindromicMatrixs<-NULL
a<-1

for(a in 1:length(every6Mer)){
  all6merPlaindromicMatrixs<-rbind(all6merPlaindromicMatrixs ,
                              cbind(PWM(every6Mer[a,]), 0.25, 0.25,0.25, 0.25, 
                                    PWM(reverseComplement(complement(every6Mer[a,])))))
  a<-a+1
}


numberOf6MersPlaindromine<-NULL
k<-1
for(i in 1:length(every6Mer)){
  grange6MersPLaindrome4Space<-matchPWM(mer6Plaindromine[k:(k+3),], 
                                        BSgenome.Mmusculus.UCSC.mm9,
                                        "100%")
  numberOf6MersPlaindromine<-rbind(numberOf6MersPlaindromine, 
                                   length(grange6MersPLaindrome4Space))
  k<-k+4
}

t.test(numberOf6MersPlaindromine)

save.image("~/Scripts/May/Arx6mers1.RData")

