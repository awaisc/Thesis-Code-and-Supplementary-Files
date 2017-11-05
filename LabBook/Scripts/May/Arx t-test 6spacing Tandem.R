
library(BSgenome.Mmusculus.UCSC.mm9)
library(Biostrings)
library(magrittr)
library(parallel)

##generating random 6mers
every5Mer <- unique(DNAStringSet(
  sapply(
    sample(c(120), 120),
    function(size)
      paste(sample(DNA_BASES, 5, replace=TRUE), collapse="")
  )
))


##genearing the the spacing of the 6mers

all5merTandemMatrixs<-NULL
a<-1

for(a in 1:length(every6Mer)){
  all5merTandemMatrixs[[a]]<- cbind(round(PWM(every6Mer[a,])*7),
                                    0.25,
                                    round(PWM(every6Mer[a,])*7))
  a<-a+1
}

##COnverting to matrix beucase matchPWM only takes matrices
all5merTandemMatrixs<-lapply(all5merTandemMatrixs, as.matrix)


##matchPWM function
matchPWMFunction<-function(x){
  matchPWM(x, BSgenome.Mmusculus.UCSC.mm9, "100")
}
i<-1
countsOf5MerAt1Space<-NULL
for(i in 1:length(all5merTandemMatrixs)){
  countsOf5MerAt1Space[i]<-length(matchPWM(all5merTandemMatrixs[[i]], BSgenome.Mmusculus.UCSC.mm9, "100%"))
  i<-i+1
}
t.test(countsOf6MerAt6Spaces)
