
library(BSgenome.Mmusculus.UCSC.mm9)
library(Biostrings)
library(magrittr)
library(parallel)

##generating random 6mers
every6Mer <- unique(DNAStringSet(
  sapply(
    sample(c(120), 120),
    function(size)
      paste(sample(DNA_BASES, 6, replace=TRUE), collapse="")
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

conversion6to5Mer<-function(x){cbind(x[,1:5], x[,7], x[,9:13])}
  


##Converting the 6mers to 5mers with the last base lobbed off and the starting based off
all5merTandemMatrixs<-lapply(all5merTandemMatrixs, conversion6to5Mer)

##COnverting to matrix beucase matchPWM only takes matrices
all5merTandemMatrixs<-lapply(all5merTandemMatrixs, as.matrix)


##matchPWM function

i<-1
countsOf5MerAt1Space<-NULL
for(i in 1:length(all5merTandemMatrixs)){
  countsOf5MerAt1Space[i]<-length(matchPWM(all5merTandemMatrixs[[i]], BSgenome.Mmusculus.UCSC.mm9, "100%"))
  i<-i+1
}
t.test(countsOf5MerAt1Space, mu = 22346)
