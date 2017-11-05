library(BSgenome.Hsapiens.UCSC.hg19)
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

all6merTandemMatrixs<-NULL
a<-1

for(a in 1:length(every6Mer)){
  all6merTandemMatrixs[[a]]<- cbind(round(PWM(every6Mer[a,])*7),
                                    0.25, 0.25,
                                    round(PWM(every6Mer[a,])*7))
  a<-a+1
}

##COnverting to matrix beucase matchPWM only takes matrices
all6merTandemMatrixs<-lapply(all6merTandemMatrixs, as.matrix)


##matchPWM function
matchPWMFunction<-function(x){
  matchPWM(x, BSgenome.Hsapiens.UCSC.hg19, "100")
}
i<-1
countsOf6MerAt2Space<-NULL
for(i in 1:length(all6merTandemMatrixs)){
  countsOf6MerAt2Space[i]<-length(matchPWM(all6merTandemMatrixs[[i]], BSgenome.Hsapiens.UCSC.hg19, "100%"))
  i<-i+1
}



t.test(countsOf6MerAt2Space)
