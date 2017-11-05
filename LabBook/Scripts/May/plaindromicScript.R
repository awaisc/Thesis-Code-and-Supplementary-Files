

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



all6merPlaindromicMatrixs<-NULL
a<-1



for(a in 1:length(every6Mer)){
  all6merPlaindromicMatrixs[[a]]<- cbind(round(PWM(every6Mer[a,])*7),
                                    0.25, 0.25,0.25, 0.25, 
                                    round(PWM(reverseComplement(complement(every6Mer[a,])))*7))
  a<-a+1
}


grange6MersPLaindrome4Space<-NULL
i<-1
for(i in 1:length(all6merTandemMatrixs)){
  grange6MersPLaindrome4Space[i]<-length(matchPWM(all6merPlaindromicMatrixs[[i]], BSgenome.Mmusculus.UCSC.mm9, "100%"))
  i<-i+1
}

test<-function(x){matchPWM(x, BSgenome.Mmusculus.UCSC.mm9, "100%")}

mclapply(all6merPlaindromicMatrixs, FUN = test, mc.cores = 4)
t.test(grange6MersPLaindrome4Space)

save.image("~/Scripts/May/Arx6mers2.RData")



