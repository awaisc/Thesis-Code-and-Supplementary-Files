##generating a position weight matrix for every time
library(BSgenome.Mmusculus.UCSC.mm9)
library(Biostrings)
library(magrittr)
library(parallel)


##selecting for 4 random 6mers
every6Mer <- unique(DNAStringSet(
  sapply(
    sample(c(120), 120),
    function(size)
      paste(sample(DNA_BASES, 6, replace=TRUE), collapse="")
  )
))
## Putting the 6mers into a tandem orientation DataFrame
every6Mer<-every6Mer%>%as.matrix()
MerTandem2Space<-cbind(every6Mer, "..", every6Mer)



##generating astrings from the tandem orientation dataFrame
i<-1
listOf6Mers<-NULL

for(i in 1:dim(MerTandem2Space)[1]){
 listOf6Mers<- rbind(listOf6Mers, paste0(MerTandem2Space[i,1], MerTandem2Space[i,3]))
 i<-i+1
}



##Generating the 6mers in the palaindormic Arrangement!
##Getting the plaindromic sequence

plaindromic6mer<-every6Mer%>%as.matrix%>%DNAStringSet()%>%reverseComplement()%>%as.data.frame

#Putting it into a dataFrame step
MerPlaindrome2Space<-cbind(every6Mer, "NN", plaindromic6mer)

##Getting the DNA strings 
i<-1
listOf6MersPlaindrome<-NULL
for(i in 1:dim(MerPlaindrome2Space)[1]){
  listOf6MersPlaindrome<- rbind(listOf6MersPlaindrome, paste0(MerPlaindrome2Space[i,1], "NN", MerPlaindrome2Space[i,3]))
  i<-i+1
}

##combining it all into a single list, a DNAstringSet wont work with vmatch. 
##the first one is Arx 6mer with 2 spaces.
all2Space<-rbind("TAATTANNATTAAT",listOf6Mers, listOf6MersPlaindrome)%>%DNAStringSet()

##Setting the parameters for the function

Storage<-function(x){
  #Note 2 mistmatches in vmatch.
Matches<-vmatchPattern(x, getSeq(BSgenome.Mmusculus.UCSC.mm9), max.mismatch = 2 )%>%as.data.frame()

##Renaming the list to the relevant chromosomes
Matches[Matches=="1"]<-"chr1"
Matches[Matches=="2"]<-"chr2"
Matches[Matches=="3"]<-"chr3"
Matches[Matches=="4"]<-"chr4"
Matches[Matches=="5"]<-"chr5"
Matches[Matches=="6"]<-"chr6"
Matches[Matches=="7"]<-"chr7"
Matches[Matches=="8"]<-"chr8"
Matches[Matches=="9"]<-"chr9"
Matches[Matches=="10"]<-"chr10"
Matches[Matches=="11"]<-"chr11"
Matches[Matches=="12"]<-"chr12"
Matches[Matches=="13"]<-"chr13"
Matches[Matches=="14"]<-"chr14"
Matches[Matches=="15"]<-"chr15"
Matches[Matches=="16"]<-"chr16"
Matches[Matches=="17"]<-"chr17"
Matches[Matches=="18"]<-"chr18"
Matches[Matches=="19"]<-"chr19"
Matches[Matches=="20"]<-"chrX"
Matches[Matches=="21"]<-"chrY"
Matches[Matches=="22"]<-"chrM"
Matches[Matches=="23"]<-"chr1_random"
Matches[Matches=="24"]<-"chr3_random"
Matches[Matches=="25"]<-"chr4_random"
Matches[Matches=="26"]<-"chr5_random"
Matches[Matches=="27"]<-"chr7_random"
Matches[Matches=="28"]<-"chr8_random"
Matches[Matches=="29"]<-"chr9_random"
Matches[Matches=="30"]<-"chr13_random"
Matches[Matches=="31"]<-"chr16_random"
Matches[Matches=="32"]<-"chr17_random"
Matches[Matches=="33"]<-"chrX_random"
Matches[Matches=="34"]<-"chrY_random"
Matches[Matches=="35"]<-"chrU_random"
return(makeGRangesFromDataFrame(Matches,
                         keep.extra.columns=FALSE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid", "group"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE))
}

##Trying this with a forloop
grangeAll2Sapces<-NULL
for(i in 1:length(all2Space)){
  grangeAll2Sapces<-cbind(grangeAll2Sapces,Storage(all2Space[i]))
  i<-i+1
  
}

##Trying this as a lapply function

grangeAll2Sapces<-mclapply(all2Space, Storage, mc.cores=7)




##fullforloop
j<-1
numberof2spaced<-NULL

for(j in 1:length(all2Space))
  #Note 2 mistmatches in vmatch.
  Matches<-vmatchPattern(all2Space[j], getSeq(BSgenome.Mmusculus.UCSC.mm9), max.mismatch = 2 )%>%as.data.frame()
  
  ##Renaming the list to the relevant chromosomes
  Matches[Matches=="1"]<-"chr1"
  Matches[Matches=="2"]<-"chr2"
  Matches[Matches=="3"]<-"chr3"
  Matches[Matches=="4"]<-"chr4"
  Matches[Matches=="5"]<-"chr5"
  Matches[Matches=="6"]<-"chr6"
  Matches[Matches=="7"]<-"chr7"
  Matches[Matches=="8"]<-"chr8"
  Matches[Matches=="9"]<-"chr9"
  Matches[Matches=="10"]<-"chr10"
  Matches[Matches=="11"]<-"chr11"
  Matches[Matches=="12"]<-"chr12"
  Matches[Matches=="13"]<-"chr13"
  Matches[Matches=="14"]<-"chr14"
  Matches[Matches=="15"]<-"chr15"
  Matches[Matches=="16"]<-"chr16"
  Matches[Matches=="17"]<-"chr17"
  Matches[Matches=="18"]<-"chr18"
  Matches[Matches=="19"]<-"chr19"
  Matches[Matches=="20"]<-"chrX"
  Matches[Matches=="21"]<-"chrY"
  Matches[Matches=="22"]<-"chrM"
  Matches[Matches=="23"]<-"chr1_random"
  Matches[Matches=="24"]<-"chr3_random"
  Matches[Matches=="25"]<-"chr4_random"
  Matches[Matches=="26"]<-"chr5_random"
  Matches[Matches=="27"]<-"chr7_random"
  Matches[Matches=="28"]<-"chr8_random"
  Matches[Matches=="29"]<-"chr9_random"
  Matches[Matches=="30"]<-"chr13_random"
  Matches[Matches=="31"]<-"chr16_random"
  Matches[Matches=="32"]<-"chr17_random"
  Matches[Matches=="33"]<-"chrX_random"
  Matches[Matches=="34"]<-"chrY_random"
  Matches[Matches=="35"]<-"chrU_random"
  numberof2spaced<-rbind(numberof2spaced, length(makeGRangesFromDataFrame(Matches,
                                                           keep.extra.columns=FALSE,
                                                           ignore.strand=TRUE,
                                                           seqinfo=NULL,
                                                           seqnames.field=c("seqnames", "seqname",
                                                                            "chromosome", "chrom",
                                                                            "chr", "chromosome_name",
                                                                            "seqid", "group"),
                                                           start.field="start",
                                                           end.field=c("end", "stop"),
                                                           strand.field="strand",
                                                           starts.in.df.are.0based=FALSE))
  )
  
j<-j+1





## starting again


every6Mer <- unique(DNAStringSet(
  sapply(
    sample(c(20), 20),
    function(size)
      paste(sample(DNA_BASES, 6, replace=TRUE), collapse="")
  )
))







###StackOverflow

matchPWM(PWM(all2Space[1], type = c("log2probratio", "prob"),
    prior.params = c(A=0.25, C=0.25, G=0.25, T=0.25)), BSgenome.Mmusculus.UCSC.mm9, "100%")


pwm <- all2Space[DNA_BASES, ]
matchPWM(pwm, subject)



#Create a matrix where 1 can be in all four positions (in 4 different columns)
m1 = sapply(1:12, function(i)
  replace(x = rep(0, 4), list = i, values = 1))

#Create a vector for the third column   
v1 = rep(x = 0.25, times = 4)

#Get all combinations of indices between 1 and 4 (with repetition)
x = 1:12
indices = data.frame(t(expand.grid(x, x, x)))

#Get the matrices
lapply(X = indices,
       function(a) cbind(m1[,a[1:2]],v1,m1[,a[3]]))
