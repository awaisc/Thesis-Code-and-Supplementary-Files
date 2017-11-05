## This script requires a script bash before and bash script after

#!bin/bash
awk '/>/{n++}{print >"common" n ".txt" }' CommonStrings.matrix.w2 

R
library(readr)


##Import Data into R
read.csvFunction<- function(x){ read_delim(x, "\t", escape_double = FALSE, na = "NA", trim_ws = TRUE, skip = 1,col_names = FALSE)}

##Select for type of File i want

files <- list.files(path="~/DataFiles/ChIPseq//",
                         pattern="common")
##WD
setwd("~/DataFiles/ChIPseq/")

##Import now
PWM<-lapply(files, read.csvFunction)

##Covert to Numeric
asNumericFunction<-function(x){apply(x, 2, as.numeric)}
PWM<-lapply(PWM, asNumericFunction)

##Convert to PFM
roundFunction<-function(x){round(100*x)}
PFM<-lapply(PWM, roundFunction)

##ADD header
for(i in 1:length(PFM)){
  names(PFM)[[i]] <- paste0(">",i)
  i<-i+1
}
##Removve the Column of NAs
subsetter<-function(x){x[,2:(dim(x)[2]-1)]}
PFM<-lapply(PFM, subsetter)


##Transmutate to fit format
PFM<-lapply(PFM, t)
PFM<- lapply(PFM, as.data.frame)

##Remove colnames
i<-1
for(i in 1:length(PFM)){
  colnames(PFM[[i]])<-NULL  
i<-i+1
}

##Remove rownames <- doesnt work because R keeps them 
i<-1
for(i in 1:length(PFM)){
  rownames(PFM[[i]])<-NULL
  i<-i+1
}


## exporting the data mainting the Structure

sink("output.txt")
PFM
sink()

## Removing the First column of Row names
#!bin/bash
cat output.txt | awk '{print $2 "\t" $3 "\t" $4"\t" $5 "\t" $6 "\t" $7}'> commonMatrices.txt
sed