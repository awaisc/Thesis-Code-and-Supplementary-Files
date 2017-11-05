## Trying to cocer the WEEDER output into R and then to a readable format by STAMP. 

temp = list.files(pattern="file.txt")

delim<-function(x){ read_delim(x, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 1, na = "NA")}
roundFunction<-function(x){ round(x[,1:length(x)-1]*100)}
colnamesFunction<-function(x){colnames(x)<-NULL}
binder<-function(x){rbind(">ABC", x)}
myfiles = lapply(temp, delim)
myfiles1= lapply(myfiles, na.omit)

myfiles=lapply(myfiles, as.numeric)
myfiles= lapply(myfiles1, roundFunction)
myfiles=lapply(myfiles, t)

names(myfiles)=c("DE 1",
                 "DE 2",
                 "DE 3",
                 "DE 4",
                 "DE 5",
                 "DE 6",
                 "DE 7",
                 "DE 8",
                 "DE 9"
                 )
write.table(myfiles, row.names = FALSE, col.names = FALSE )


lapply(myfiles, write.table, row.names=FALSE, col.names= FALSE)%>%write.csv()

