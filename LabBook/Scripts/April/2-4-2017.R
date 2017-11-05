### Combinding into new  Table
library(readxl)
tessaData <- read_excel("~/Research Proposal/Supplementary_tables_1-3.xlsx", 
                        skip = 1)%>% as.data.frame()

ArxTable3DataFrame<- as.data.frame(ArxTable3)
rnaSeqData <- left_join(ArxTable3DataFrame, tessaData, by = c("Gene"="Gene ID"))
newTable<-cbind(rnaSeqData, "absoulte of log FC"= abs(rnaSeqData$logFC) )

##plotting the new Table and seeing that there is not correlation!
library(ggplot2)
ggplot(data= newTable, aes(x="Distance to closest Arx TFBS", y="absoulte of log FC" ) )+
  geom_jitter()+
  geom_point()+
  theme_bw()
#no relationship between the distance and the differential expression!