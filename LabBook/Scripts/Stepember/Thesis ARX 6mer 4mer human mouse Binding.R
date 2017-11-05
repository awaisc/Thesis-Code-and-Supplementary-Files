arx6merTFBS




ARX6Human<-readRDS("~/DataFiles/ChIPseq/Human/ARX6merHg19Sites")

arx6mer<-readRDS("~/DataFiles/ChIPseq/Mouse/ARX6mermm9Sites")
enhancers<-import("~/DataFiles/Enhancer Tracks/Mouse/mouse_permissive_enhancers_phase_1_and_2.bed")
promoters<-import("~/DataFiles/Gene Tracks/Mouse/mm9.bed")%>%promoters()
genes<-import("~/DataFiles/Gene Tracks/Mouse/mm9.bed")

cbind.data.frame(
"Enhancers"=subsetByOverlaps(arx6mer, enhancers)%>%length(),
"Genes"=subsetByOverlaps(arx6mer, genes)%>%length(),
"Promoters"=subsetByOverlaps(arx6mer, promoters)%>%length(),
"Other"=length(arx6mer)-sum(subsetByOverlaps(arx6mer, enhancers)%>%length(),
                    subsetByOverlaps(arx6mer, genes)%>%length(),
                    subsetByOverlaps(arx6mer, promoters)%>%length())
)
arx6merTFBS




ARX6Human<-readRDS("~/DataFiles/ChIPseq/Human/ARX6merHg19Sites")

arx6mer<-readRDS("~/DataFiles/ChIPseq/Mouse/ARX6mermm9Sites")
enhancers<-import("~/DataFiles/Enhancer Tracks/Mouse/mouse_permissive_enhancers_phase_1_and_2.bed")
promoters<-import("~/DataFiles/Gene Tracks/Mouse/mm9.bed")%>%promoters()
genes<-import("~/DataFiles/Gene Tracks/Mouse/mm9.bed")

cbind.data.frame(
  "Enhancers"=subsetByOverlaps(arx6mer, enhancers)%>%length(),
  "Genes"=subsetByOverlaps(arx6mer, genes)%>%length(),
  "Promoters"=subsetByOverlaps(arx6mer, promoters)%>%length(),
  "Other"=length(arx6mer)-sum(subsetByOverlaps(arx6mer, enhancers)%>%length(),
                              subsetByOverlaps(arx6mer, genes)%>%length(),
                              subsetByOverlaps(arx6mer, promoters)%>%length())
)
