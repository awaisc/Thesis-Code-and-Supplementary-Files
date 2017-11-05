library(Gviz)
library(biomaRt)
EnhancerConservedMotifs<-readRDS(file = "~/DataFiles/Enhancer Tracks/Human/ARX Motif Enhancers")
fantom5Enhancers<- import("~/DataFiles/Enhancer Tracks/Human/human_permissive_enhancers_phase_1_and_2.bed")
imrCellTADS<-import("~/DataFiles/HiC/Human/IMR90_domains_hg19.bed")
hg19Genes<-import("~/DataFiles/Gene Tracks/Human/hg.bed")
hg19Promoters<-promoters(hg19Genes)

ARXEnrichedPromoters<-subset(hg19Promoters,countOverlaps(hg19Promoters, arx6merTFBS)>=4)



ARXInPromotersEnriched<-subsetByOverlaps(arx6merTFBS,ARXEnrichedPromoters)

ARXImage<-AnnotationTrack(ARXInPromotersEnriched,name ="ARX Motifs" )
subset(ARXInPromotersEnriched, countOverlaps(ARXInPromotersEnriched, phyloPTrack)>=4)



phylopARXScores<-import("~/DataFiles/Conservation/Human/hg19.100way.phyloP100way.bw",which =ARXEnrichedPromoters)



##Subset for only the most conserved regions
#we do a 1.3 instead of 0.9 as phylop is scores between -14 and 3 hence, 90% conservation is 15.3
phyloPTrack<-subset(phylopARXScores, score>=1.3)

# ##See which Motifs fall into these regions
# polyPConserved<-subset(arx6merTFBS,countOverlaps( arx6merTFBS,phyloPTrack)>=6)

##Phast Con Scores Same as above
phastConARXScores<-import("~/DataFiles/Conservation/Human/hg19.100way.phastCons.bw",which =ARXEnrichedPromoters)

PhastConTrack<-subset(phastConARXScores, score>=0.7)

##See which Motifs fall into these regions
phastConserved<-subset(arx6merTFBS,countOverlaps(arx6merTFBS, PhastConTrack)>=6)  


backgroundColor<-"#8C867A"


ARXImage<-AnnotationTrack(ARXInPromotersEnriched,
                          name ="ARX Motifs" ,     
                          fill="#3FB8AF", 
                          background.track = "#FFFEDB", 
                          background.title= backgroundColor, 
                          fontcolor="black", 
                          cex.title = 0.72, 
                          rotation.title = 0, 
                          lwd.title = 0.5, title.width = 1.5, 
                          cex.main = 1,
                          fontcolor.title = "white",
                          color.scale= "black")
ARXEnrichedPromotersTrack<-AnnotationTrack(ARXEnrichedPromoters, name= "Promoters", 
                                           fill="#FF3D7F",
                                           col.baseline= "red", 
                                           ylim= c(-1,3), 
                                           background.track = backgroundColor, 
                                           background.title= backgroundColor, 
                                           fontcolor="black", 
                                           cex.title = 0.72, 
                                           rotation.title = 0, 
                                           lwd.title = 0.5, title.width = 1.5, 
                                           cex.main = 1,
                                           fontcolor.title = "white",
                                           color.scale= "black")


PhastConGVIZTrack<-phastConARXScores%>%DataTrack(type= c("histogram"),
                                                 name= "PhastCon",
                                                 fill="#AB9F77", 
                                                 background.track = "#FFFEDB", 
                                                 background.title= backgroundColor, 
                                                 fontcolor="black", 
                                                 cex.title = 0.72, 
                                                 rotation.title = 90, 
                                                 lwd.title = 0.5, title.width = 1.5, 
                                                 cex.main = 1,
                                                 fontcolor.title = "white",
                                                 color.scale= "black")




phylopGVIZTrack<-phylopARXScores%>%DataTrack(type= c( "histogram"), name= "Phylop",
                                             fill="#B4A3CF",
                                             col.baseline= "red", 
                                             ylim= c(-1,3), 
                                             background.track = backgroundColor, 
                                             background.title= backgroundColor, 
                                             fontcolor="black", 
                                             cex.title = 0.72, 
                                             rotation.title = 90, 
                                             lwd.title = 0.5, title.width = 1.5, 
                                             cex.main = 1,
                                             fontcolor.title = "white",
                                             color.scale= "black")





Itrack <-IdeogramTrack(genome = "hg19",
                       name= "Ideogram",
                       chromosome = "chr7", 
                       background.track = "#FFFEDB", 
                       background.title= backgroundColor, 
                       fontcolor="black", 
                       cex.title = 0.72, 
                       rotation.title = 0, 
                       lwd.title = 0.5, title.width = 1.5, 
                       cex.main = 1,
                       fontcolor.title = "white",
                       color.scale= "black")
Gtrack<-GenomeAxisTrack(range = ARXEnrichedPromoters, 
                        name= "Genome Axis",
                        background.track = "#FFFEDB", 
                        background.title= backgroundColor, 
                        fontcolor="black", 
                        cex.title = 0.72, 
                        rotation.title = 0, 
                        lwd.title = 0.5, title.width = 1.5, 
                        cex.main = 1,
                        fontcolor.title = "white",
                        color.scale= "black")
plotTracks(c(Itrack, Gtrack, ARXEnrichedPromotersTrack ,ARXImage, phylopGVIZTrack, PhastConGVIZTrack),
           from=28473877, to=28473953, chromosome = "chr7", legend= TRUE, 
           main = "Comparison of Conservation methods to ARX MotifsARX Motifs",
           fontcolor="black", 
           cex.title = 0.72,
           lwd.title = 0.5,
           title.width = 1.5, 
           cex.main = 1,
           fontcolor.title = "white",
           color.scale= "black", sizes = c(1.5,1.5,2,2,4,4))
