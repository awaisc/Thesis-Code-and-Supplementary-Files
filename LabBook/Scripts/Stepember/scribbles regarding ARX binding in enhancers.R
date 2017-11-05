
ARXMouse6Mer<-readRDS("~/DataFiles/ChIPseq/Mouse/ARX6mermm9Sites")
mcols(ARXMouse6Mer)<-cbind.data.frame("Model"="6Mer")

ARXTandem2Mouse<-readRDS("~/DataFiles/ChIPseq/Mouse/ARXTande2SpacedSites")
mcols(ARXTandem2Mouse)<-cbind.data.frame("Model"="ARXTandem2")

ARXJolmaMouse<-readRDS("~/DataFiles/ChIPseq/Mouse/Jolmamm9Sites")
mcols(ARXJolmaMouse)<-cbind.data.frame("Model"="Jolma")

Plaindromic4SpacedTFBSMouse<-readRDS("~/DataFiles/ChIPseq/Mouse/Plaindromic4Spacedmm9")
mcols(Plaindromic4SpacedTFBSMouse)<-cbind.data.frame("Model"="ARXPlaindromic4Spaced")


ARXMotifModelsMouse<-c(ARXMouse6Mer, ARXTandem2Mouse, ARXJolmaMouse, Plaindromic4SpacedTFBSMouse)%>%unlist()

MouseEnhancerMotifs<-subsetByOverlaps(ARXMotifModelsMouse,EnhancersMouse)

ARXEnhancerMotifsMouse<-MouseEnhancerMotifs%>%AnnotationTrack(genome = "mm9", stacking = "dense", strand= "*",
                                                         col.line="black", feature= (mcols(MouseEnhancerMotifs))$Model)

displayPars(ARXEnhancerMotifsMouse) <- list(`6Mer` = "#FF0000", `ARXTandem2` = "#FF6E00", 
                           `Jolma` = "#32CD32", `ARXPlaindromic4Spaced` = "#99CD32")




plotTracks(trackList = c(contactProbabilitiesMouse,
                         EnhancersMouseChromosomeSpecificMouse,
                         ARXEnhancerMotifsMouse,
                         Arx6merMouseTrack, 
                         EnhancersMouseChromosomeSpecificMouse, 
                         promotertrackChromosomeSpecificMouse, 
                         knownGenesMouse,
                         chromHMM_RoadMapAllMouse), 
           from =90000000, 
           to= 91000000,
           chromosome= "chrX",
           cex.title = 0.72, 
           rotation.title = 0, 
           showAxis = FALSE, 
           background.title = "white",
           lwd.title = 2, 
           title.width = 2, 
           cex.main = 5, 
           col = NULL, 
           fontcolor.title = "black")

EnhancersHuman<-import("~/DataFiles/Enhancer Tracks/Human/human_permissive_enhancers_phase_1_and_2.bed")

#Motifs In Enhancers
ARXHuman6Mer<-readRDS("~/DataFiles/ChIPseq/Human/ARX6merHg19Sites")
mcols(ARXHuman6Mer)<-cbind.data.frame("Model"="6Mer")
ARXTandem2<-readRDS("~/DataFiles/ChIPseq/Human/ARXTande2SpacedSites")
mcols(ARXTandem2)<-cbind.data.frame("Model"="ARXTandem2")
ARXHumanJolma<-readRDS("~/DataFiles/ChIPseq/Human/JolmaTFBS")
mcols(ARXHumanJolma)<-cbind.data.frame("Model"="Jolma")
Plaindromic4SpacedTFBS<-readRDS("~/DataFiles/ChIPseq/Human/Plaindromic4SpacedTFBS")
mcols(Plaindromic4SpacedTFBS)<-cbind.data.frame("Model"="ARXPlaindromic4Spaced")

ARXMotifModels<-c(ARXHuman6Mer, ARXTandem2, ARXHumanJolma, Plaindromic4SpacedTFBS)%>%unlist()

HumanEnhancerMotifs<-subsetByOverlaps(ARXMotifModels,EnhancersHuman)

ARXEnhancerMotifs<-HumanEnhancerMotifs%>%AnnotationTrack(genome = "hg19", stacking = "dense", strand= "*",
                                                         col.line="black", feature= (mcols(HumanEnhancerMotifs))$Model,
                                                         name="ARX Motifs In Enhancers")
#ColoringTrack
displayPars(ARXEnhancerMotifs) <- list(`6Mer` = "#FF0000", `ARXTandem2` = "#FF6E00", 
                                       `Jolma` = "#32CD32", `ARXPlaindromic4Spaced` = "#99CD32")



plotTracks(trackList = c(contactProbabilitiesMouse,
                         EnhancersMouseChromosomeSpecificMouse,
                         ARXEnhancerMotifsMouse,
                         Arx6merMouseTrack, 
                         EnhancersMouseChromosomeSpecificMouse, 
                         promotertrackChromosomeSpecificMouse, 
                         knownGenesMouse,
                         chromHMM_RoadMapAllMouse), 
           from =90000000, 
           to= 91000000,
           chromosome= "chrX",
           cex.title = 0.72, 
           rotation.title = 0, 
           showAxis = FALSE, 
           background.title = "white",
           lwd.title = 2, 
           title.width = 2, 
           cex.main = 5, 
           col = NULL, 
           fontcolor.title = "black")

