#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(gridExtra)
library(Gviz)
library(GenomicInteractions)
library(rtracklayer)
library(magrittr)
library(parallel)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(org.Hs.eg.db)
library(org.Mm.eg.db)



shinyServer(function(input, output) {
 


  ##############################################
  #Render Human GVIZ plot 1
  #############################################33
  output$HumangvizPlot <- renderPlot({

if(!exists("chrM")){
  
  chrM<-input$chrM
  assign("chrM", chrM, .GlobalEnv)
  
  humanIdeogramTrack<-IdeogramTrack(chromosome = input$chrM, genome="hg19",name= "Ideogram")
  gHumanTrack<-GenomeAxisTrack(name= "Axis")
  
  assign("humanIdeogramTrack", humanIdeogramTrack, .GlobalEnv)
  assign("gHumanTrack", gHumanTrack, .GlobalEnv)
  ########################################
  ###HumanInputs
  #######################################
  
  ## Human ChormHMM Tracks
  Pancreas="DataFiles/ChromHMM/human/coMET/E098_15_coreMarks_mnemonics.bed.gz"%>%import()
  PancreasIslets="DataFiles/ChromHMM/human/coMET/E093_15_coreMarks_mnemonics.bed.gz"%>%import()
  fetalBrainFemale="DataFiles/ChromHMM/human/coMET/E082_15_coreMarks_mnemonics.bed.gz"%>%import()
  fetalBrainMale="DataFiles/ChromHMM/human/coMET/E081_15_coreMarks_mnemonics.bed.gz"%>%import()
  H9NeuronCells="DataFiles/ChromHMM/human/coMET/E010_15_coreMarks_mnemonics.bed.gz"%>%import()
  H9NeuronProgenitorCells="DataFiles/ChromHMM/human/coMET/E009_15_coreMarks_mnemonics.bed.gz"%>%import()
  
  
  assign("Pancreas", Pancreas, .GlobalEnv)
  assign("PancreasIslets", PancreasIslets, .GlobalEnv)
  assign("fetalBrainFemale", fetalBrainFemale, .GlobalEnv)
  assign("fetalBrainMale", fetalBrainMale, .GlobalEnv)
  assign("H9NeuronCells", H9NeuronCells, .GlobalEnv)
  assign("H9NeuronProgenitorCells", H9NeuronProgenitorCells, .GlobalEnv)
  
  
  #Enhancers
  EnhancersHuman<-import("DataFiles/Enhancers Track/Human/human_permissive_enhancers_phase_1_and_2.bed.gz")
  
  #All Arx Motif Track!!!! 
  arx6MerTFBSHg19<-readRDS("DataFiles/ChIPseq/Human/ARX6merHg19Sites")
  mcols(arx6MerTFBSHg19)<-cbind.data.frame("Model"= "6mer")
  ArxP6mer4Human<-readRDS("DataFiles/ChIPseq/Human/Plaindromic4SpacedTFBS")
  mcols(ArxP6mer4Human)<-cbind.data.frame("Model"= "P6mer4")
  ArxJolmaHuman<-readRDS("DataFiles/ChIPseq/Human/JolmaTFBS")
  mcols(ArxJolmaHuman)<-cbind.data.frame("Model"= "Jolma")
  ArxT6mer2Human<-readRDS("DataFiles/ChIPseq/Human/ARXTande2SpacedSites")
  mcols(ArxT6mer2Human)<-cbind.data.frame("Model"= "T6mer2")
  
  ## See latter down for impletation of this track
  arxMotifsHumanRaw<-c(arx6MerTFBSHg19,
                       ArxP6mer4Human,
                       ArxJolmaHuman,
                       ArxT6mer2Human)%>%unlist()
  ### Predicted TFBS
  ARXMotifModels<-readRDS("DataFiles/ChIPseq/Human/PredictedARXTFBS")
  
  
  ARXEnhancerMotifs<-ARXMotifModels%>%AnnotationTrack(genome = "hg19", 
                                                      stacking = "dense",
                                                      strand= "*",
                                                           col.line="black", 
                                                      feature= (mcols(ARXMotifModels))$model,
                                                           name="Predicted TFBS")
  #ColoringTrack
  displayPars(ARXEnhancerMotifs) <- list(`6mer` = "#e6194b",
                                         `Tandem2Spaced` = "#3cb44b", 
                                         `Jolma` = "#0082c8", 
                                         `Palindromic 4 Spaced` = "#008080")
  
  
  assign("ARXEnhancerMotifs", ARXEnhancerMotifs, .GlobalEnv)
  assign("EnhancersHuman", EnhancersHuman, .GlobalEnv)
  assign("arxMotifsHumanRaw", arxMotifsHumanRaw, .GlobalEnv)
  
  
  geneTracks<-import("DataFiles/Gene Tracks/Human/hg.bed.gz")
  promoterTracks<-geneTracks%>%promoters()%>%GRanges()
  
  

  
  #Human HiC Data
  interactionsHumanBrain<-readRDS(file= "DataFiles/HiC/Human/SignificantInteractionsBetweenEnhancersContainingARX")%>%InteractionTrack(
    name= "ARX Significant Interactions")
  
  contactProbabilities<- readRDS(file="DataFiles/HiC/Human/contactProbabilitiesHuman")%>%InteractionTrack(
    name= "Contact Probabilities"
  )
  
  
  #Coloring the HiC data
  displayPars(contactProbabilities) = list(col.interactions="red",
                                           col.anchors.line = "black",
                                           interaction.dimension="height", 
                                           interaction.measure ="counts",
                                           plot.trans=FALSE,
                                           plot.outside = TRUE, 
                                           col.outside="lightblue", 
                                           anchor.height = 0.1)
  
  displayPars(interactionsHumanBrain) = list(col.interactions="red",
                                             col.anchors.line = "black",
                                             interaction.dimension="height", 
                                             interaction.measure ="counts",
                                             plot.trans=FALSE,
                                             plot.outside = TRUE, 
                                             col.outside="lightblue", 
                                             anchor.height = 0.1)
  
  assign("promoterTracks", promoterTracks, .GlobalEnv)
  assign("arxMotifsHumanRaw", arxMotifsHumanRaw, .GlobalEnv)
  assign("interactionsHumanBrain", interactionsHumanBrain, .GlobalEnv)
  assign("contactProbabilities", contactProbabilities, .GlobalEnv)

  
  
  
  
  ##ChromHMM Track Generator specifically for humans
  chromHMMTrackGenerator<-function (gen = "hg19",
                                    chr, 
                                    from, 
                                    to,
                                    bedFile, 
                                    featureDisplay = featureDisplay, 
                                    colorcase = "roadmap15") 
  {
    desiredRegion <- subset(get(bedFile), end > from & 
                              start < to & seqnames == chr)
    track <- AnnotationTrack(desiredRegion, 
                             stacking = "dense",
                             col.line="black",
                             feature = (mcols(desiredRegion))$name,
                             genome = "hg19",
                             strand= "*",
                             name = paste(bedFile))
    
    displayPars(track) <- list(`1_TssA` = "#FF0000", `2_TssAFlnk` = "#FF6E00", 
                               `3_TxFlnk` = "#32CD32", `4_Tx` = "#008000", `5_TxWk` = "#006400", 
                               `6_EnhG` = "#C2E105", `7_Enh` = "#FFFF00", `8_ZNF/Rpts` = "#66CDAA", 
                               `9_Het` = "#8A91D0", `10_TssBiv` = "#CD5C5C", `11_BivFlnk` = "#E9967A", 
                               `12_EnhBiv` = "#BDB76B", `13_ReprPC` = "#3A3838", 
                               `14_ReprPCWk` = "#808080", `15_Quies` = "#DCDCDC", 
                               Empty = "#ffffff")
    return(track)
  }
  
  chromHMM_RoadMapAll<-lapply(c("Pancreas",
                                "PancreasIslets",
                                "fetalBrainFemale",
                                "fetalBrainMale",
                                "H9NeuronCells",
                                "H9NeuronProgenitorCells"), function(x){chromHMMTrackGenerator(gen="hg19",
                                                                                                  chr= input$chrM,
                                                                                                  from  = input$fromM,
                                                                                                  to = input$toM,
                                                                                                  bedFile = x,
                                                                                                  featureDisplay = "all", 
                                                                                                  colorcase='roadmap15' )})
  
  # Gene Track with symbols :D
  knownGenes <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                genome="hg19", 
                                chromosome="chrX", 
                                showId=TRUE,
                                geneSymbol=TRUE, 
                                name="UCSC")
  
  symbols <- unlist(mapIds(org.Hs.eg.db, gene(knownGenes),
                           "SYMBOL", "ENTREZID", 
                           multiVals = "first"))
  symbol(knownGenes) <- symbols[gene(knownGenes)]
  
  #Promoter and Motif Track
  promotertrackChromosomeSpecific<-promoterTracks%>%subset(. , 
                                                           seqnames==input$chrM)%>%AnnotationTrack(., name= "PromoterTrack", 
                                                                                                   genome="hg19")
  geneTrackChromosomeSpecific<-knownGenes
  EnhancersHumanChromosomeSpecific<-EnhancersHuman%>%subset(. ,
                                                            seqnames==input$chrM)%>%AnnotationTrack(., name = "Enhancers",
                                                                                                    genome = "hg19")
  
  # Arx All motifs Track
  Arx6merHumanTrack<-subset(arxMotifsHumanRaw, 
                             seqnames==input$chrM & start > input$fromM & end< input$toM)%>%AnnotationTrack(.,genome = "hg19",
                                                                                                            stacking = "dense", 
                                                                                                            col.line="black",
                                                                                                            name="ALL ARX Motifs",
                                                                                                            feature= (mcols(.))$Model)
  
  displayPars(Arx6merHumanTrack) <- list(`6mer` = "#e6194b", 
                                         `T6mer2` = "#3cb44b", 
                                         `Jolma` = "#0082c8", 
                                         `P6mer4` = "#008080")
  

  assign("chromHMM_RoadMapAll", chromHMM_RoadMapAll, .GlobalEnv)
  assign("Arx6merHumanTrack", Arx6merHumanTrack, .GlobalEnv)
  assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
  assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
  assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)
  assign("knownGenes", knownGenes, .GlobalEnv)
  assign("chromHMMTrackGenerator", chromHMMTrackGenerator, .GlobalEnv)
  
  
  
  
  
  if(input$contactProbabilities==TRUE){
    plotTracks(trackList = c(humanIdeogramTrack,
                             gHumanTrack,
                             contactProbabilities, 
                             EnhancersHumanChromosomeSpecific,
                             ARXEnhancerMotifs,
                             Arx6merHumanTrack, 
                             promotertrackChromosomeSpecific, 
                             geneTrackChromosomeSpecific,
                             chromHMM_RoadMapAll), 
               sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
               from =input$fromM, 
               to= input$toM,
               chromosome= input$chrM,
               cex.title = 0.72, 
               rotation.title = 0, 
               showAxis = FALSE, 
               background.title = "white",
               lwd.title = 2, 
               title.width = 2, 
               cex.main = 5, 
               col = NULL, 
               fontcolor.title = "black")
  } else{
    
    plotTracks(trackList =c(humanIdeogramTrack,
                            gHumanTrack,
                            interactionsHumanBrain, 
                            EnhancersHumanChromosomeSpecific,
                            ARXEnhancerMotifs,
                            Arx6merHumanTrack, 
                            promotertrackChromosomeSpecific, 
                            geneTrackChromosomeSpecific,
                            chromHMM_RoadMapAll), 
               sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
               from =input$fromM, 
               to= input$toM,
               chromosome= input$chrM,
               cex.title = 0.72, 
               rotation.title = 0, 
               showAxis = FALSE, 
               background.title = "white",
               lwd.title = 2, 
               title.width = 2, 
               cex.main = 5, 
               col = NULL, 
               fontcolor.title = "black")
  
}}else if(!chrM==input$chrM){ 
  chrM<-input$chrM
  assign("chrM", chrM, .GlobalEnv)
  humanIdeogramTrack<-IdeogramTrack(chromosome = input$chrM, genome="hg19",name= "Ideogram")
  
  
  chromHMM_RoadMapAll<-lapply(c("Pancreas",
                                "PancreasIslets",
                                "fetalBrainFemale",
                                "fetalBrainMale",
                                "H9NeuronCells",
                                "H9NeuronProgenitorCells"), function(x){chromHMMTrackGenerator(gen="hg19",
                                                                                               chr=input$chrM,
                                                                                               from = input$fromM,
                                                                                               to = input$toM,
                                                                                               bedFile = x,
                                                                                               featureDisplay = "all", 
                                                                                               colorcase='roadmap15' )})
 
  # Gene Track with symbols :D
  knownGenes <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                genome="hg19", 
                                chromosome=input$chrM, 
                                showId=TRUE,
                                geneSymbol=TRUE, 
                                name="UCSC")
  symbols <- unlist(mapIds(org.Hs.eg.db, gene(knownGenes),
                           "SYMBOL", "ENTREZID", 
                           multiVals = "first"))
  symbol(knownGenes) <- symbols[gene(knownGenes)]
  
  #Promoter and Motif Track
    promotertrackChromosomeSpecific<-promoterTracks%>%subset(. , 
                                                             seqnames==input$chrM)%>%AnnotationTrack(., name= "PromoterTrack", 
                                                                                                     genome="hg19")
    geneTrackChromosomeSpecific<-knownGenes
    EnhancersHumanChromosomeSpecific<-EnhancersHuman%>%subset(. ,
                                                              seqnames==input$chrM)%>%AnnotationTrack(., name = "Enhancers",
                                                                                                      genome = "hg19")
    
    Arx6merHumanTrack<-subset(arxMotifsHumanRaw, 
                              seqnames==input$chrM & start > input$fromM & end< input$toM)%>%AnnotationTrack(.,genome = "hg19",
                                                                                                             stacking = "dense", 
                                                                                                             col.line="black",
                                                                                                             name="ALL ARX Motifs",
                                                                                                             feature= (mcols(.))$Model)
    displayPars(Arx6merHumanTrack) <- list(`6mer` = "#e6194b", 
                                           `T6mer2` = "#3cb44b", 
                                           `Jolma` = "#0082c8", 
                                           `P6mer4` = "#008080")
    

    
    
    assign("humanIdeogramTrack", humanIdeogramTrack, .GlobalEnv)
    assign("chromHMM_RoadMapAll", chromHMM_RoadMapAll, .GlobalEnv)
    assign("Arx6merHumanTrack", Arx6merHumanTrack, .GlobalEnv)
    assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
    assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
    assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)
    assign("knownGenes", knownGenes, .GlobalEnv)
    
    
    
    
    if(input$contactProbabilities==TRUE){
    plotTracks(trackList = c(humanIdeogramTrack,
                             gHumanTrack,
                             contactProbabilities, 
                             EnhancersHumanChromosomeSpecific,
                             ARXEnhancerMotifs,
                             Arx6merHumanTrack, 
                             promotertrackChromosomeSpecific, 
                             geneTrackChromosomeSpecific,
                             chromHMM_RoadMapAll), 
               sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
               from =input$fromM, 
               to= input$toM,
               chromosome= input$chrM,
               cex.title = 0.72, 
               rotation.title = 0, 
               showAxis = FALSE, 
               background.title = "white",
               lwd.title = 2, 
               title.width = 2, 
               cex.main = 5, 
               col = NULL, 
               fontcolor.title = "black")
      } else{
        
        plotTracks(trackList =c(humanIdeogramTrack,
                                gHumanTrack,
                                interactionsHumanBrain, 
                                EnhancersHumanChromosomeSpecific,
                                ARXEnhancerMotifs,
                                Arx6merHumanTrack, 
                                promotertrackChromosomeSpecific, 
                                geneTrackChromosomeSpecific,
                                chromHMM_RoadMapAll), 
                   sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
                   from =input$fromM, 
                   to= input$toM,
                   chromosome= input$chrM,
                   cex.title = 0.72, 
                   rotation.title = 0, 
                   showAxis = FALSE, 
                   background.title = "white",
                   lwd.title = 2, 
                   title.width = 2, 
                   cex.main = 5, 
                   col = NULL, 
                   fontcolor.title = "black")
      }} else if(input$contactProbabilities==TRUE) {
  
        
        chromHMM_RoadMapAll<-lapply(c("Pancreas",
                                      "PancreasIslets",
                                      "fetalBrainFemale",
                                      "fetalBrainMale",
                                      "H9NeuronCells",
                                      "H9NeuronProgenitorCells"), function(x){chromHMMTrackGenerator(gen="hg19",
                                                                                                        chr=input$chrM,
                                                                                                        from = input$fromM,
                                                                                                        to = input$toM,
                                                                                                        bedFile = x,
                                                                                                        featureDisplay = "all", 
                                                                                                        colorcase='roadmap15' )})
        # Arx All motifs Track
        
        Arx6merHumanTrack<-subset(arxMotifsHumanRaw, 
                                  seqnames==input$chrM & start > input$fromM & end< input$toM)%>%AnnotationTrack(.,genome = "hg19",
                                                                                                                 stacking = "dense", 
                                                                                                                 col.line="black",
                                                                                                                 name="ALL ARX Motifs",
                                                                                                                 feature= (mcols(.))$Model)
        
        displayPars(Arx6merHumanTrack) <- list(`6mer` = "#e6194b", 
                                               `T6mer2` = "#3cb44b", 
                                               `Jolma` = "#0082c8", 
                                               `P6mer4` = "#008080")
        
        
        
        assign("chromHMM_RoadMapAll", chromHMM_RoadMapAll, .GlobalEnv)
        assign("Arx6merHumanTrack", Arx6merHumanTrack, .GlobalEnv)
        assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
        assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
        assign("knownGenes", knownGenes, .GlobalEnv)
        
     
        
  plotTracks(trackList = c(humanIdeogramTrack,
                           gHumanTrack,
                           contactProbabilities, 
                           EnhancersHumanChromosomeSpecific,
                           ARXEnhancerMotifs,
                           Arx6merHumanTrack, 
                           promotertrackChromosomeSpecific, 
                           geneTrackChromosomeSpecific,
                           chromHMM_RoadMapAll), 
             sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
             from =input$fromM, 
             to= input$toM,
             chromosome= input$chrM,
             cex.title = 0.72, 
             rotation.title = 0, 
             showAxis = FALSE, 
             background.title = "white",
             lwd.title = 2, 
             title.width = 2, 
             cex.main = 5, 
             col = NULL, 
             fontcolor.title = "black")
  
}else {
  
  chromHMM_RoadMapAll<-lapply(c("Pancreas",
                                "PancreasIslets",
                                "fetalBrainFemale",
                                "fetalBrainMale",
                                "H9NeuronCells",
                                "H9NeuronProgenitorCells"), function(x){chromHMMTrackGenerator(gen="hg19",
                                                                                               chr=input$chrM,
                                                                                               from = input$fromM,
                                                                                               to = input$toM,
                                                                                               bedFile = x,
                                                                                               featureDisplay = "all", 
                                                                                               colorcase='roadmap15' )})
  # Arx All motifs Track
  
  Arx6merHumanTrack<-subset(arxMotifsHumanRaw, 
                            seqnames==input$chrM & start > input$fromM & end< input$toM)%>%AnnotationTrack(.,genome = "hg19",
                                                                                                           stacking = "dense", 
                                                                                                           col.line="black",
                                                                                                           name="ALL ARX Motifs",
                                                                                                           feature= (mcols(.))$Model)
  
  displayPars(Arx6merHumanTrack) <- list(`6mer` = "#e6194b", 
                                         `T6mer2` = "#3cb44b", 
                                         `Jolma` = "#0082c8", 
                                         `P6mer4` = "#008080")
  
  
  
  plotTracks(trackList = c(humanIdeogramTrack,
                           gHumanTrack,
                           interactionsHumanBrain, 
                           EnhancersHumanChromosomeSpecific,
                           ARXEnhancerMotifs,
                           Arx6merHumanTrack, 
                           promotertrackChromosomeSpecific, 
                           geneTrackChromosomeSpecific,
                           chromHMM_RoadMapAll), 
             sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
             from =input$fromM, 
             to= input$toM,
             chromosome= input$chrM,
             cex.title = 0.72, 
             rotation.title = 0, 
             showAxis = FALSE, 
             background.title = "white",
             lwd.title = 2, 
             title.width = 2, 
             cex.main = 5, 
             col = NULL, 
             fontcolor.title = "black")
  }
    },height = 850,width =  1600)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 
  
  
  ######################################################
  ############MOuse GVIZ PLot 1
  ######################################################

  output$MousegvizPlot <-  renderPlot({
    if(!exists("testesMouse")){
      
      
      mouseIdeogramTrack<-IdeogramTrack(chromosome = input$chrM, genome="mm9",name= "Ideogram")
      gmouseTrack<-GenomeAxisTrack(name= "Axis")
      
      assign("mouseIdeogramTrack", mouseIdeogramTrack, .GlobalEnv)
      assign("gmouseTrack", gmouseTrack, .GlobalEnv)
      ####################################################
      #####MOUSE DATA
      ###########################################################3
      #Mouse ChromHMM inputs
      testesMouse="DataFiles/ChromHMM/mouse/testes_cStates_HMM.bed.gz"%>%import()
      brainMouse= "DataFiles/ChromHMM/mouse/brain_cStates_HMM.bed.gz"%>%import()
      thymusMouse="DataFiles/ChromHMM/mouse/thymus_cStates_HMM.bed.gz"%>%import()
      heartMouse="DataFiles/ChromHMM/mouse/heart_cStates_HMM.bed.gz"%>%import()
      mESCMouse="DataFiles/ChromHMM/mouse/mESC_cStates_HMM.bed.gz"%>%import()
      intestineMouse="DataFiles/ChromHMM/mouse/intestine_cStates_HMM.bed.gz"%>%import()
      
      assign("testesMouse", testesMouse, .GlobalEnv)
      assign("brainMouse", brainMouse, .GlobalEnv)
      assign("thymusMouse", thymusMouse, .GlobalEnv)
      assign("heartMouse", heartMouse, .GlobalEnv)
      assign("mESCMouse", mESCMouse, .GlobalEnv)
      assign("intestineMouse", intestineMouse, .GlobalEnv)
      
      #Mouse Inputs
      Arx6merMouse<-readRDS("DataFiles/ChIPseq/Mouse/ARX6mermm9Sites")
      mcols(Arx6merMouse)<-cbind.data.frame("Model"= "6mer")
      ArxP6mer4Mouse<-readRDS("DataFiles/ChIPseq/Mouse/Plaindromic4Spacedmm9")
      mcols(ArxP6mer4Mouse)<-cbind.data.frame("Model"= "P6mer4")
      ArxJolmaMouse<-readRDS("DataFiles/ChIPseq/Mouse/JolmaTFBS")
      mcols(ArxJolmaMouse)<-cbind.data.frame("Model"= "Jolma")
      ArxT6mer2Mouse<-readRDS("DataFiles/ChIPseq/Mouse/ARXTande2SpacedSites")
      mcols(ArxT6mer2Mouse)<-cbind.data.frame("Model"= "T6mer2")
      
      
      Arx6merMouseRaw<-c(Arx6merMouse, 
                      ArxP6mer4Mouse,
                      ArxT6mer2Mouse,
                      ArxJolmaMouse)%>%unlist()
      
      EnhancersMouse<-import("DataFiles/Enhancers Track/Mouse/mouse_permissive_enhancers_phase_1_and_2.bed.gz")
      geneTracksMouse<-import("DataFiles/Gene Tracks/Mouse/mm9.bed.gz")
      promoterTracksMouse<-promoters(geneTracksMouse)%>%GRanges()
      
      assign("promoterTracksMouse", promoterTracksMouse, .GlobalEnv)
      assign("EnhancersMouse", EnhancersMouse, .GlobalEnv)
      
      
      
      #Hic Data
      contactProbabilitiesMouse<-readRDS("DataFiles/HiC/Mouse/contactProbabilitiesMouse")%>%InteractionTrack(name= "Contact Probabilities")
      interactionBrainMouse<-readRDS("DataFiles/HiC/Mouse/mm9StasticallySignificantInteractions")%>%InteractionTrack(name = "Significant Interactions")
      
      #Coloring the tracks
      displayPars(contactProbabilitiesMouse) = list(col.interactions="red",
                                                    col.anchors.line = "black",
                                                    interaction.dimension="height", 
                                                    interaction.measure ="counts",
                                                    plot.trans=FALSE,
                                                    plot.outside = TRUE, 
                                                    col.outside="lightblue", 
                                                    anchor.height = 0.1)
      displayPars(interactionBrainMouse) = list(col.interactions="red",
                                                col.anchors.line = "black",
                                                interaction.dimension="height", 
                                                interaction.measure ="counts",
                                                plot.trans=FALSE,
                                                plot.outside = TRUE, 
                                                col.outside="lightblue", 
                                                anchor.height = 0.1)
      
      
      
      assign("contactProbabilitiesMouse", contactProbabilitiesMouse, .GlobalEnv)
      assign("interactionBrainMouse", interactionBrainMouse, .GlobalEnv)
      assign("promoterTracksMouse", promoterTracksMouse, .GlobalEnv)
      assign("Arx6merMouseRaw", Arx6merMouseRaw, .GlobalEnv)

      
      #Predicted Arx Motif Track!!!! 
      
      ARXMotifModelsMouse<-readRDS("DataFiles/ChIPseq/Mouse/ArxPredictedTFBS")
      
      ARXEnhancerMotifsMouse<-ARXMotifModelsMouse%>%AnnotationTrack(genome = "mm9", 
                                                                    stacking = "dense", 
                                                                    strand= "*",
                                                            col.line="black", 
                                                            feature= (mcols(ARXMotifModelsMouse))$model,
                                                            name="Predicted Arx TFBS")
      #ColoringTrack
      displayPars(ARXEnhancerMotifsMouse) <- list(`6mer` = "#e6194b", 
                                                  `Tandem2Spaced` = "#3cb44b", 
                                               `Jolma` = "#0082c8", 
                                               `Palindromic 4 Spaced` = "#008080")
      
      

      
      assign("ARXEnhancerMotifsMouse",ARXEnhancerMotifsMouse, .GlobalEnv)
      assign("interactionBrainMouse", interactionBrainMouse, .GlobalEnv)
      
      
      chromHMMTrackGeneratorMouse<-function (gen = "mm9", chr, from, to, bedFile, featureDisplay = featureDisplay, 
                                             colorcase = "roadmap15") 
      {
        desiredRegion <- subset(get(bedFile), end > from & 
                                  start < to & seqnames == chr)
        mcols(desiredRegion)<-cbind.data.frame("name"=(mcols(desiredRegion))$name)
        track<-AnnotationTrack(desiredRegion, stacking = "dense", col.line="black", feature = 
                                 (mcols(desiredRegion))$name, genome = "mm9", strand= "*",
                               name = paste(bedFile))
        if (colorcase == "roadmap15") {
          displayPars(track) <- list(`1_Txn_Elongation` = "#FF0000", `2_Weak_Txn` = "#FF6E00", 
                                     `9_Strong_Enhancer` = "#32CD32", `4_Poised_Enhancer` = "#008000", `5_Active_Promoter` = "#006400", 
                                     `6_Strong_Enhancer` = "#C2E105", `7_Active_Promoter` = "#FFFF00", `8_Strong_Enhancer` = "#66CDAA", 
                                     `9_Txn_Transition` = "#8A91D0", `10_Poised_Promoter` = "#CD5C5C", `11_Repressed` = "#E9967A", 
                                     `15_Insulator` = "#BDB76B", `12_Heterochrom` = "#3A3838", 
                                     `14_Heterochrom` = "#808080", `13_Heterochrom` = "#DCDCDC", 
                                     Empty = "#ffffff")
        } else {
          stop("Invalid in function RoadMap :color choice invalid :\n")
        }
        track
      }
      
      assign("chromHMMTrackGeneratorMouse", chromHMMTrackGeneratorMouse, .GlobalEnv)

      
      #Gene Symbol Track 
      knownGenesMouse <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm9.knownGene,
                                         genome="mm9",
                                         chromosome=input$chrM,
                                         showId=TRUE,
                                         geneSymbol=TRUE, 
                                         name="UCSC")
      symbolsMouse <- unlist(mapIds(org.Mm.eg.db, 
                                    gene(knownGenesMouse),
                                    "SYMBOL", 
                                    "ENTREZID",
                                    multiVals = "first"))
      symbol(knownGenesMouse) <- symbolsMouse[gene(knownGenesMouse)]
      
      
      promotertrackChromosomeSpecificMouse<-promoterTracksMouse%>%subset(. , 
                                                                         seqnames==input$chrM)%>%AnnotationTrack(., name= "PromoterTrack",
                                                                                                                 genome= "mm9", 
                                                                                                                 stacking= "dense")
      geneTrackChromosomeSpecificMouse<-knownGenesMouse
      
      EnhancersMouseChromosomeSpecificMouse<-subset(EnhancersMouse, seqnames==input$chrM)%>%AnnotationTrack(name= "Enhancer Track", 
                                                                                                            stacking= "dense",
                                                                                                            genome= "mm9")
      
      assign("promotertrackChromosomeSpecificMouse", promotertrackChromosomeSpecificMouse , .GlobalEnv)
      assign("geneTrackChromosomeSpecificMouse", geneTrackChromosomeSpecificMouse, .GlobalEnv)
      assign("EnhancersMouseChromosomeSpecificMouse", EnhancersMouseChromosomeSpecificMouse, .GlobalEnv)
      assign("chromHMMTrackGeneratorMouse", chromHMMTrackGeneratorMouse, .GlobalEnv)
      assign("knownGenesMouse", knownGenesMouse, .GlobalEnv)
      
      
      #Base pair and Chormosome specific Tracks
      chromHMM_RoadMapAllMouse<-lapply(c("testesMouse",
                                         "brainMouse",
                                         "thymusMouse",
                                         "heartMouse",
                                         "mESCMouse",
                                         "intestineMouse"), function(x){chromHMMTrackGeneratorMouse(gen="mm9",
                                                                                                       chr=input$chrM,
                                                                                                       from = input$fromM,
                                                                                                       to = input$toM,
                                                                                                       bedFile = x,
                                                                                                       featureDisplay = "all", 
                                                                                                       colorcase='roadmap15' )})
      
      
      
      
      ### Raw Motif Traack

      Arx6merMouseTrack<-subset(Arx6merMouseRaw, seqnames==input$chrM & start > input$fromM & end< input$toM)%>%
        AnnotationTrack(.,genome = "mm9", 
                        stacking = "dense",
                        strand= "*",
                        name= "All Motifs",
                        feature= (mcols(.))$Model)
      
      
      
      displayPars(Arx6merMouseTrack) <- list(`6mer` = "#e6194b", 
                                             `T6mer2` = "#3cb44b", 
                                             `Jolma` = "#0082c8", 
                                             `P6mer4` = "#008080")
      
      
      
      assign( "chromHMM_RoadMapAllMouse",chromHMM_RoadMapAllMouse, .GlobalEnv)
      assign( "Arx6merMouseTrack",Arx6merMouseTrack, .GlobalEnv)
      
      
      if(input$contactProbabilities==TRUE){
        
        plotTracks(trackList = c(mouseIdeogramTrack,
                                 gmouseTrack,
                                 contactProbabilitiesMouse, 
                                 EnhancersMouseChromosomeSpecificMouse,
                                 ARXEnhancerMotifsMouse,
                                 Arx6merMouseTrack,
                                 promotertrackChromosomeSpecificMouse, 
                                 geneTrackChromosomeSpecificMouse,
                                 chromHMM_RoadMapAllMouse), 
                   sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
                   from =input$fromM, 
                   to= input$toM,
                   chromosome= input$chrM,
                   cex.title = 0.72, 
                   rotation.title = 0, 
                   showAxis = FALSE, 
                   background.title = "white",
                   lwd.title = 2, 
                   title.width = 2, 
                   cex.main = 5, 
                   col = NULL, 
                   fontcolor.title = "black")
      } else{
        plotTracks(trackList =  c(mouseIdeogramTrack,
                                  gmouseTrack,
                                  interactionBrainMouse, 
                                  EnhancersMouseChromosomeSpecificMouse,
                                  ARXEnhancerMotifsMouse,
                                  Arx6merMouseTrack,
                                  promotertrackChromosomeSpecificMouse, 
                                  geneTrackChromosomeSpecificMouse,
                                  chromHMM_RoadMapAllMouse), 
                   sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
                   from =input$fromM, 
                   to= input$toM,
                   chromosome= input$chrM,
                   cex.title = 0.72, 
                   rotation.title = 0, 
                   showAxis = FALSE, 
                   background.title = "white",
                   lwd.title = 2, 
                   title.width = 2, 
                   cex.main = 5, 
                   col = NULL, 
                   fontcolor.title = "black")
        
        
      
   }
      }else if(!chrM==input$chrM ){ 
      chrM<- input$chrM
      assign("chrM", chrM, .GlobalEnv)
      
      #Gene Symbol Track 
      knownGenesMouse <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm9.knownGene,
                                         genome="mm9",
                                    chromosome=input$chrM,
                                    showId=TRUE,
                                    geneSymbol=TRUE, 
                                    name="UCSC")
      symbolsMouse <- unlist(mapIds(org.Mm.eg.db, 
                                    gene(knownGenesMouse),
                               "SYMBOL", 
                               "ENTREZID",
                               multiVals = "first"))
      symbol(knownGenesMouse) <- symbolsMouse[gene(knownGenesMouse)]
      
      
      promotertrackChromosomeSpecificMouse<-promoterTracksMouse%>%subset(. , 
                                                                         seqnames==input$chrM)%>%AnnotationTrack(., name= "PromoterTrack",
                                                                                                                 genome= "mm9", 
                                                                                                                 stacking= "dense")
      geneTrackChromosomeSpecificMouse<-knownGenesMouse
      
      EnhancersMouseChromosomeSpecificMouse<-subset(EnhancersMouse, seqnames==input$chrM)%>%AnnotationTrack(name= "Enhancer Track", 
                                                                                                            stacking= "dense",
                                                                                                            genome= "mm9")
      
      
      #Base pair and Chormosome specific Tracks
      
      #Base pair and Chormosome specific Tracks
      chromHMM_RoadMapAllMouse<-lapply(c("testesMouse",
                                         "brainMouse",
                                         "thymusMouse",
                                         "heartMouse",
                                         "mESCMouse",
                                         "intestineMouse"), function(x){chromHMMTrackGeneratorMouse(gen="mm9",
                                                                                                    chr=input$chrM,
                                                                                                    from = input$fromM,
                                                                                                    to = input$toM,
                                                                                                    bedFile = x,
                                                                                                    featureDisplay = "all", 
                                                                                                    colorcase='roadmap15' )})
      
      Arx6merMouseTrack<-subset(Arx6merMouseRaw, seqnames==input$chrM & start > input$fromM & end< input$toM)%>%
        AnnotationTrack(.,genome = "mm9", 
                        stacking = "dense",
                        strand= "*",
                        name= "All Motifs",
                        feature= (mcols(.))$Model)
      
      
      
      displayPars(Arx6merMouseTrack) <- list(`6mer` = "#e6194b", 
                                             `T6mer2` = "#3cb44b", 
                                             `Jolma` = "#0082c8",
                                             `P6mer4` = "#008080")
      
      
      
      
      
      assign( "chromHMM_RoadMapAllMouse",chromHMM_RoadMapAllMouse, .GlobalEnv)
      assign( "Arx6merMouseTrack",Arx6merMouseTrack, .GlobalEnv)
      assign("promotertrackChromosomeSpecificMouse", promotertrackChromosomeSpecificMouse , .GlobalEnv)
      assign("geneTrackChromosomeSpecificMouse", geneTrackChromosomeSpecificMouse, .GlobalEnv)
      assign("EnhancersMouseChromosomeSpecificMouse", EnhancersMouseChromosomeSpecificMouse, .GlobalEnv)
      assign("knownGenesMouse", knownGenesMouse, .GlobalEnv)
      
      if(input$contactProbabilities==TRUE){
        
        plotTracks(trackList = c(mouseIdeogramTrack,
                                 gmouseTrack,
                                 contactProbabilitiesMouse, 
                                 EnhancersMouseChromosomeSpecificMouse,
                                 ARXEnhancerMotifsMouse,
                                 Arx6merMouseTrack,
                                 promotertrackChromosomeSpecificMouse, 
                                 geneTrackChromosomeSpecificMouse,
                                 chromHMM_RoadMapAllMouse), 
                   sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
                   from =input$fromM, 
                   to= input$toM,
                   chromosome= input$chrM,
                   cex.title = 0.72, 
                   rotation.title = 0, 
                   showAxis = FALSE, 
                   background.title = "white",
                   lwd.title = 2, 
                   title.width = 2, 
                   cex.main = 5, 
                   col = NULL, 
                   fontcolor.title = "black")
      } else{
        plotTracks(trackList = c(mouseIdeogramTrack,
                                 gmouseTrack,
                                 contactProbabilitiesMouse, 
                                 EnhancersMouseChromosomeSpecificMouse,
                                 ARXEnhancerMotifsMouse,
                                 Arx6merMouseTrack,
                                 promotertrackChromosomeSpecificMouse, 
                                 geneTrackChromosomeSpecificMouse,
                                 chromHMM_RoadMapAllMouse), 
                   sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
                   from =input$fromM, 
                   to= input$toM,
                   chromosome= input$chrM,
                   cex.title = 0.72, 
                   rotation.title = 0, 
                   showAxis = FALSE, 
                   background.title = "white",
                   lwd.title = 2, 
                   title.width = 2, 
                   cex.main = 5, 
                   col = NULL, 
                   fontcolor.title = "black")
        
        
        
        
        
        
      }} else if(input$contactProbabilities==TRUE) {
        #If there is just a base pair change
        
        
        
        
        #Base pair and Chormosome specific Tracks
        chromHMM_RoadMapAllMouse<-lapply(c("testesMouse",
                                           "brainMouse",
                                           "thymusMouse",
                                           "heartMouse",
                                           "mESCMouse",
                                           "intestineMouse"), function(x){chromHMMTrackGeneratorMouse(gen="mm9",
                                                                                                      chr=input$chrM,
                                                                                                      from = input$fromM,
                                                                                                      to = input$toM,
                                                                                                      bedFile = x,
                                                                                                      featureDisplay = "all", 
                                                                                                      colorcase='roadmap15' )})
        
        Arx6merMouseTrack<-subset(Arx6merMouseRaw, seqnames==input$chrM & start > input$fromM & end< input$toM)%>%
          AnnotationTrack(.,genome = "mm9", 
                          stacking = "dense",
                          strand= "*",
                          name= "All Motifs",
                          feature= (mcols(.))$Model)



displayPars(Arx6merMouseTrack) <- list(`6mer` = "#e6194b", 
                                       `T6mer2` = "#3cb44b", 
                                       `Jolma` = "#0082c8", 
                                       `P6mer4` = "#008080")


        assign( "chromHMM_RoadMapAllMouse",chromHMM_RoadMapAllMouse, .GlobalEnv)
        assign( "Arx6merMouseTrack",Arx6merMouseTrack, .GlobalEnv)
        assign("promotertrackChromosomeSpecificMouse", promotertrackChromosomeSpecificMouse , .GlobalEnv)
        assign("geneTrackChromosomeSpecificMouse", geneTrackChromosomeSpecificMouse, .GlobalEnv)
        assign("EnhancersMouseChromosomeSpecificMouse", EnhancersMouseChromosomeSpecificMouse, .GlobalEnv)
        assign("knownGenesMouse", knownGenesMouse, .GlobalEnv)
        

        plotTracks(trackList = c(mouseIdeogramTrack,
                                 gmouseTrack,
                                 contactProbabilitiesMouse, 
                                 EnhancersMouseChromosomeSpecificMouse,
                                 ARXEnhancerMotifsMouse,
                                 Arx6merMouseTrack,
                                 promotertrackChromosomeSpecificMouse, 
                                 geneTrackChromosomeSpecificMouse,
                                 chromHMM_RoadMapAllMouse), 
                   sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
                   from =input$fromM, 
                   to= input$toM,
                   chromosome= input$chrM,
                   cex.title = 0.72, 
                   rotation.title = 0, 
                   showAxis = FALSE, 
                   background.title = "white",
                   lwd.title = 2, 
                   title.width = 2, 
                   cex.main = 5, 
                   col = NULL, 
                   fontcolor.title = "black")
        
      }else {
        
        #Base pair and Chormosome specific Tracks
        chromHMM_RoadMapAllMouse<-lapply(c("testesMouse",
                                           "brainMouse",
                                           "thymusMouse",
                                           "heartMouse",
                                           "mESCMouse",
                                           "intestineMouse"), function(x){chromHMMTrackGeneratorMouse(gen="mm9",
                                                                                                      chr=input$chrM,
                                                                                                      from = input$fromM,
                                                                                                      to = input$toM,
                                                                                                      bedFile = x,
                                                                                                      featureDisplay = "all", 
                                                                                                      colorcase='roadmap15' )})
        
        
        Arx6merMouseTrack<-subset(Arx6merMouseRaw, seqnames==input$chrM & start > input$fromM & end< input$toM)%>%
          AnnotationTrack(.,genome = "mm9", 
                          stacking = "dense",
                          strand= "*",
                          name= "All Motifs",
                          feature= (mcols(.))$Model)



displayPars(Arx6merMouseTrack) <- list(`6mer` = "#e6194b",
                                       `T6mer2` = "#3cb44b", 
                                       `Jolma` = "#0082c8", 
                                       `P6mer4` = "#008080")
        
        assign( "chromHMM_RoadMapAllMouse",chromHMM_RoadMapAllMouse, .GlobalEnv)
        assign( "Arx6merMouseTrack",Arx6merMouseTrack, .GlobalEnv)
        
        
        
        
        
        plotTracks(trackList =c(mouseIdeogramTrack,
                                gmouseTrack,
                                interactionBrainMouse, 
                                EnhancersMouseChromosomeSpecificMouse,
                                ARXEnhancerMotifsMouse,
                                Arx6merMouseTrack,
                                promotertrackChromosomeSpecificMouse, 
                                geneTrackChromosomeSpecificMouse,
                                chromHMM_RoadMapAllMouse), 
                   sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
                   from =input$fromM, 
                   to= input$toM,
                   chromosome= input$chrM,
                   cex.title = 0.72, 
                   rotation.title = 0, 
                   showAxis = FALSE, 
                   background.title = "white",
                   lwd.title = 2, 
                   title.width = 2, 
                   cex.main = 5, 
                   col = NULL, 
                   fontcolor.title = "black")
      }
  },height = 850,width =  1600)
  
  
  output$LegendsPlot<- renderImage({
    
    list(
      src = "www/EpigenomicsRoadMapLegendHMM.jpeg",
      contentType = "image/jpeg",
      alt = "Human/Epigenomics Road Map Legend"
    )
    
  }, deleteFile = FALSE)
  
  
  output$LegendsPlotMouse<- renderImage({
    
    list(
      src = "www/mm9ChromHMMStates.jpeg",
      contentType = "image/jpeg",
      alt = "Human/Epigenomics Road Map Legend"
    )
    
  }, deleteFile = FALSE)
  
  
})



