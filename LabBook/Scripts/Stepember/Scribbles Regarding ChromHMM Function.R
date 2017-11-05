chr<-"chr1"
start<-1
end<-100000000
bedFilePath<-"PancreasHMM"
featureDisplay<-"all"
gen<-"hg19"
chromHMMTrackGenerator<-function (gen = "hg19", chr, start, end, bedFilePath, featureDisplay = featureDisplay, 
                                  colorcase = "roadmap15") 
{
  if (is.null(gen)) {
    stop("Invalid function chromHMM_RoadMap :gen null:\n")
  }
  if (is.null(chr)) {
    stop("Invalid function chromHMM_RoadMap :chr null:\n")
  }
  if (is.null(start)) {
    stop("Invalid function chromHMM_RoadMap :start null:\n")
  }
  if (is.null(end)) {
    stop("Invalid function chromHMM_RoadMap :end null:\n")
  }
  if (is.null(bedFilePath)) {
    stop("Invalid function chromHMM_RoadMap :bedFilePath null:\n")
  }
  chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  bedFile <-import(get(bedFilePath))
  desiredRegion <- subset(bedFile, end > from & 
                            start < to & chromosome_name == chr)
  desiredRegionDisplay <- desiredRegion
  if (!("all" %in% featureDisplay)) {
    desiredRegionDisplay <- desiredRegion[which(desiredRegion$feature_type_name %in% 
                                                  featureDisplay), ]
  }
  if (nrow(desiredRegionDisplay) == 0) {
    desiredRegionDisplay <- data.frame(nrow = 1, ncol = 4)
    desiredRegionDisplay[1, 1] <- chr
    desiredRegionDisplay[1, 2] <- start
    desiredRegionDisplay[1, 3] <- end
    desiredRegionDisplay[1, 4] <- "Empty"
  }
  track <- AnnotationTrack(genome = gen, chromosome = chr, 
                           strand = "*", start = desiredRegionDisplay[, 2], 
                           end = desiredRegionDisplay[, 3], 
                           feature = desiredRegionDisplay[, 4], 
                           name = paste(bedFilePath, " chromHMM RoadMap"), stacking = "dense", col.line = "black", 
                           col = NULL, collapse = FALSE)
  if (colorcase == "roadmap15") {
    displayPars(track) <- list(`1_TssA` = "#FF0000", `2_TssAFlnk` = "#FF6E00", 
                               `3_TxFlnk` = "#32CD32", `4_Tx` = "#008000", `5_TxWk` = "#006400", 
                               `6_EnhG` = "#C2E105", `7_Enh` = "#FFFF00", `8_ZNF/Rpts` = "#66CDAA", 
                               `9_Het` = "#8A91D0", `10_TssBiv` = "#CD5C5C", `11_BivFlnk` = "#E9967A", 
                               `12_EnhBiv` = "#BDB76B", `13_ReprPC` = "#3A3838", 
                               `14_ReprPCWk` = "#808080", `15_Quies` = "#DCDCDC", 
                               Empty = "#ffffff")
  }
  else if (colorcase == "roadmap18") {
    displayPars(data_trackfunc) <- list(`1_TssA` = "#FF0000", 
                                        `2_TssFlnk` = "#FF4500", `3_TssFlnkU` = "#FF4500", 
                                        `4_TssFlnkD` = "#FF4500", `5_Tx` = "#008000", `6_TxWk` = "#006400", 
                                        `7_EnhG1` = "#C2FF05", `8_EnhG2` = "#C2FF05", `9_EnhA1` = "#FFC34D", 
                                        `10_EnhA2` = "#FFC34D", `11_EnhWk` = "#FFFF00", `12_ZNF/Rpts` = "#66CDAA", 
                                        `13_Het` = "#8A91D0", `14_TssBiv` = "#CD5C5C", `15_EnhBiv` = "#BDB76B", 
                                        `16_ReprPC` = "#808080", `17_ReprPC` = "#C0C0C0", 
                                        `18_Quies` = "#FFFFFF")
  }
  else if (colorcase == "comet18") {
    displayPars(data_trackfunc) <- list(`1_TssA` = "#FF0000", 
                                        `2_TssFlnk` = "#FF6E00", `3_TssFlnkU` = "#FF9300", 
                                        `4_TssFlnkD` = "#DA7B08", `5_Tx` = "#008000", `6_TxWk` = "#006400", 
                                        `7_EnhG1` = "#C2FF05", `8_EnhG2` = "#C2FFBD", `9_EnhA1` = "#FE00DB", 
                                        `10_EnhA2` = "#FFA7D6", `11_EnhWk` = "#FFFF00", `12_ZNF/Rpts` = "#66CDAA", 
                                        `13_Het` = "#8A91D0", `14_TssBiv` = "#CD5C5C", `15_EnhBiv` = "#BDB76B", 
                                        `16_ReprPC` = "#323232", `17_ReprPC` = "#AFAFAF", 
                                        `18_Quies` = "#DCDCDC")
  }
  else if (colorcase == "roadmap25") {
    displayPars(data_trackfunc) <- list(`1_TssA` = "#FF0000", 
                                        `2_PromU` = "#FF4500", `3_PromD1` = "#FF4500", `4_PromD2` = "#FF4500", 
                                        `5_Tx5???` = "#008000", `6_Tx` = "#008000", `7_Tx3???` = "#008000", 
                                        `8_TxWk` = "#009600", `9_TxReg` = "#C2FF05", `10_TxEnh5???` = "#C2FF05", 
                                        `11_TxEnh3???` = "#C2FF05", `12_TxEnhW` = "#C2FF05", 
                                        `13_EnhA1` = "#FFC34D", `14_EnhA2` = "#FFC34D", `15_EnhAF` = "#FFC34D", 
                                        `16_EnhW1` = "#FFFF00", `17_EnhW2` = "#FFFF00", `18_EnhAc` = "#FFFF00", 
                                        `19_DNase` = "#FFFF66", `20_ZNF/Rpts` = "#66CDAA", 
                                        `21_Het` = "#8A91D0", `22_PromP` = "#E6B8B7", `23_PromBiv` = "#7030A0", 
                                        `24_ReprPC` = "#808080", `25_Quies` = "#FFFFFF")
  }
  else if (colorcase == "comet25") {
    displayPars(data_trackfunc) <- list(`1_TssA` = "#FF0000", 
                                        `2_PromU` = "#FC6D00", `3_PromD1` = "#DD8100", `4_PromD2` = "#AD7622", 
                                        `5_Tx5???` = "#008000", `6_Tx` = "#004D00", `7_Tx3???` = "#009462", 
                                        `8_TxWk` = "#00FE00", `9_TxReg` = "#00FFFF", `10_TxEnh5???` = "#009FFF", 
                                        `11_TxEnh3???` = "#0028FF", `12_TxEnhW` = "#0000AE", 
                                        `13_EnhA1` = "#FF00FF", `14_EnhA2` = "#FFB2FF", `15_EnhAF` = "#FFD8FF", 
                                        `16_EnhW1` = "#FFFF00", `17_EnhW2` = "#E3FF8C", `18_EnhAc` = "#FFD500", 
                                        `19_DNase` = "#FFFFC2", `20_ZNF/Rpts` = "#66CDAA", 
                                        `21_Het` = "#8A91D0", `22_PromP` = "#E6B8B7", `23_PromBiv` = "#7030A0", 
                                        `24_ReprPC` = "#646464", `25_Quies` = "#DCDCDC")
  }
  else {
    stop("Invalid in function RoadMap :color choice invalid :\n")
  }
  track
}



chromHMMTrackGeneratorMouse<-function (gen = "mm9", chr, from, to, bedFilePath, featureDisplay = featureDisplay, 
                                       colorcase = "roadmap15") 
{
  bedFile <-import(get(bedFilePath))
  desiredRegion <- subset(bedFile, end > from & 
                            start < to & seqnames == chr)
  
  if (nrow(desiredRegionDisplay) == 0) {
    desiredRegionDisplay <- data.frame(nrow = 1, ncol = 4)
    desiredRegionDisplay[1, 1] <- chr
    desiredRegionDisplay[1, 2] <- start
    desiredRegionDisplay[1, 3] <- end
    desiredRegionDisplay[1, 4] <- "Empty"
  }
  track<-AnnotationTrack(desiredRegion, stacking = "dense", col.line="black", feature = 
                           (mcols(desiredRegion))$name, genome = "mm9", strand= "*",
                         name = paste(bedFilePath, " chromHMM RoadMap"))
  if (colorcase == "roadmap15") {
    displayPars(track) <- list(`1_Txn_Elongation` = "#FF0000", `2_Weak_Txn` = "#FF6E00", 
                               `9_Strong_Enhancer` = "#32CD32", `4_Poised_Enhancer` = "#008000", `5_Active_Promoter` = "#006400", 
                               `6_Strong_Enhancer` = "#C2E105", `7_Active_Promoter` = "#FFFF00", `8_Strong_Enhancer` = "#66CDAA", 
                               `9_Txn_Transition` = "#8A91D0", `10_Poised_Promoter` = "#CD5C5C", `11_Repressed` = "#E9967A", 
                               `12_Heterochrom` = "#BDB76B", `13_Heterochrom` = "#3A3838", 
                               `14_Heterochrom` = "#808080", `15_Insulator` = "#DCDCDC", 
                               Empty = "#ffffff")
  }
  else if (colorcase == "roadmap18") {
    displayPars(data_trackfunc) <- list(`1_TssA` = "#FF0000", 
                                        `2_TssFlnk` = "#FF4500", `3_TssFlnkU` = "#FF4500", 
                                        `4_TssFlnkD` = "#FF4500", `5_Tx` = "#008000", `6_TxWk` = "#006400", 
                                        `7_EnhG1` = "#C2FF05", `8_EnhG2` = "#C2FF05", `9_EnhA1` = "#FFC34D", 
                                        `10_EnhA2` = "#FFC34D", `11_EnhWk` = "#FFFF00", `12_ZNF/Rpts` = "#66CDAA", 
                                        `13_Het` = "#8A91D0", `14_TssBiv` = "#CD5C5C", `15_EnhBiv` = "#BDB76B", 
                                        `16_ReprPC` = "#808080", `17_ReprPC` = "#C0C0C0", 
                                        `18_Quies` = "#FFFFFF")
  }
  else if (colorcase == "comet18") {
    displayPars(data_trackfunc) <- list(`1_TssA` = "#FF0000", 
                                        `2_TssFlnk` = "#FF6E00", `3_TssFlnkU` = "#FF9300", 
                                        `4_TssFlnkD` = "#DA7B08", `5_Tx` = "#008000", `6_TxWk` = "#006400", 
                                        `7_EnhG1` = "#C2FF05", `8_EnhG2` = "#C2FFBD", `9_EnhA1` = "#FE00DB", 
                                        `10_EnhA2` = "#FFA7D6", `11_EnhWk` = "#FFFF00", `12_ZNF/Rpts` = "#66CDAA", 
                                        `13_Het` = "#8A91D0", `14_TssBiv` = "#CD5C5C", `15_EnhBiv` = "#BDB76B", 
                                        `16_ReprPC` = "#323232", `17_ReprPC` = "#AFAFAF", 
                                        `18_Quies` = "#DCDCDC")
  }
  else if (colorcase == "roadmap25") {
    displayPars(data_trackfunc) <- list(`1_TssA` = "#FF0000", 
                                        `2_PromU` = "#FF4500", `3_PromD1` = "#FF4500", `4_PromD2` = "#FF4500", 
                                        `5_Tx5???` = "#008000", `6_Tx` = "#008000", `7_Tx3???` = "#008000", 
                                        `8_TxWk` = "#009600", `9_TxReg` = "#C2FF05", `10_TxEnh5???` = "#C2FF05", 
                                        `11_TxEnh3???` = "#C2FF05", `12_TxEnhW` = "#C2FF05", 
                                        `13_EnhA1` = "#FFC34D", `14_EnhA2` = "#FFC34D", `15_EnhAF` = "#FFC34D", 
                                        `16_EnhW1` = "#FFFF00", `17_EnhW2` = "#FFFF00", `18_EnhAc` = "#FFFF00", 
                                        `19_DNase` = "#FFFF66", `20_ZNF/Rpts` = "#66CDAA", 
                                        `21_Het` = "#8A91D0", `22_PromP` = "#E6B8B7", `23_PromBiv` = "#7030A0", 
                                        `24_ReprPC` = "#808080", `25_Quies` = "#FFFFFF")
  }
  else if (colorcase == "comet25") {
    displayPars(data_trackfunc) <- list(`1_TssA` = "#FF0000", 
                                        `2_PromU` = "#FC6D00", `3_PromD1` = "#DD8100", `4_PromD2` = "#AD7622", 
                                        `5_Tx5???` = "#008000", `6_Tx` = "#004D00", `7_Tx3???` = "#009462", 
                                        `8_TxWk` = "#00FE00", `9_TxReg` = "#00FFFF", `10_TxEnh5???` = "#009FFF", 
                                        `11_TxEnh3???` = "#0028FF", `12_TxEnhW` = "#0000AE", 
                                        `13_EnhA1` = "#FF00FF", `14_EnhA2` = "#FFB2FF", `15_EnhAF` = "#FFD8FF", 
                                        `16_EnhW1` = "#FFFF00", `17_EnhW2` = "#E3FF8C", `18_EnhAc` = "#FFD500", 
                                        `19_DNase` = "#FFFFC2", `20_ZNF/Rpts` = "#66CDAA", 
                                        `21_Het` = "#8A91D0", `22_PromP` = "#E6B8B7", `23_PromBiv` = "#7030A0", 
                                        `24_ReprPC` = "#646464", `25_Quies` = "#DCDCDC")
  }
  else {
    stop("Invalid in function RoadMap :color choice invalid :\n")
  }
  track
}

chromHMM_RoadMapAll<-lapply(c("PancreasHMM",
                              "PancreasIsletsHMM",
                              "fetalBrainFemaleHMM",
                              "fetalBrainMaleHMM",
                              "H9NeuronCellsHMM",
                              "H9NeuronProgenitorCellsHMM"), function(x){chromHMMTrackGenerator(gen="hg19",
                                                                                                "chr1",
                                                                                                start = 1000,
                                                                                                end = 248956422,
                                                                                                bedFilePath = x,
                                                                                                featureDisplay = "all", 
                                                                                                colorcase='roadmap15' )})


(mcols(desiredRegion))$name
