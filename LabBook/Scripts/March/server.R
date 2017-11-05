
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyServer(function(input, output) {


  

    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)

    # draw the histogram with the specified number of bins
  output$distPlot <- renderPlot({ genes<- genes(txdb, filter=list(gene_id = "170302"), single.strand.genes.only = TRUE)
  sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg19, genes) #getting ARX sequence
  characterSequence <- as.character(sequence)
  DNAstring1 <- DNAStringSet(characterSequence) #ARX dnastring
  windowsize <- Views(DNAstring1, start = 1:12244, width = 10) #chunks of 10
  windowSizeDataFrame <- as.data.frame(windowsize)#convert to dataframe
  baseFrequency <- alphabetFrequency(windowsize)
  windowSizebaseFrequency <- cbind(windowSizeDataFrame, baseFrequency) #combing fre
  GC <- rowSums(windowSizebaseFrequency[c("C","G")]) / rowSums(windowSizebaseFrequency[c("A", "T", "C", "G")])
  
  GCProp <- GRanges(seqnames = rep("chrX", length(windowsize)),
                    ranges = IRanges(
                      start = start(genes) + start(windowsize) - 1,
                      width = width(windowsize)),
                    strand = "-",
                    sequence = as.character(windowsize),
                    GC = GC,
                    seqinfo = seqinfo(genes))
  gcDataTrack<-DataTrack(GCProp)
  plotTracks(gcDataTrack, type = "b", from =25021813 ,to = 25034065)
  })

})
