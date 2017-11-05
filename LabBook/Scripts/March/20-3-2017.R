library(MotifDb)
library(seqLogo)
library(Biostrings)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(MotIV)
library(magrittr)


#ARX annotationTRACK
mdb <- MotifDb
matrices.human <- MotifDb::query(mdb, 'hsapiens')
arx <- "170302"
pfm.arx_jolma2013.jaspar <- MotifDb::query(mdb, 'Hsapiens-jolma2013-ARX')[[1]]
seqLogo(pfm.arx_jolma2013.jaspar)
chromosomal.loc <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="gene")[arx]
pcm.arx_jolma2013.jaspar <- round(100 * pfm.arx_jolma2013.jaspar)


arxMotifs <- matchPWM(pcm.arx_jolma2013.jaspar, BSgenome.Hsapiens.UCSC.hg19[["chrX"]], "100%")
ranges<-cbind.data.frame(as.data.frame(ranges(arxMotifs)), as.data.frame(arxMotifs), seqnames = c("chrX"))
ARXMotifsGRangeForChrX<- makeGRangesFromDataFrame(ranges,
                                                  keep.extra.columns=TRUE,
                                                  ignore.strand=FALSE,
                                                  seqinfo=NULL,
                                                  seqnames.field="seqnames",
                                                  start.field="start",
                                                  end.field=c("end", "stop"),
                                                  strand.field="strand",
                                                  starts.in.df.are.0based=FALSE)

arxTrack<- AnnotationTrack(ARXMotifsGRangeForChrX, name= "ARX motif")

##GeneTrack
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
geneTrack <- GeneRegionTrack(txdb, 
                             chromsome = "chrX",
                             name = "Gene Tracks")
plotTracks(geneTrack,
           collapseTranscripts = "longest",
           shape = "arrow", 
           transcriptiAnnontation = "symbol")
#above requires more tweaking its not plotting together or some shit

##DataTrack
genes<- genes(txdb, filter=list(gene_id = "170302"), single.strand.genes.only = TRUE)
sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg19, genes)%>% #getting ARX sequence
characterSequence <- as.character(sequence)
DNAstring1 <- DNAString(characterSequence) #ARX dnastring
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


