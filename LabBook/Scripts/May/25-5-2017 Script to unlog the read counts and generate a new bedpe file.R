library(readr)
library(GenomicInteractions)
interInterchromosomes <- read_delim("~/interInterchromosomes.bedpe", 
                                    "\t", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE)%>%as.data.frame
interInterchromosomes$X8<-10^(interInterchromosomes$X8)
exportb("interInterchromosomesInteractions.bedpe")
interInterchromosomes
write.table(x = interInterchromosomes,file = "properReadCounts.bedpe", quote= FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
makeGenomicInteractionsFromFile("/home/a1649239/properReadCounts.bedpe",  type = "bedpe", experiment_name = "Draft HiC Mouse Embyronic", description = "mouseBrain")
