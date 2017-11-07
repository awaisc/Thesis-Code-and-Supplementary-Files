In this primary fold you will find a series of folders.

Gene Ontology and Phenotype Enrichement Folder: 
Contains the full list of GO terms for gene ontology conducted on gene lists that were predicted to be regulated by Arx TFBS

HumanShinyAPP Folder:

This contains all the dataFiles required to load up the shiny App. 

Please download it and open the server/ui script in R.

PLease set the working directory to the location of the folder with the function
setwd()
Then the app should run, if there are any crashes upon the first load up of the app, please clear the environment and reload the app.


LabBook:
Contains knitted word docs for our R and some bash scripts over the year.

N2a Cells ChIP+Diff
Contains the raw homer output files for our promoter analysis on geens that were differenitally expressed in N2a cells and contianed a ChIP enriched region in the promoter


Non4merContainingPeaks:
The full denovo analysis results from Quille et al's (2011) Arx ChIP-chip peaks not containg a 6mer motif, 4mer or Jolma's model peaks

SubPal ChIP + Diff:
The raw output from Homers promoter analysis of promoters of genes with a Chip-Enriched region in the promoter and showewd differential expression in the sub pallium knockout in mice

Unbaised Predicted Arx TFBS:
The bedfiles containing the unbiasedly predicted ARX Arx TFBS in the mouse and human genome for eahc motif model.
These were conducted in the mm9 and hg19 genome respsectively.

memesuite Output:
Memesuite output on NOn4mer containing Peaks




