# Analysis of binding site enrichment

This code will conduct Gene Ontology Analysis via EnrichR in R, using a spreadsheet where each column is an individual list of genes to query.

Input: Transcription factors whose binding motifs are overrepresented in the promoters of selected genes (relative to random control sequences) as identified using PSCAN in which p<0.05

The data is arranged such that each column is the list of genes which have over represented bidning sites in the selected promoter. the name of the column corresponds to the promoter under examination. E.g. column 1 entitled "Xist" is the list of transcription factors that are likely to bind to the Xist promoter.

This code loops through the whole dataset, each column corresponds to the enriched TFs for each gene promoter

The code takes each column individually and searches these genes in all 192 databases available via EnrichR, returning a final object final_significant_res that combines all of the significant results.


Time taken to run is dependent on internet connection - this took 57 mins 12 s running via Eduroam

## Code

### Required packages

```
library(dplyr)
library(enrichR)
library(base)
```

### Perform the analysis

```{r}
PSCAN_output_up_age <- read_csv("PSCAN_output_up_age.csv") # load in your dataset
final_res <- cbind("1","2","3","4","5") # set a blank object to add EnrichR output into
  colnames(final_res) <- c("gene name","term","padj","genes","database")

for (h in colnames(PSCAN_output_up_age)){ # set up a for loop that runs this code for each column, h, of the spreadsheet PSCANPoutput_up_age
  query_genes <- PSCAN_output_up_age[[h]]
  
listEnrichrSites()
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("Genome_Browser_PWMs" ) # generates a starting object from the results of analysis in the first database
if (websiteLive) {
  enriched <- enrichr(query_genes, dbs)
  enrichment_res <- as.data.frame(enriched)
  results <- cbind(enrichment_res[,c(1,4,9)],
                   "Genome_Browser_PWMs")
  colnames(results)[c(1:4)] <- c("term","padj","genes","database")
  all_res <- results
}
listEnrichrSites()
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)
#query_genes <- PSCAN_output_up_age$Xist

for (i in dbs$libraryName){ # establishmnet of a second, nested, for loop that loops through each individual database, i, examining the set of genes defined by the variable h.
  if (websiteLive) {
  enriched <- enrichr(query_genes, i)}
  enrichment_res <- as.data.frame(enriched)
  
  if (nrow(enrichment_res)<1){} # skip over if there are no results for the gene set in the queried database - otherwise the code is terminated
  else {
  results <- cbind(enrichment_res[,c(1,4,9)], paste0(i))
  colnames(results)[c(1:4)] <- c("term","padj","genes","database")
  all_res <- rbind(all_res, results) }

}
all_res$database <- as.factor(all_res$database)
summary(all_res$database)
allres_databases <- levels(all_res$database)
ind_res <- cbind(paste0(h),all_res)
colnames(ind_res)[1] <- "gene name"
final_res <- rbind(final_res,ind_res)
}
  
final_significant_res <- final_res
final_significant_res$padj <- as.numeric(final_significant_res$padj) # format to numeric so that it can be filtered using numeric criteria
final_significant_res <- final_significant_res %>% dplyr::filter(final_significant_res$padj < 0.05) # filter to significant results
```

## Example output

<img width="955" alt="Screenshot 2022-03-27 at 20 26 37" src="https://user-images.githubusercontent.com/67189202/160297444-ad11f8dd-103b-46f3-b69d-78919f4dd2c6.png">

## References

Chen EY, Tan CM, Kou Y, Duan Q, Wang Z, Meirelles GV, Clark NR, Ma'ayan A.
Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics. 2013; 128(14).

Kuleshov MV, Jones MR, Rouillard AD, Fernandez NF, Duan Q, Wang Z, Koplev S, Jenkins SL, Jagodnik KM, Lachmann A, McDermott MG, Monteiro CD, Gundersen GW, Ma'ayan A.
Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research. 2016; gkw377 .

Xie Z, Bailey A, Kuleshov MV, Clarke DJB., Evangelista JE, Jenkins SL, Lachmann A, Wojciechowicz ML, Kropiwnicki E, Jagodnik KM, Jeon M, & Ma’ayan A.
Gene set knowledge discovery with Enrichr. Current Protocols, 1, e90. 2021. doi: 10.1002/cpz1.90
