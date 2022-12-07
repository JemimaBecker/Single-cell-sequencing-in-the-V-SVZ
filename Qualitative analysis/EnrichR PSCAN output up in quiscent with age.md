---
title: "Enrichr PSCAN output up in quiescent with age"
author: "Jemima Becker"
date: "2022-08-27"
output: html_document
---
```{r}

library(dplyr)
library(enrichR)
library(base)
library(tictoc)
```

PSCAN: 1317.53 sec elapsed


```{r}
tic("PSCAN")
PSCAN_output_up_age <- read_csv("up in quiescent PSCAN output.csv")


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
toc()
```

```{r}
final_significant_res$term_DB <- paste(final_significant_res$term, final_significant_res$database)
write.csv(final_significant_res,"final_significant_res.csv")
```

there are 43643 levels
```{r}

final_significant_res$term_DB <- as.factor(final_significant_res$term_DB)
summ <- final_significant_res %>% count(final_significant_res$term_DB)
all_summ <- summ %>% dplyr::filter(summ$n >= 9)
write.csv(all_summ,"all_summ.csv")
```
