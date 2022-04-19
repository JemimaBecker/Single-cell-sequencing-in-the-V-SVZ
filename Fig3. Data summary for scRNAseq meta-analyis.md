# Fig 3: Data summary for scRNAseq meta-analysis 

<img width="639" alt="Screenshot 2022-04-19 at 13 57 50" src="https://user-images.githubusercontent.com/67189202/164009216-30ee741b-264d-4d9c-8446-53f16fcbf2ef.png">

> Figure 3: Data summary for scRNAseq meta-analysis. 

> A: Total number of studies in which each lncRNA is differentially expressed across all categories, data is expressed as log2(number of lncRNAs) vs the number of studies that each gene was differentially expressed in. 

> B: Overlap between lncRNAs differentially expressed with stem cell activation (n=750) and those that are downregulated with ageing (n=28). 

> C: EnrichR gene ontology results for transcription factors that bind the promoters of lncRNAs downregulated with age. EnrichR was run on PSCAN outputs to find ontologies of transcription factors that bind each promoter, the number of lncRNA promoters that were associated with each GO term is shown.

Bar graphs were plotted in Numbers, and Venn diagram made in PowerPoint.

## Fig3C step 1: Conduct PSCAN

PCAN was conducted via the web server (http://159.149.160.88/pscan/) to examine which transcription factor binding motifs are enriched in the promoters of lncRNAs that are downregulated with age.
This was conducted in a window -450bp upstream and +50bp downstream of the transcriptional start site for each gene.
A transcription factor binding site was treated as "enriched" with a PSCAN p<0.05.

## Fig3C step 2: Conduct GO via Enrichr

Enrichr was run individually for the transcritption factors that bound each promoter. 

```{r}
library(dplyr)
library(enrichR)
library(base)
```
```{r}
PSCAN_output_up_age <- read_csv("PSCAN_output_up_age.csv") # load in the dataset
final_res <- cbind("1","2","3","4","5") # set a blank object to add EnrichR output into
  colnames(final_res) <- c("gene name","term","padj","genes","database")

for (h in colnames(PSCAN_output_up_age)){ # set up a for loop that runs this code for each column, h, of the spreadsheet PSCANPoutput_up_age. Each column name is the lncRNA whose promoter binding sites are being examined
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

    for (i in dbs$libraryName){ # establishmnet of a second, nested, for loop that loops through each individual database, i, examining the set of genes    defined by the variable h.
        if (websiteLive) {
        enriched <- enrichr(query_genes, i)}
        enrichment_res <- as.data.frame(enriched)
        if (nrow(enrichment_res)<1){} # skip over if there are no results for the gene set in the queried database - otherwise the code is terminated
        else {
        results <- cbind(enrichment_res[,c(1,4,9)], paste0(i))
        colnames(results)[c(1:4)] <- c("term","padj","genes","database")
        all_res <- rbind(all_res, results) 
        }

      }
  all_res$database <- as.factor(all_res$database)
  summary(all_res$database)
  allres_databases <- levels(all_res$database)
  ind_res <- cbind(paste0(h),all_res)
  colnames(ind_res)[1] <- "gene name"
  final_res <- rbind(final_res,ind_res)
}

# save files and filter
final_significant_res <- final_res
final_significant_res$padj <- as.numeric(final_significant_res$padj) # format to numeric so that it can be filtered using numeric criteria
final_significant_res <- final_significant_res %>% dplyr::filter(final_significant_res$padj < 0.05) # filter to significant results
```
## Fig3C step 3: Plot

The number of times that each individual GO term was significant was summed to examine whether there were any terms that were associated with multiple lncRNA promoters.
This was plotted as a bar graph in Numbers.
