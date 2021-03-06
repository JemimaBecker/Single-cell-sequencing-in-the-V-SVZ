# Data summary for scRNAseq meta-analysis 

## Fig 5: ScRNAseq data summary
| Figure 5: Data summary for scRNAseq meta-analysis. |
| --- |
| <img width="639" alt="Screenshot 2022-04-19 at 13 57 50" src="https://user-images.githubusercontent.com/67189202/164009216-30ee741b-264d-4d9c-8446-53f16fcbf2ef.png"> |
| A: Total number of studies in which each lncRNA is differentially expressed across all categories, data is expressed as log2(number of lncRNAs) vs the number of studies that each gene was differentially expressed in. |
| B: Overlap between lncRNAs differentially expressed with stem cell activation (n=750) and those that are downregulated with ageing (n=28). |
|C: EnrichR gene ontology results for transcription factors that bind the promoters of lncRNAs downregulated with age. EnrichR was run on PSCAN outputs to find ontologies of transcription factors that bind each promoter, the number of lncRNA promoters that were associated with each GO term is shown. |


Bar graphs were plotted in Numbers, and Venn diagram made in PowerPoint.

### Conduct PSCAN

- PCAN was conducted via the web server (http://159.149.160.88/pscan/) to examine which transcription factor binding motifs are enriched in the promoters of lncRNA that is downregulated with age.
- This was conducted in a window -450bp upstream and +50bp downstream of the transcriptional start site for each gene.
- If a transcription factor binding site was enriched in the selected sequence with P<0.05 relative to a scrambled control sequence it was treated as significant.

Enrichr was run individually for the transcritption factors that bound each promoter. 

### Conduct EnrichR
```{r}
#Load relevant packages
library(dplyr)
library(enrichR)
library(base)
PSCAN_output_up_age <- read_csv("PSCAN_output_up_age.csv") # load in the dataset
final_res <- cbind("1","2","3","4","5") # set a blank object to add EnrichR output into
  colnames(final_res) <- c("gene name","term","padj","genes","database")
  
# set up a for loop that runs this code for each column, h, of the spreadsheet PSCANPoutput_up_age. Each column name is the lncRNA whose promoter binding sites are being examined
  
for (h in colnames(PSCAN_output_up_age)){ 
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
    for (i in dbs$libraryName){ 
    
# establishmnet of a second, nested, for loop that loops through each individual database, i, examining the set of genes defined by the variable h.
    
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

### Plot

- The number of times that each individual GO term was significant was summed to examine whether there were any terms that were associated with multiple lncRNA promoters.
- This was plotted as a bar graph in Numbers.

## Example output

<img width="955" alt="Screenshot 2022-03-27 at 20 26 37" src="https://user-images.githubusercontent.com/67189202/160297444-ad11f8dd-103b-46f3-b69d-78919f4dd2c6.png">


## References

- Chen EY, Tan CM, Kou Y, Duan Q, Wang Z, Meirelles GV, Clark NR, Ma'ayan A. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics. 2013; 128(14).
- F.Zambelli, G.Pesole, G.Pavesi, Pscan: Finding Over-represented Transcription Factor Binding Site Motifs in Sequences from Co-Regulated or Co-Expressed Genes. Nucleic Acids Research 2009 37(Web Server issue):W247-W252.
- Kuleshov MV, Jones MR, Rouillard AD, Fernandez NF, Duan Q, Wang Z, Koplev S, Jenkins SL, Jagodnik KM, Lachmann A, McDermott MG, Monteiro CD, Gundersen GW, Ma'ayan A. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research. 2016; gkw377 .
- Xie Z, Bailey A, Kuleshov MV, Clarke DJB., Evangelista JE, Jenkins SL, Lachmann A, Wojciechowicz ML, Kropiwnicki E, Jagodnik KM, Jeon M, & Ma???ayan A. Gene set knowledge discovery with Enrichr. Current Protocols, 1, e90. 2021. doi: 10.1002/cpz1.90
