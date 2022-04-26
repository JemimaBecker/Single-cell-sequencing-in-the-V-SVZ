# Sox2ot case study

## Figure 6: Sox2ot case study
|Figure 6: Sox2ot case study|
| --- |
| <img width="686" alt="Screenshot 2022-04-19 at 15 27 33" src="https://user-images.githubusercontent.com/67189202/164027237-f680add8-17f8-46da-86ae-012bf713197b.png"> |
| A: Expression of Sox2ot in mouse brain, determined by ISH. Red=high, green=low expression. Data from Allen Brain Atlas. |
| B: Expression of Sox2ot in different cell types. Data from Mouse Organogenesis Cell Atlas. |
| C: EnrichR output for WikiPathways 2019 Mouse (left) and WikiPathways 2021 Human (right), for genes with overrepresented binding sites in the mouse or human Sox2ot promoters. Binding site identified by PSCAN p<0.05. |

## Binding site enrichment

- Enriched binding sites were identified by running PSCAN via webserver (http://159.149.160.88/pscan/) to find transcription factors whose binding sites were overrepresented in a window 450bp upstream and 50bp downstream of the transcription start site for the Human SOX2-OT and Mouse Sox2ot.
- EnrichR was run on transcription factors with overrepresented binding site motifs (p<0.05).

1: Install packages
```{r}
library(devtools)
library(readr)
install_github("wjawaid/enrichR")
install.packages("enrichR")
library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)
```
2: Load data
- This can also be done by reading in a .csv file of the gene names
```{r}
query_genes<-c("OVOL2"	,"HOXA10"	,"CDX4"	,"ZNF684"	,"HOXD9"	,"Rara"	,"CDX1"	,"GATA1"	,"SOX10"	,"FOXE1"	,"RARA"	,"PITX2"	,"OTX2"	,"CDX2"	,"SOX13"	,"SOX9"	,"SOX18"	,"PROX1"	,"SRF"	,"SOX2"	,"Zic1"	,"PITX1"	,"HOXC12"	,"EWSR1-FLI1"	,"Sox6"	,"HNF4A"	,"Foxd3"	,"RXRA"	,"Sox3"	,"Nfe2l2"	,"IRF9"	,"HOXD10"	,"ZNF384"	,"HOXC10"	,"CLOCK"	,"Zic2"	,"RARA")
if (websiteLive) {
  enriched <- enrichr(query_genes, dbs)
}
#this returns all of the potential tables as a massive object - can select individually or use code below to put them all together
names <- dbs$libraryName
```
Add column to specify the database that each term is from
```{r}
for (i in names){
  enriched[[i]] <- cbind(enriched[[i]], enriched$i)
  colnames(enriched[[i]])[10] <- "database"
}
```
Collate all of the results to a single data frame:
```{r}
enrichment_results <- as.data.frame(if (websiteLive) enriched[["Genome_Browser_PWMs"]])
for (i in names)
  {
  enrichment_results <- rbind(enrichment_results, as.data.frame(if (websiteLive) enriched[[i]]))}
 ```
Filter to significant results and plot
```{r}
significant_enrichment_results <- subset(enrichment_results, "Adjusted p value" < 0.05)
visualisation for individual results e.g.:
dbs <- c("WikiPathways_2019_Mouse" )
if (websiteLive) {
  enriched <- enrichr(query_genes, dbs)
}
if (websiteLive) plotEnrich(enriched[["WikiPathways_2019_Mouse"]], 
                            showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",
                            xlab = "Enriched terms (WikiPathways_2019_Mouse)",
                            ylab = "Gene count",
                            title = "Enrichment analysis for Sox2ot binding motifs in Mouse"
)
```
Repeat for human
```{r}
query_genes<-c("POU4F2"	,"ARGFX"	,"CUX1"	,"ALX3"	,"CUX2"	,"LBX1"	,"LHX9"	,"mix-a"	,"PRRX2"	,"PAX4"	,"HSF2"	,"PPARA","RXRA"	,"LHX5"	,"ZNF140"	,"POU4F3"	,"GSX1"	,"VSX2"	,"MIXL1"	,"NR4A1"	,"LMX1A"	,"SOX21"	,"PAX7"	,"Alx4"	,"Lhx4"	,"NKX6-1"	,"EMX1"	,"MNX1"	,"ISX"	,"RAX2"	,"UNCX"	,"MYB"	,"PAX3"	,"Sox11"	,"LMX1B"	,"Shox2"	,"PRRX1"	,"VSX1"		,"Lhx8"	,"Dmbx1"	,"LBX2"	,"HOXD13"	,"HNF1B"	,"SHOX"	,"MEOX2"	,"HOXA2"	,"Hmx3"	,"MEOX1"	,"NKX6-2"	,"VAX2"	,"TLX2"	,"DRGX"	,"VAX1"	,"ESX1"	,"EN1"	,"DLX6"	,"Alx1"	,"NOTO"	,"Twist2"	,"GSX2"	,"Hmx2"	,"ZNF652"	,"BACH2"	,"HSF1"	,"Hmx1"	,"HOXD3")
if (websiteLive) {
  enriched <- enrichr(query_genes, dbs)
}
#this returns all of the potential tables as a massive object - can select individually
dbs <- c("WikiPathway_2021_Human" )
if (websiteLive) {
  enriched <- enrichr(query_genes, dbs)
}
if (websiteLive) plotEnrich(enriched[["WikiPathway_2021_Human"]], 
                            showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",
                            xlab = "Enriched terms (WikiPathway_2021_Human)",
                            ylab = "Gene count",
                            title = "Enrichment analysis for SOX2OT binding motifs in Human"
)
```

## References

- Allen Institute for Brain Science (2004). Allen Mouse Brain Atlas [dataset]. Available from mouse.brain-map.org.Allen Institute for Brain Science (2011).
- Cao, J., Spielmann, M., Qiu, X. et al. The single-cell transcriptional landscape of mammalian organogenesis. Nature 566, 496–502 (2019). https://doi.org/10.1038/s41586-019-0969-x
- Chen EY, Tan CM, Kou Y, Duan Q, Wang Z, Meirelles GV, Clark NR, Ma'ayan A. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics. 2013; 128(14).
- Lein, E., Hawrylycz, M., Ao, N. et al. Genome-wide atlas of gene expression in the adult mouse brain. Nature 445, 168–176 (2007). https://doi.org/10.1038/nature05453
- Kuleshov MV, Jones MR, Rouillard AD, Fernandez NF, Duan Q, Wang Z, Koplev S, Jenkins SL, Jagodnik KM, Lachmann A, McDermott MG, Monteiro CD, Gundersen GW, Ma'ayan A. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research. 2016; gkw377 .
- Xie Z, Bailey A, Kuleshov MV, Clarke DJB., Evangelista JE, Jenkins SL, Lachmann A, Wojciechowicz ML, Kropiwnicki E, Jagodnik KM, Jeon M, & Ma’ayan A. Gene set knowledge discovery with Enrichr. Current Protocols, 1, e90. 2021. doi: 10.1002/cpz1.90 
- F.Zambelli, G.Pesole, G.Pavesi Pscan: Finding Over-represented Transcription Factor Binding Site Motifs in Sequences from Co-Regulated or Co-Expressed Genes. Nucleic Acids Research 2009 37(Web Server issue):W247-W252.
