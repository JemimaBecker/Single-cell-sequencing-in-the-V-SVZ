# Sox2ot Case study

1: Setup
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
```{r}
#enter your data - in this case the output genes from PSCAN looking at the Sox2ot promoter
query_genes<-c("OVOL2","HOXA10","CDX4","ZNF684","HOXD9","Rara","CDX1","GATA1","SOX10","FOXE1","RARA","PITX2","OTX2","CDX2","SOX13","SOX9","SOX18","PROX1","SRF","SOX2"	,"Zic1","PITX1","HOXC12","EWSR1-FLI1","Sox6","HNF4A","Foxd3","RXRA","Sox3","Nfe2l2","IRF9","HOXD10","ZNF384","HOXC10","CLOCK","Zic2","RARA")
```

```{r}
if (websiteLive) {
  enriched <- enrichr(query_genes, dbs)
}
#this returns all of the potential tables as a massive object - can select individually or use code below to put them all together
names <- dbs$libraryName
```
3: Add column to specify the database that each term is from
```{r}
for (i in names){
  enriched[[i]] <- cbind(enriched[[i]], enriched$i)
  colnames(enriched[[i]])[10] <- "database"
}
# collate all of the results to a single data frame:

enrichment_results <- as.data.frame(if (websiteLive) enriched[["Genome_Browser_PWMs"]])
for (i in names)
  {
  enrichment_results <- rbind(enrichment_results, as.data.frame(if (websiteLive) enriched[[i]]))}
# filter to significant results !! need to edit name of col, cant do atm as enrichr is down so cant make object

significant_enrichment_results <- subset(enrichment_results, "Adjusted p value" < 0.05)
visualisation for individual results e.g.:
```{r}
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
Repeat with PSCAN results from human
```{r}
query_genes<-c("POU4F2","ARGFX","CUX1","ALX3","CUX2","LBX1","LHX9","mix-a","PRRX2","PAX4","HSF2","PPARA","RXRA","LHX5","ZNF140","POU4F3","GSX1","VSX2","MIXL1","NR4A1","LMX1A","SOX21","PAX7","Alx4","Lhx4","NKX6-1","EMX1","MNX1","ISX","RAX2","UNCX","MYB","PAX3","Sox11","LMX1B","Shox2","PRRX1","VSX1","Lhx8","Dmbx1","LBX2","HOXD13","HNF1B","SHOX","MEOX2","HOXA2","Hmx3","MEOX1","NKX6-2","VAX2","TLX2","DRGX","VAX1","ESX1","EN1","DLX6","Alx1","NOTO","Twist2","GSX2","Hmx2","ZNF652","BACH2","HSF1","Hmx1","HOXD3")
if (websiteLive) {
  enriched <- enrichr(query_genes, dbs)
}
#this returns all of the potential tables as a massive object - can select individually or use code below to put them all together
#names <- dbs$libraryName

dbs <- c("BioPlanet_2019" )
if (websiteLive) {
  enriched <- enrichr(query_genes, dbs)
}

if (websiteLive) plotEnrich(enriched[["BioPlanet_2019"]], 
                            showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",
                            xlab = "Enriched terms (BioPlanet_2019)",
                            ylab = "Gene count",
                            title = "Enrichment analysis for SOX2OT binding motifs in Human"
)

dbs <- c("CellMarker_Augmented_2021" )
if (websiteLive) {
  enriched <- enrichr(query_genes, dbs)
}

if (websiteLive) plotEnrich(enriched[["CellMarker_Augmented_2021"]], 
                            showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",
                            xlab = "Enriched terms (CellMarker_Augmented_2021)",
                            ylab = "Gene count",
                            title = "Enrichment analysis for SOX2OT binding motifs in Human"
)

dbs <- c("PanglaoDB_Augmented_2021" )
if (websiteLive) {
  enriched <- enrichr(query_genes, dbs)
}

if (websiteLive) plotEnrich(enriched[["PanglaoDB_Augmented_2021"]], 
                            showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",
                            xlab = "Enriched terms (PanglaoDB_Augmented_2021)",
                            ylab = "Gene count",
                            title = "Enrichment analysis for SOX2OT binding motifs in Human"
)

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

dbs <- c("WikiPathways_2019_Mouse" )
if (websiteLive) {
  enriched <- enrichr(query_genes, dbs)
}

if (websiteLive) plotEnrich(enriched[["WikiPathways_2019_Mouse"]], 
                            showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",
                            xlab = "Enriched terms (WikiPathways_2019_Mouse)",
                            ylab = "Gene count",
                            title = "Enrichment analysis for SOX2OT binding motifs in Human"
)

```
