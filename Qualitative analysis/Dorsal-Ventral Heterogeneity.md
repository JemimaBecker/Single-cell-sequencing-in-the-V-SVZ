# Dorsal-Ventral Heterogeneity in V-SVZ lncRNA expression

## Key results

- Note that not all of the genes could be annotated - I included genes with no annotated biotype, hence the inclusion of the histone variants in the plot.
- Importantly, many lncRNAs that are known to be regulators of neurogenesis (Meg3, Sox2ot, Pax6os1) show statistically significant variance in the distribution of expression within the niche
- This is consitent with a heterogenetiy in neurogenic capacity and distribution of NSCs of different classes within the niche

| Plot | Plot1 (zoomed in)|
|---|---|
| <img width="972" alt="Screenshot 2022-02-22 at 16 34 12" src="https://user-images.githubusercontent.com/67189202/155176574-f06272cc-a0e5-4b35-8915-96aa80ef7d03.png"> | <img width="971" alt="Screenshot 2022-02-22 at 16 35 36" src="https://user-images.githubusercontent.com/67189202/155176827-d675bcc1-ff7b-418e-9a3e-0b959367f39c.png"> |

## Papers

Cebrian-Silla, A., Nascimento, M.A., Redmond, S.A., Mansky, B., Wu, D., Obernier, K., Romero Rodriguez, R., Gonzalez-Granero, S., García-Verdugo, J.M., Lim, D.A., et al. (2021). Single-cell analysis of the ventricular-subventricular zone reveals signatures of dorsal and ventral adult neurogenesis. ELife 10, e67436.

## Code


load the data
```{r setup, include=FALSE}
library(dplyr)
library(readr)
CS2021Lineage_DE_genes_Table_1 <- read_csv("Cebrian-Silla et al 2021 /Supplementary file 5 Dorsal and ventral lineages DE genes and GO terms/Lineage DE genes-Table 1.csv")
CS2021Lineage_DE_genes_Table_1 <- CS2021Lineage_DE_genes_Table_1 %>% dplyr::filter(CS2021Lineage_DE_genes_Table_1$p_val_adj < 0.05)
CS_2021_DV_DE <- CS2021Lineage_DE_genes_Table_1[,c(1,3,6)]
colnames(CS_2021_DV_DE) <- c("external_gene_name","log2FC","CS 2021 Padj D-V DE")
```
Load ensembl
```{r}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)
ensembl    <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
                                 'chromosome_name', 'start_position', 'end_position', 
                                 'strand', 'gene_biotype','transcript_count'), 
                    mart = ensembl)

```
Annotate the genes and filter to biotypes of interest
```{r}
CS_2021_DV_DE <- merge(CS_2021_DV_DE,ensEMBL2id,
by.x="external_gene_name",
by.y="external_gene_name",
all.x=TRUE,
all.y=FALSE)
CS_2021_DV_DE$gene_biotype <- as.factor(CS_2021_DV_DE$gene_biotype)
DV_DE_ncRNA <- CS_2021_DV_DE %>% dplyr::filter(CS_2021_DV_DE$gene_biotype == "lncRNA"|                      gene_biotype == "processed_pseudogene"|        
          gene_biotype == "transcribed_processed_pseudogene"|
          is.na(CS_2021_DV_DE$gene_biotype))
```
volcano plot
```{r}
 if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')


  BiocManager::install('EnhancedVolcano')
  devtools::install_github('kevinblighe/EnhancedVolcano')
```
```{r}
  library(EnhancedVolcano)
  
plot <-   EnhancedVolcano(DV_DE_ncRNA,
    lab = DV_DE_ncRNA$external_gene_name,
    title = "Dorsal-Ventral differences in lncRNA expression",
    x = 'log2FC',
    y = 'CS 2021 Padj D-V DE',
    xlim = c(-0.8,1.6),
    ylim=c(0,175),
    pCutoff = 0.05,
  FCcutoff = 0.1,
  labSize = 4.0,
    legendPosition = 'right',
    legendLabSize = 12,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.75)
    
plot1 <-   EnhancedVolcano(DV_DE_ncRNA,
    lab = DV_DE_ncRNA$external_gene_name,
    title = "Dorsal-Ventral differences in lncRNA expression",
    x = 'log2FC',
    y = 'CS 2021 Padj D-V DE',
    xlim = c(-0.8,0.75),
    ylim=c(0,75),
    pCutoff = 0.05,
  FCcutoff = 0.1,
  labSize = 4.0,
    legendPosition = 'right',
    legendLabSize = 12,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.75)
```
