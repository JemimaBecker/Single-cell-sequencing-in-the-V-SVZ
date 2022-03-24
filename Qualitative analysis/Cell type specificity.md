---
title: "Lineage and cell types"
author: "Jemima Becker"
date: "24/03/2022"
output: html_document
---
Load relevant packages
```{r}
library(readr)
library(biomaRt)
library(readxl)
library(dplyr)
ensembl    <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
                                 'chromosome_name', 'start_position', 'end_position', 
                                 'strand', 'gene_biotype'), 
                    mart = ensembl)
```

# Del Aguila 2022
using supplementary file 1; cell cluster enriched genes at different stagess

```{r}
delA_supp1 <- read_csv("supp1 del aguila 2022.csv")
colnames(delA_supp1) <- delA_supp1[2,]
delA_supp1 <- delA_supp1[-c(1,2),]
colnames(delA_supp1)[1] <- "external_gene_name"
excel_sheets("Del Aguila 2022/supp1 del aguila 2022.xlsx")

delA_supp1_E12.5 <-read_excel("Del Aguila 2022/supp1 del aguila 2022.xlsx", 
                       sheet = "Suppl Table 1-1 E12.5")
delA_supp1_E18.5 <-read_excel("Del Aguila 2022/supp1 del aguila 2022.xlsx",
                       sheet = "Suppl Table 1-2 E18.5")
delA_supp1_P14 <-read_excel("Del Aguila 2022/supp1 del aguila 2022.xlsx",
                       sheet = "Suppl Table 1-3 P14")
delA_supp1_Adult <-read_excel("Del Aguila 2022/supp1 del aguila 2022.xlsx",
                       sheet = "Suppl Table 1-4 Adult")

colnames(delA_supp1_Adult) <- delA_supp1_Adult[2,]
delA_supp1_Adult <- delA_supp1_Adult[-c(1,2),]
colnames(delA_supp1_Adult)[1] <- "external_gene_name"
delA_pvalues_celltypes <- delA_supp1_Adult[c(1,6,7)]
colnames(delA_pvalues_celltypes)[c(2:3)] <- c("delA padj","delA cluster")
delA_pvalues_celltypes <- delA_pvalues_celltypes %>% dplyr::filter(delA_pvalues_celltypes$`delA padj`< 0.05)
```

# Fu et al 2021
Table S1. Differentially expressed genes for each cluster, related to Figure 2

```{R}
fu_supp1 <- read_excel("fu supp1.xlsx")
excel_sheets("fu supp1.xlsx")
fu_supp1_data <-read_excel("fu supp1.xlsx",
                       sheet = "Table S1")
colnames(fu_supp1_data)[c(1:8)] <- c("external_gene_name" ,"Fu p_val","Fu avg_logFC", "Fu pct.1","Fu pct.2","Fu p_value_adjust","Fu cluster","Fu gene") 
fu_p_values_celltypes <- fu_supp1_data[,c(1,6,7)]
fu_p_values_celltypes <- fu_p_values_celltypes %>% dplyr::filter(fu_p_values_celltypes$`Fu p_value_adjust`<0.05)
```

```{r}
merge_cell_types <- merge(fu_p_values_celltypes, delA_pvalues_celltypes,
                          by.x="external_gene_name",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=TRUE)


merge_cell_types <- merge(merge_cell_types, ensEMBL2id,
                          by.x="external_gene_name",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=FALSE)
merge_cell_types$gene_biotype <- as.factor(merge_cell_types$gene_biotype)

merge_cell_types_ncRNA <- merge_cell_types %>% dplyr::filter(merge_cell_types$gene_biotype=="lncRNA"                             |gene_biotype=="polymorphic_pseudogene"        
 |gene_biotype== "processed_pseudogene"     
 |gene_biotype=="pseudogene"                    
 |gene_biotype== "TEC"                
 |gene_biotype=="transcribed_unitary_pseudogene"
 |gene_biotype=="unprocessed_pseudogene"        
|is.na(merge_cell_types$gene_biotype))
merge_cell_types_ncRNA$`delA padj` <- as.numeric(merge_cell_types_ncRNA$`delA padj`)
```


```{r}
rank_calc <- merge_cell_types_ncRNA
names <- rank_calc[,-c(2,4)]
rank_calc <- rank_calc[,c(2,4)]
rank_calc[!is.na(rank_calc)] <- 1
rank_calc[is.na(rank_calc)] <- 0
ranks <- cbind(names, rank_calc)
sapply(ranks, class)
chars <- sapply(rank_calc, is.character)
rank_calc[ , chars] <- as.data.frame(apply(rank_calc[ , chars], 2, as.numeric))
celltype_summed_ranks <- cbind(names,
                      rowSums(rank_calc),
                      rank_calc)
```
```{r}
write.csv(celltype_summed_ranks,"OUTPUT_celltype_ranks.csv")
```
