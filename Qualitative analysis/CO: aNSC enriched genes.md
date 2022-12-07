---
title: "findallmarkers aNSCs"
author: "Jemima Becker"
date: "2022-08-27"
output: html_document
---

This code combines datasets listing genes that have been found via FindAllMarkers to be enriched within aNSC clusters and calculates Cumulative Occurence (CO) as discussed in the paper.

1: Setup
```{r setup, include=FALSE}
library(readr)
library(readxl)
library(dplyr)
library(biomaRt)
library(tidyverse)
ensembl    <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
                                 'chromosome_name', 'start_position', 'end_position', 
                                 'strand', 'gene_biotype'), 
                    mart = ensembl)
```
2: Load in data and combine to a single dataframe
```{r}
FindAll_aNSCs <- as.data.frame(all_lncRNAs_merged[,c("external_gene_name","Del Aguila 2022 aNSCs TAPs","Mizrak 2019 aNSCs","Mizrak 2020 NESFPLO aNSC","Mizrak 2020 GCERT2 aNSC","Mizrak 2020 SCOPE aNSC")])
```

3: Calulate CO and annotate dataframe by gene names
```{r}
names <- FindAll_aNSCs[,1]
ranks <- FindAll_aNSCs[,-c(1)]
ranks[!is.na(ranks)] <- as.numeric(1)
ranks[is.na(ranks)] <- as.numeric(0)
ranks2 <- data.frame(apply(ranks, 2, function(x) as.numeric(as.character(x))))
sapply(ranks2, class)
FindAll_aNSCs_ranksum <- cbind(names,
                      rowSums(ranks2),
                      ranks2)

colnames(FindAll_aNSCs_ranksum)[2] <- "cumulative occurences FindAll_aNSCs"
FindAll_aNSCs_ranksum <- merge(FindAll_aNSCs_ranksum, ensEMBL2id,
                          by.x="names",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=FALSE)
write.csv(FindAll_aNSCs_ranksum,"FindAll_aNSCs_ranksum.csv")
```

Data used from:
del Águila, Á., Adam, M., Ullom, K., Shaw, N., Qin, S., Ehrman, J., Nardini, D., Salomone, J., Gebelein, B., Lu, Q.R., et al. (2022). Olig2 defines a subset of neural stem cells that produce specific olfactory bulb interneuron subtypes in the subventricular zone of adult mice. Development 149, dev200028. https://doi.org/10.1242/dev.200028.
Mizrak, D., Bayin, N.S., Yuan, J., Liu, Z., Suciu, R., Niphakis, M.J., Ngo, N., Lum, K.M., Cravatt, B.F., Joyner, A.L., et al. (2019a). Single-cell profiling and SCOPE-seq reveal the lineage dynamics of adult neurogenesis and NOTUM as a key V-SVZ regulator. BioRxiv 770610. https://doi.org/10.1101/770610.
Mizrak, D., Levitin, H.M., Delgado, A.C., Crotet, V., Yuan, J., Chaker, Z., Silva-Vargas, V., Sims, P.A., and Doetsch, F. (2019b). Single-Cell Analysis of Regional Differences in Adult V-SVZ Neural Stem Cell Lineages. Cell Reports 26, 394-406.e5. https://doi.org/10.1016/j.celrep.2018.12.044.
