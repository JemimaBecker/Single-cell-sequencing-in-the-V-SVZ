# Stem cell activation

Using the following papers:
Mizrak et al., 2020
Zywitza et al., 2018

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

### Mizrak et al 2020

Table S1. Binomial Test Results Showing the Cell-Type Specificity of Genes in OB Neuron Subtypes and V-SVZ Neuronal Lineage Clusters in Different Samples, Related to Figures 1, 2, 3, and S1.
need to do some more formatting here to isolate the specific datasets - many many tabs on this spreadsheet
just take the sheets that represetn aNSCs

```{r}
mizrak_supp1 <- read_excel("mizrak et al 2020 supp1.xlsx")
excel_sheets("mizrak et al 2020 supp1.xlsx")
mizrak_aNSC_NESFLPO <- read_excel("mizrak et al 2020 supp1.xlsx", sheet = "NESFLPO_aNSC"  )
mizrak_aNSC_GCERT <- read_excel("mizrak et al 2020 supp1.xlsx", sheet = "GCERT2_aNSC"   )
mizrak_aNSC_SCOPE <- read_excel("mizrak et al 2020 supp1.xlsx", sheet =  "SCOPE_aNSC"    )

colnames(mizrak_aNSC_NESFLPO)[2] <- "external_gene_name"
colnames(mizrak_aNSC_GCERT)[2] <- "external_gene_name"
colnames(mizrak_aNSC_SCOPE)[2] <- "external_gene_name"

colnames(mizrak_aNSC_NESFLPO)[3] <- "Mizrak NESFLPO FDR"
colnames(mizrak_aNSC_GCERT)[3] <- "Mizrak GCERT FDR"
colnames(mizrak_aNSC_SCOPE)[3] <- "Mizrak SCOPE FDR"

colnames(mizrak_aNSC_NESFLPO)[4] <- "Mizrak NESFLPO Log2effect"
colnames(mizrak_aNSC_GCERT)[4] <- "Mizrak GCERT Log2effect"
colnames(mizrak_aNSC_SCOPE)[4] <- "Mizrak SCOPE Log2effect"

mizrak_aNSC_NESFLPO <- mizrak_aNSC_NESFLPO[,c(2:4)]
mizrak_aNSC_GCERT<-mizrak_aNSC_GCERT[,c(2:4)]
mizrak_aNSC_SCOPE<-mizrak_aNSC_SCOPE[,c(2:4)]

mizrak_aNSC_NESFLPO <- mizrak_aNSC_NESFLPO %>% dplyr::filter(mizrak_aNSC_NESFLPO$`Mizrak NESFLPO FDR` < 0.05)
mizrak_aNSC_GCERT <- mizrak_aNSC_GCERT %>% dplyr::filter(mizrak_aNSC_GCERT$`Mizrak GCERT FDR` < 0.05)
mizrak_aNSC_SCOPE <- mizrak_aNSC_SCOPE %>% dplyr::filter(mizrak_aNSC_SCOPE$`Mizrak SCOPE FDR` < 0.05)
```

### Zywitza et al 2018

load sheets of the table individually to see sheets for aNSCs and qNSCs

```{r}
scl8_aNSCs_Table_1 <- read_csv("Zywitza et al 2018/mmc4 (1)/scl8_aNSCs-Table 1.csv")
scl7_qNSCs_Table_1 <- read_csv("Zywitza et al 2018/mmc4 (1)/scl7_qNSCs-Table 1.csv")
cl17_NSCs_Table_1 <- read_csv("Zywitza et al 2018/mmc2/cl17_NSCs-Table 1.csv")
colnames(scl8_aNSCs_Table_1)[6] <- "cluster or subcluster"
colnames(scl7_qNSCs_Table_1)[6] <- "cluster or subcluster"
colnames(cl17_NSCs_Table_1)[6] <- "cluster or subcluster"
zywitza_NSCs <- rbind(scl8_aNSCs_Table_1,scl7_qNSCs_Table_1,cl17_NSCs_Table_1)
colnames(zywitza_NSCs) <- c("external_gene_name","Zywitza Pval","Zywitza avg_diff","Zywitza Pct1","Zywitza Pct2","Zywitza cluster or subcluster")

zywitza_aNSCs <- scl8_aNSCs_Table_1
colnames(zywitza_aNSCs)[c(1:3)] <- c("external_gene_name","Zywitza p value","Zywitza avg diff")
zywitza_aNSCs <- zywitza_aNSCs[,c(1:3)]
zywitza_aNSCs <- zywitza_aNSCs %>% dplyr::filter(zywitza_aNSCs$`Zywitza p value` < 0.05)
```
merge the fules for different datasets by gene name
```{r}
merge_activation <- merge(zywitza_aNSCs, mizrak_aNSC_NESFLPO,
                          by.x="external_gene_name",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=TRUE)
merge_activation <- merge(merge_activation, mizrak_aNSC_GCERT,
                          by.x="external_gene_name",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=TRUE)
merge_activation <- merge(merge_activation, mizrak_aNSC_SCOPE,
                          by.x="external_gene_name",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=TRUE)

merge_activation <- merge(merge_activation, ensEMBL2id,
                          by.x="external_gene_name",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=FALSE)
merge_activation$gene_biotype <- as.factor(merge_activation$gene_biotype)
merge_activation_ncRNA <- merge_activation %>% dplyr::filter(merge_activation$gene_biotype=="lncRNA"                             
|gene_biotype=="miRNA"                              
|gene_biotype=="misc_RNA"                          
|gene_biotype=="polymorphic_pseudogene"            
|gene_biotype=="processed_pseudogene"                         
|gene_biotype=="pseudogene"                           
|gene_biotype=="TEC"                             
|gene_biotype=="transcribed_processed_pseudogene"
|gene_biotype=="transcribed_unitary_pseudogene"     
|gene_biotype=="transcribed_unprocessed_pseudogene" 
|gene_biotype=="unitary_pseudogene"                 
|gene_biotype=="unprocessed_pseudogene"    
|is.na(merge_activation$gene_biotype))
merge_activation_p <- merge_activation_ncRNA[,c(1,2,4,6,8,11:16)]
merge_activation_p$`Mizrak NESFLPO FDR` <- as.numeric(merge_activation_p$`Mizrak NESFLPO FDR`)
```
calculate Cumulative Occurences (CO) by summing the number of timnes that a lncRNA has a results with p<0.05. All p values greater than 0.05 are converted to NAs, and then 0; all p values less than 0.05 are converted to 1 and then summed accros the row.
```{r}
rank_calc <- merge_activation_p
names <- rank_calc[,c(1,6:11)]
rank_calc <- rank_calc[,-c(1,6:11)]
rank_calc[!is.na(rank_calc)] <- 1
rank_calc[is.na(rank_calc)] <- 0
ranks <- cbind(names, rank_calc)
sapply(ranks, class)
chars <- sapply(rank_calc, is.character)
rank_calc[ , chars] <- as.data.frame(apply(rank_calc[ , chars], 2, as.numeric))
activation_summed_ranks <- cbind(names,
                      rowSums(rank_calc),
                      rank_calc)
```
```{r}
write.csv(activation_summed_ranks,"OUTPUT_activation.csv")
```
## Citations

- Howe, K.L., Achuthan, P., Allen, J., Allen, J., Alvarez-Jarreta, J., Amode, M.R., Armean, I.M., Azov, A.G., Bennett, R., Bhai, J., et al. (2021). Ensembl 2021. Nucleic Acids Research 49, D884–D891. https://doi.org/10.1093/nar/gkaa942.
- Mizrak, D., Bayin, N.S., Yuan, J., Liu, Z., Suciu, R.M., Niphakis, M.J., Ngo, N., Lum, K.M., Cravatt, B.F., Joyner, A.L., et al. (2020). Single-Cell Profiling and SCOPE-Seq Reveal Lineage Dynamics of Adult Ventricular-Subventricular Zone Neurogenesis and NOTUM as a Key Regulator. Cell Rep 31, 107805. https://doi.org/10.1016/j.celrep.2020.107805.
- Zywitza, V., Misios, A., Bunatyan, L., Willnow, T.E., and Rajewsky, N. (2018). Single-Cell Transcriptomics Characterizes Cell Types in the Subventricular Zone and Uncovers Molecular Defects Impairing Adult Neurogenesis. Cell Reports 25, 2457-2469.e8. https://doi.org/10.1016/j.celrep.2018.11.003.
