# Stem cell ageing


## Papers

Borett et al., 2020
Xie et al., 2020

## Code 

Load data and required packages

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

### Borett et al., 2020 
Supplementary file 6; DE in adult dNSC vs E14
("UP" = increases with age, "down" = decreases with age)
```{r}
borrett_supp6 <- read_excel("Borett 2020/1-s2.0-S221112472031007X-mmc7.xlsx")
borrett_supp6 <- cbind(borrett_supp6[,1], "Borrett 2020 UP", borrett_supp6[,5],"Borrett 2020 DOWN")
colnames(borrett_supp6) <- c("external_gene_name","Borett Comparison","external_gene_name","Borett Comparison")
borrett_supp6 <- borrett_supp6[-c(1:4),]
borrett_UP <- borrett_supp6[,c(1,2)]
borrett_down <- borrett_supp6[,c(3,4)]
borrett_supp6 <- rbind(borrett_supp6[,c(1,2)],borrett_supp6[,c(3,4)])
```
### Xie et al 2020

Xie Up_DEGs_2M_H1_H2 = down with age
Xie Up_DEGs_12M_H1_H2 = up with age
Xie Up_DEGs_2M_H3_L1 = down with age
Xie Up_DEGs_12M_H3_L1 = up with age

```{r setup, include=FALSE}
library(dplyr)
library(readr)
Sheet_1_TS3E_Up_DEGs_2M_H1_H2_Fig_7G <- read_csv("Xie et al 2020/TS3/Sheet 1-TS3E Up DEGs 2M_H1-H2, Fig.7G.csv")
Sheet_1_TS3F_Up_DEGs_12M_H1_H2_Fig_7G <- read_csv("Xie et al 2020/TS3/Sheet 1-TS3F. Up DEGs 12M_H1-H2, Fig.7G.csv")
Sheet_1_TS3H_Up_DEGs_2M_H3_L1_Fig_7G <- read_csv("Xie et al 2020/TS3/Sheet 1-TS3H. Up DEGs 2M H3-L1, Fig.7G.csv")
Sheet_1_TS3I_Up_DEGs_12M_H3_L1_Fig_7G <- read_csv("Xie et al 2020/TS3/Sheet 1-TS3I. Up DEGs 12M H3-L1, Fig.7G.csv")

colnames(Sheet_1_TS3E_Up_DEGs_2M_H1_H2_Fig_7G)[c(2,3)] <- c("Xie Up_DEGs_2M_H1_H2 log2FC","Xie Up_DEGs_2M_H1_H2 FDR p")
colnames(Sheet_1_TS3F_Up_DEGs_12M_H1_H2_Fig_7G)[c(2,3)] <- c("Xie Up_DEGs_12M_H1_H2 log2FC","Xie Up_DEGs_12M_H1_H2 FDR p")
colnames(Sheet_1_TS3H_Up_DEGs_2M_H3_L1_Fig_7G)[c(2,3)] <- c("Xie Up_DEGs_2M_H3_L1 log2FC","Xie Up_DEGs_2M_H3_L1 FDR p") 
colnames(Sheet_1_TS3I_Up_DEGs_12M_H3_L1_Fig_7G)[c(2,3)] <- c("Xie Up_DEGs_12M_H3_L1 log2FC","Xie Up_DEGs_12M_H3_L1 FDR p")

xie_DEGs <- merge(Sheet_1_TS3E_Up_DEGs_2M_H1_H2_Fig_7G, Sheet_1_TS3F_Up_DEGs_12M_H1_H2_Fig_7G,
                  by.x="GeneID",
                  by.y="GeneID",
                  all.x=TRUE,
                  all.y=TRUE)
xie_DEGs <- merge(xie_DEGs, Sheet_1_TS3H_Up_DEGs_2M_H3_L1_Fig_7G,
                  by.x="GeneID",
                  by.y="GeneID",
                  all.x=TRUE,
                  all.y=TRUE)
xie_DEGs <- merge(xie_DEGs, Sheet_1_TS3I_Up_DEGs_12M_H3_L1_Fig_7G,
                  by.x="GeneID",
                  by.y="GeneID",
                  all.x=TRUE,
                  all.y=TRUE)
colnames(xie_DEGs)[1] <- "external_gene_name"
```
Merge files by gene name
```{r}
all_ncRNA_ageing <- merge(borrett_supp6,xie_DEGs,
                          by.x="external_gene_name",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=TRUE)
```
```{r}
xie_increasingwithage <- xie_DEGs[,c(1,5,9)]
xie_decreasingwithage <- xie_DEGs[,c(1,3,7)]
all_increasing <- merge(xie_increasingwithage,borrett_UP,
                        by.x="external_gene_name",
                        by.y="external_gene_name",
                        all.x=TRUE,
                        all.y=TRUE)
all_decreasing <- merge(xie_decreasingwithage,borrett_down,
                        by.x="external_gene_name",
                        by.y="external_gene_name",
                        all.x=TRUE,
                        all.y=TRUE)
```
Calculate CO - after filtering so only all p values are either NA (for p>0.05) or p<0.05
NAs converted to 0
any p<0.05 converted to 1 and then summed per gene
```{r}
rank_calc <- all_increasing
names <- rank_calc[,1]
rank_calc <- rank_calc[,-c(1)]
rank_calc[!is.na(rank_calc)] <- 1
rank_calc[is.na(rank_calc)] <- 0
rank_calc$`Borett Comparison` <- as.numeric(rank_calc$`Borett Comparison`)

ranks <- cbind(names, rank_calc)
sapply(ranks, class)
chars <- sapply(rank_calc, is.character)
rank_calc[ , chars] <- as.data.frame(apply(rank_calc[ , chars], 2, as.numeric))
increasing_summed_ranks <- cbind(names,
                      rowSums(rank_calc),
                      rank_calc)
colnames(increasing_summed_ranks)[1] <- "external_gene_name"


rank_calc <- all_decreasing
names <- rank_calc[,1]
rank_calc <- rank_calc[,-c(1)]
rank_calc[!is.na(rank_calc)] <- 1
rank_calc[is.na(rank_calc)] <- 0
rank_calc$`Borett Comparison` <- as.numeric(rank_calc$`Borett Comparison`)

ranks <- cbind(names, rank_calc)
sapply(ranks, class)
chars <- sapply(rank_calc, is.character)
rank_calc[ , chars] <- as.data.frame(apply(rank_calc[ , chars], 2, as.numeric))
decreasing_summed_ranks <- cbind(names,
                      rowSums(rank_calc),
                      rank_calc)
colnames(decreasing_summed_ranks)[1] <- "external_gene_name"
```

```{r}

increasing_summed_ranks <- merge(increasing_summed_ranks, ensEMBL2id,
                        by.x="external_gene_name",
                        by.y="external_gene_name",
                        all.x=TRUE,
                        all.y=FALSE)
decreasing_summed_ranks <- merge(decreasing_summed_ranks, ensEMBL2id,
                        by.x="external_gene_name",
                        by.y="external_gene_name",
                        all.x=TRUE,
                        all.y=FALSE)
increasing_summed_ranks$gene_biotype <- as.factor(increasing_summed_ranks$gene_biotype)
decreasing_summed_ranks$gene_biotype <- as.factor(decreasing_summed_ranks$gene_biotype)
increasing_summed_ranks <- increasing_summed_ranks %>% dplyr::filter(increasing_summed_ranks$gene_biotype=="lncRNA"
|gene_biotype=="processed_pseudogene"
|gene_biotype=="TEC"
|gene_biotype=="transcribed_processed_pseudogene"
|gene_biotype=="transcribed_unprocessed_pseudogene"
|gene_biotype=="unprocessed_pseudogene"
|is.na(increasing_summed_ranks$gene_biotype))
decreasing_summed_ranks <- decreasing_summed_ranks %>% dplyr::filter(decreasing_summed_ranks$gene_biotype=="lncRNA"
|gene_biotype=="processed_pseudogene"
|gene_biotype=="TEC"
|gene_biotype=="transcribed_processed_pseudogene"
|gene_biotype=="transcribed_unprocessed_pseudogene"
|gene_biotype=="unprocessed_pseudogene"
|is.na(decreasing_summed_ranks$gene_biotype))
```
```{r}
write.csv(increasing_summed_ranks,"OUTPUT_increasing_with_age.csv")
write.csv(decreasing_summed_ranks,"OUTPUT_decreasing_with_age.csv")

```
## Citations
- Borrett, M.J., Innes, B.T., Jeong, D., Tahmasian, N., Storer, M.A., Bader, G.D., Kaplan, D.R., and Miller, F.D. (2020). Single-Cell Profiling Shows Murine Forebrain Neural Stem Cells Reacquire a Developmental State when Activated for Adult Neurogenesis. Cell Reports 32, 108022. https://doi.org/10.1016/j.celrep.2020.108022.
- Howe, K.L., Achuthan, P., Allen, J., Allen, J., Alvarez-Jarreta, J., Amode, M.R., Armean, I.M., Azov, A.G., Bennett, R., Bhai, J., et al. (2021). Ensembl 2021. Nucleic Acids Research 49, D884???D891. https://doi.org/10.1093/nar/gkaa942.
- Xie, X.P., Laks, D.R., Sun, D., Poran, A., Laughney, A.M., Wang, Z., Sam, J., Belenguer, G., Fari??as, I., Elemento, O., et al. (2020). High-resolution mouse subventricular zone stem-cell niche transcriptome reveals features of lineage, anatomy, and aging. Proceedings of the National Academy of Sciences 117, 31448???31458. https://doi.org/10.1073/pnas.2014389117.
