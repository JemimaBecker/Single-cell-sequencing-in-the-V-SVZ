---
title: "Differential expression with ageing"
author: "Jemima Becker"
date: "22/02/2022"
output: word_document
---
## Papers
Kalamakis, G., Brüne, D., Ravichandran, S., Bolz, J., Fan, W., Ziebell, F., Stiehl, T., Catalá-Martinez, F., Kupke, J., Zhao, S., et al. (2019). Quiescence Modulates Stem Cell Maintenance and Regenerative Capacity in the Aging Brain. Cell 176, 1407-1419.e14.
Xie, X.P., Laks, D.R., Sun, D., Poran, A., Laughney, A.M., Wang, Z., Sam, J., Belenguer, G., Fariñas, I., Elemento, O., et al. (2020). High-resolution mouse subventricular zone stem-cell niche transcriptome reveals features of lineage, anatomy, and aging. PNAS 117, 31448–31458.

## Code

output files: 
xie_DEGs
kalamakis_ncRNAs
all_ncRNAs_ageing2

### load data and required packages

Xie et al 2020
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
```

Kalamakis et al 2018
```{r}
X23_vs_2_MO_NSCs_Smartseq2_DESe_Table_1 <- read_csv("Kalamakis et al 2018/Kalamakis et al 2018/23 vs 2 MO NSCs Smartseq2 (DESe-Table 1.csv")
kalamakis_sig <- X23_vs_2_MO_NSCs_Smartseq2_DESe_Table_1 %>% dplyr::filter(X23_vs_2_MO_NSCs_Smartseq2_DESe_Table_1$pvalue < 0.05)
kalamakis_sig$gene_biotype <- as.factor(kalamakis_sig$gene_biotype)
levels(kalamakis_sig$gene_biotype)
kalamakis_ncRNAs <- kalamakis_sig %>% dplyr::filter(kalamakis_sig$gene_biotype ==
"antisense"                       
| gene_biotype== "bidirectional_promoter_lncRNA"      
| gene_biotype== "lincRNA"                           
| gene_biotype== "misc_RNA"                           
| gene_biotype== "polymorphic_pseudogene"            
| gene_biotype== "processed_pseudogene"               
| gene_biotype== "processed_transcript"                  
| gene_biotype== "pseudogene"                         
| gene_biotype== "sense_intronic"                     
| gene_biotype== "TEC"                               
| gene_biotype== "transcribed_processed_pseudogene"   
| gene_biotype== "transcribed_unitary_pseudogene"     
| gene_biotype== "transcribed_unprocessed_pseudogene"
| gene_biotype== "unprocessed_pseudogene"   
| is.na(kalamakis_sig$gene_biotype))
colnames(kalamakis_ncRNAs)[c(3,6)] <- c("Kalamakis LogFC","Kalamakis PValue")
```

add them together 
```{r}
all_ncRNA_ageing <- merge(kalamakis_ncRNAs,xie_DEGs,
                          by.x="gene_symbol",
                          by.y="GeneID",
                          all.x=TRUE,
                          all.y=TRUE)
all_ncRNA_ageing2 <- all_ncRNA_ageing %>% dplyr::filter(!is.na(all_ncRNA_ageing$gene_symbol)) #removes the terms that could not be annotated because they have depreciated ensembl IDs and are no longer in the database

ageing_pvalues <- as.data.frame(cbind(all_ncRNA_ageing2$gene_symbol,
all_ncRNA_ageing2$`Kalamakis PValue`,
all_ncRNA_ageing2$`Xie Up_DEGs_2M_H1_H2 FDR p` ,
all_ncRNA_ageing2$`Xie Up_DEGs_12M_H1_H2 FDR p`,
all_ncRNA_ageing2$`Xie Up_DEGs_2M_H3_L1 FDR p` ,
all_ncRNA_ageing2$`Xie Up_DEGs_12M_H3_L1 FDR p` ))


colnames(ageing_pvalues) <- c("external_gene_name",
                              "Kalamakis PValue",
                              "Xie Up_DEGs_2M_H1_H2 FDR p",
                              "Xie Up_DEGs_12M_H1_H2 FDR p",
                              "Xie Up_DEGs_2M_H3_L1 FDR p",
                              "Xie Up_DEGs_12M_H3_L1 FDR p")

ageing_lncRNAs_pvalues2 <- ageing_pvalues[,-c(1)]

ageing_lncRNAs_pvalues2$`Kalamakis PValue`<- as.numeric(ageing_lncRNAs_pvalues2$`Kalamakis PValue`)
ageing_lncRNAs_pvalues2$`Xie Up_DEGs_2M_H1_H2 FDR p` <- as.numeric(ageing_lncRNAs_pvalues2$`Xie Up_DEGs_2M_H1_H2 FDR p`)
ageing_lncRNAs_pvalues2$`Xie Up_DEGs_12M_H1_H2 FDR p`<- as.numeric(ageing_lncRNAs_pvalues2$`Xie Up_DEGs_12M_H1_H2 FDR p`)
ageing_lncRNAs_pvalues2$`Xie Up_DEGs_2M_H3_L1 FDR p` <- as.numeric(ageing_lncRNAs_pvalues2$`Xie Up_DEGs_2M_H3_L1 FDR p`)
ageing_lncRNAs_pvalues2$`Xie Up_DEGs_12M_H3_L1 FDR p` <- as.numeric(ageing_lncRNAs_pvalues2$`Xie Up_DEGs_12M_H3_L1 FDR p`)

ageing_lncRNAs_pvalues_mean <- rowMeans(ageing_lncRNAs_pvalues2, na.rm=TRUE)
ageing_lncRNAs_means <- cbind(ageing_pvalues$external_gene_name,ageing_lncRNAs_pvalues_mean ,ageing_lncRNAs_pvalues2)  

rank_calc <- ageing_pvalues
names <- rank_calc[,1]
rank_calc <- rank_calc[,-c(1)]
rank_calc[!is.na(rank_calc)] <- 1
rank_calc[is.na(rank_calc)] <- 0
ranks <- cbind(names, rank_calc)
sapply(ranks, class)
chars <- sapply(rank_calc, is.character)
rank_calc[ , chars] <- as.data.frame(apply(rank_calc[ , chars], 2, as.numeric))
summed_ranks <- cbind(names,
                      rowSums(rank_calc),
                      rank_calc)
```
```{r}
colnames(ageing_lncRNAs_means)[1] <- "external_gene_name"
ageing_lncRNAs_means <- cbind(ageing_lncRNAs_means$external_gene_name,summed_ranks$`rowSums(rank_calc)`,ageing_lncRNAs_means[,-c(1)])
```

```{r}
colnames(ageing_lncRNAs_means)[c(1,2,3)] <- c("external_gene_name","summed_rank","mean_pvalue")
ageing_lncRNAs_means$summed_rank <- as.numeric(ageing_lncRNAs_means$summed_rank)
ageing_lncRNAs_means$mean_pvalue <- as.numeric(ageing_lncRNAs_means$mean_pvalue)
log_pvalues <- -log10(ageing_lncRNAs_means$mean_pvalue)
multiplied <- ageing_lncRNAs_means$summed_rank*log_pvalues
ageing_lncRNA_ranks <- cbind(ageing_lncRNAs_means$external_gene_name,multiplied,log_pvalues,ageing_lncRNAs_means[,-c(1)])
colnames(ageing_lncRNA_ranks)[1:5] <- c("external_gene_name","Ageing multiplied","Ageing -log10(mean pvalue)","Ageing summed ranks","Ageing mean pvalues")
write.csv(ageing_lncRNA_ranks,"Ageing_lncRNA_ranks_24-02-22.csv")
write.csv(all_ncRNA_ageing2,"Ageing_lncRNAs_24-02-22.csv")


```
