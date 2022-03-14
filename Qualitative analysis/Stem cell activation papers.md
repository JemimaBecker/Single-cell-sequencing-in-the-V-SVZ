---
title: "Stem cell activation"
author: "Jemima Becker"
date: "23/02/2022"
output: html_document
---

## Papers
Luo, Y., Coskun, V., Liang, A., Yu, J., Cheng, L., Ge, W., Shi, Z., Zhang, K., Li, C., Cui, Y., et al. (2015). Single-Cell Transcriptome Analyses Reveal Signals to Activate Dormant Neural Stem Cells. Cell 161, 1175–1186.
Donega, V., Geest, A. van der, Sluijs, J., Dijk, R.V., Wang, C.C., Basak, O., Pasterkamp, R.J., and Hol, E. (2022). Single-cell profiling of human subventricular zone progenitors identifies SFRP1 as a target for stimulating progenitor activation.

- Borrett, M.J., Innes, B.T., Jeong, D., Tahmasian, N., Storer, M.A., Bader, G.D., Kaplan, D.R., and Miller, F.D. (2020). Single-Cell Profiling Shows Murine Forebrain Neural Stem Cells Reacquire a Developmental State when Activated for Adult Neurogenesis. Cell Reports 32, 108022.
- Zywitza, V., Misios, A., Bunatyan, L., Willnow, T.E., and Rajewsky, N. (2018). Single-Cell Transcriptomics Characterizes Cell Types in the Subventricular Zone and Uncovers Molecular Defects Impairing Adult Neurogenesis. Cell Reports 25, 2457-2469.e8.
- Llorens-Bobadilla, E., Zhao, S., Baser, A., Saiz-Castro, G., Zwadlo, K., and Martin-Villalba, A. (2015). Single-Cell Transcriptomics Reveals a Population of Dormant Neural Stem Cells that Become Activated upon Brain Injury. Cell Stem Cell 17, 329–340.
- Mizrak, D., Bayin, N.S., Yuan, J., Liu, Z., Suciu, R., Niphakis, M.J., Ngo, N., Lum, K.M., Cravatt, B.F., Joyner, A.L., et al. (2019). Single-cell profiling and SCOPE-seq reveal the lineage dynamics of adult neurogenesis and NOTUM as a key V-SVZ regulator. BioRxiv 770610.
- Delgado, A.C., Maldonado-Soto, A.R., Silva-Vargas, V., Mizrak, D., von Känel, T., Tan, K.R., Paul, A., Madar, A., Cuervo, H., Kitajewski, J., et al. (2021). Release of stem cells from quiescence reveals gliogenic domains in the adult mouse brain. Science 372, 1205–1209.
```{r}
setwd("~/Library/CloudStorage/OneDrive-Nexus365/1. V-SVZ lncRNAs (main)/1.1 scRNAseq_paper/Stem cell activation and lineage progression")
```

## Output files
LB_ncRNAs
zywitza_NSCs
Mizrak_all
delgado_all
borett_all

## Code
```{r}
library(readr)
library(readxl)
library(dplyr)
```

Borett et al 2020
```{r setup, include=FALSE}
Sheet_1_E14_DE_genes_between_cortex_and_GE_RPs <- read_csv("Borett et al 2020/1-s2.0-S221112472031007X-mmc4/Sheet 1-E14 DE genes between cortex and GE RPs.csv")
sheet_P6_P7_DE_genes_between_cortical_and_GE_dNSCs <- read_csv("Borett et al 2020/1-s2.0-S221112472031007X-mmc4/sheet-P6 P7 DE genes between cortical and GE dNSCs.csv")
colnames(Sheet_1_E14_DE_genes_between_cortex_and_GE_RPs)[c(2,3,4)] <- c("Borett E14 DE cortex vs GE RPs avg diff","Borett E14 DE cortex vs GE RPs FWER","Borett E14 DE cortex vs GE RPs E14+ve cluster")
colnames(sheet_P6_P7_DE_genes_between_cortical_and_GE_dNSCs)[c(2,3,4)] <- c("Borett P6/P7 DE cortical vs GE dNSC avg diff","Borett P6/P7 DE cortical vs GE dNSC FWER","Borett P6/P7 DE cortical vs GE dNSC P6/P7+ve cluster")
borett_all <- merge(Sheet_1_E14_DE_genes_between_cortex_and_GE_RPs,sheet_P6_P7_DE_genes_between_cortical_and_GE_dNSCs,
                    by.x="Gene",
                    by.y="Gene",
                    all.x=TRUE,
                    all.y=TRUE)
colnames(borett_all)[1] <- "external_gene_name"
```

Couturier et al 2020
```{r setup, include=FALSE}
```

Delgado et al 2021
nb that there are many differnet comparisons made in this study - using three files comparising specific stage transitions
```{r setup, include=FALSE}
delgadocomp3_Pβ_CD133_vsPβ_Table_1 <- read_csv("Delgado et al 2021/DE_tables2/Pβ+CD133+vsPβ+-Table 1.csv")
delgadocomp5_Pβ_EGFR_vsPβ_CD133_Table_1 <- read_csv("Delgado et al 2021/DE_tables2/Pβ+EGFR+vsPβ+CD133+-Table 1.csv")
delgadocomp8_Pβ_EGFR_vsPβ_EGFR_Table_1 <- read_csv("Delgado et al 2021/DE_tables2/Pβ+EGFR+vsPβ-EGFR+-Table 1.csv")

delgadocomp3_Pβ_CD133_vsPβ_Table_1 <- delgadocomp3_Pβ_CD133_vsPβ_Table_1[,-c(5)]
colnames(delgadocomp3_Pβ_CD133_vsPβ_Table_1) <- c("external_gene_name","Delgado comp3 logFC","Delgado comp3 PValue","Delgado comp3 FDR")
delgadocomp5_Pβ_EGFR_vsPβ_CD133_Table_1 <- delgadocomp5_Pβ_EGFR_vsPβ_CD133_Table_1[,-c(5)]
colnames(delgadocomp5_Pβ_EGFR_vsPβ_CD133_Table_1) <- c("external_gene_name","Delgado comp5 logFC","Delgado comp5 PValue","Delgado comp5 FDR")
delgadocomp8_Pβ_EGFR_vsPβ_EGFR_Table_1 <- delgadocomp8_Pβ_EGFR_vsPβ_EGFR_Table_1[,-c(5)]
colnames(delgadocomp8_Pβ_EGFR_vsPβ_EGFR_Table_1) <- c("external_gene_name","Delgado comp8 logFC","Delgado comp8 PValue","Delgado comp8 FDR")

library(dplyr)
delgadocomp5_Pβ_EGFR_vsPβ_CD133_Table_1$`Delgado comp5 PValue`<- 
  as.numeric(delgadocomp5_Pβ_EGFR_vsPβ_CD133_Table_1$`Delgado comp5 PValue`)
delgadocomp3_Pβ_CD133_vsPβ_Table_1$`Delgado comp3 PValue`<-
  as.numeric(delgadocomp3_Pβ_CD133_vsPβ_Table_1$`Delgado comp3 PValue`)
delgadocomp8_Pβ_EGFR_vsPβ_EGFR_Table_1$`Delgado comp8 PValue`<- 
  as.numeric(delgadocomp8_Pβ_EGFR_vsPβ_EGFR_Table_1$`Delgado comp8 PValue`)

delgadocomp3_Pβ_CD133_vsPβ_Table_1 <-delgadocomp3_Pβ_CD133_vsPβ_Table_1 %>%
  dplyr::filter(delgadocomp3_Pβ_CD133_vsPβ_Table_1$`Delgado comp3 PValue` < 0.05)
delgadocomp5_Pβ_EGFR_vsPβ_CD133_Table_1 <- delgadocomp5_Pβ_EGFR_vsPβ_CD133_Table_1 %>%   
  dplyr::filter(delgadocomp5_Pβ_EGFR_vsPβ_CD133_Table_1$`Delgado comp5 PValue`< 0.05)

delgadocomp8_Pβ_EGFR_vsPβ_EGFR_Table_1 <- delgadocomp8_Pβ_EGFR_vsPβ_EGFR_Table_1 %>% 
  dplyr::filter(delgadocomp8_Pβ_EGFR_vsPβ_EGFR_Table_1$`Delgado comp8 PValue`< 0.05)

delgado_all <- merge(delgadocomp3_Pβ_CD133_vsPβ_Table_1,delgadocomp5_Pβ_EGFR_vsPβ_CD133_Table_1,
                     by.x="external_gene_name",
                     by.y="external_gene_name",
                     all.x=TRUE,
                     all.y=TRUE)
delgado_all <- merge(delgado_all,delgadocomp8_Pβ_EGFR_vsPβ_EGFR_Table_1,
                     by.x="external_gene_name",
                     by.y="external_gene_name",
                     all.x=TRUE,
                     all.y=TRUE)
```

Mizrak et al 2019
there are three datasets here - processed by authors to find genes important for these transitions
```{r setup, include=FALSE}
Mizrak_SCOPE_aNSC_Table_1 <- read_csv("Mizrak et al 2019/mmc2/SCOPE_aNSC-Table 1.csv")
Mizrak_GCERT2_aNSC_Table_1 <- read_csv("Mizrak et al 2019/mmc2/GCERT2_aNSC-Table 1.csv")
Mizrak_NESFLPO_aNSC_Table_1 <- read_csv("Mizrak et al 2019/mmc2/NESFLPO_aNSC-Table 1.csv")
colnames(Mizrak_SCOPE_aNSC_Table_1)[c(3,4)] <- c("Mizrak SCOPE_aNSC FDR","Mizrak SCOPE_aNSC log2 Effect")
colnames(Mizrak_GCERT2_aNSC_Table_1)[c(3,4)] <- c("Mizrak GCERT2_aNSC FDR","Mizrak GCERT2_aNSC log2 Effect")
colnames(Mizrak_NESFLPO_aNSC_Table_1)[c(3,4)] <- c("Mizrak NESFLPO_aNSC FDR","Mizrak NESFLPO_aNSC log2 Effect")

Mizrak_SCOPE_aNSC_Table_1$`Mizrak SCOPE_aNSC FDR` <- as.numeric(Mizrak_SCOPE_aNSC_Table_1$`Mizrak SCOPE_aNSC FDR`)
Mizrak_GCERT2_aNSC_Table_1$`Mizrak GCERT2_aNSC FDR` <- as.numeric(Mizrak_GCERT2_aNSC_Table_1$`Mizrak GCERT2_aNSC FDR`)
Mizrak_NESFLPO_aNSC_Table_1$`Mizrak NESFLPO_aNSC FDR` <- as.numeric(Mizrak_NESFLPO_aNSC_Table_1$`Mizrak NESFLPO_aNSC FDR`)

Mizrak_SCOPE_aNSC_Table_1 <- Mizrak_SCOPE_aNSC_Table_1 %>% dplyr::filter(Mizrak_SCOPE_aNSC_Table_1$`Mizrak SCOPE_aNSC FDR` < 0.05)
Mizrak_GCERT2_aNSC_Table_1 <- Mizrak_GCERT2_aNSC_Table_1 %>% dplyr::filter(Mizrak_GCERT2_aNSC_Table_1$`Mizrak GCERT2_aNSC FDR` < 0.05)
Mizrak_NESFLPO_aNSC_Table_1 <- Mizrak_NESFLPO_aNSC_Table_1 %>% dplyr::filter(Mizrak_NESFLPO_aNSC_Table_1$`Mizrak NESFLPO_aNSC FDR` < 0.05)

Mizrak_SCOPE_aNSC_Table_1 <- Mizrak_SCOPE_aNSC_Table_1[,c(2:4)]
Mizrak_GCERT2_aNSC_Table_1 <- Mizrak_GCERT2_aNSC_Table_1[,c(2:4)]
Mizrak_NESFLPO_aNSC_Table_1 <- Mizrak_NESFLPO_aNSC_Table_1[,c(2:4)]
Mizrak_all <- merge(Mizrak_SCOPE_aNSC_Table_1, Mizrak_GCERT2_aNSC_Table_1,
                    by.x="gene",
                    by.y="gene",
                    all.x=TRUE,
                    all.y=TRUE)
Mizrak_all <- merge(Mizrak_all, Mizrak_NESFLPO_aNSC_Table_1,
                    by.x="gene",
                    by.y="gene",
                    all.x=TRUE,
                    all.y=TRUE)
colnames(Mizrak_all)[1] <- "external_gene_name"
```

Zywitza et al 2018
```{r setup, include=FALSE}
scl8_aNSCs_Table_1 <- read_csv("Zywitza et al 2018/mmc4 (1)/scl8_aNSCs-Table 1.csv")
scl7_qNSCs_Table_1 <- read_csv("Zywitza et al 2018/mmc4 (1)/scl7_qNSCs-Table 1.csv")
cl17_NSCs_Table_1 <- read_csv("Zywitza et al 2018/mmc2/cl17_NSCs-Table 1.csv")
colnames(scl8_aNSCs_Table_1)[6] <- "cluster or subcluster"
colnames(scl7_qNSCs_Table_1)[6] <- "cluster or subcluster"
colnames(cl17_NSCs_Table_1)[6] <- "cluster or subcluster"
zywitza_NSCs <- rbind(scl8_aNSCs_Table_1,scl7_qNSCs_Table_1,cl17_NSCs_Table_1)
colnames(zywitza_NSCs) <- c("external_gene_name","Zywitza Pval","Zywitza avg_diff","Zywitza Pct1","Zywitza Pct2","Zywitza cluster or subcluster")
```

Llorens-babadilla et al 2015 NB that this activation upon induction of injury
```{r}
library(readxl)
library(dplyr)
LB <- read_excel("Llorens-Bobadilla et al 2015/mmc4 (1).xlsx")
LB$gene_biotype <- as.factor(LB$gene_biotype)
levels(LB$gene_biotype)
LB_ncRNAs <- LB %>% dplyr::filter(LB$gene_biotype == "antisense"|
                                    gene_biotype == "lincRNA"|
                                    gene_biotype == "misc_RNA"|
                                    gene_biotype == "polymorphic_pseudogene"|
                                    gene_biotype == "processed_transcript"|
                                    gene_biotype == "pseudogene")
LB_ncRNAs <- LB_ncRNAs[,c(8,2:7)]
colnames(LB_ncRNAs) <- c("external_gene_name","LB basemean","LB log2FC","LB IFcSE","LB stat","LB pvalue","LB padj")
```
merge and annotate
```{r}
activation_merged <- merge(borett_all, delgado_all,
                           by.x="external_gene_name",
                           by.y="external_gene_name",
                           all.x=TRUE,
                           all.y=TRUE)
activation_merged <- merge(activation_merged, Mizrak_all,
                           by.x="external_gene_name",
                           by.y="external_gene_name",
                           all.x=TRUE,
                           all.y=TRUE)
activation_merged <- merge(activation_merged, zywitza_NSCs,
                           by.x="external_gene_name",
                           by.y="external_gene_name",
                           all.x=TRUE,
                           all.y=TRUE)
activation_merged <- merge(activation_merged, LB_ncRNAs,
                           by.x="external_gene_name",
                           by.y="external_gene_name",
                           all.x=TRUE,
                           all.y=TRUE)

activation_merged_ann <- merge(activation_merged, ensEMBL2id,
                               by.x="external_gene_name",
                               by.y="external_gene_name",
                               all.x=TRUE,
                               all.y=FALSE)
activation_merged_ann$gene_biotype <- as.factor(activation_merged_ann$gene_biotype)
levels(activation_merged_ann$gene_biotype)
library(dplyr)
activation_lncRNAs <- activation_merged_ann %>% dplyr::filter(activation_merged_ann$gene_biotype == "lncRNA"
| gene_biotype == "misc_RNA"                          
| gene_biotype == "polymorphic_pseudogene"             
| gene_biotype == "processed_pseudogene"
| gene_biotype == "pseudogene"                        
| gene_biotype == "ribozyme"                             
| gene_biotype == "TEC"
| gene_biotype == "transcribed_processed_pseudogene"  
| gene_biotype == "transcribed_unitary_pseudogene"     
| gene_biotype == "transcribed_unprocessed_pseudogene"
| gene_biotype == "unitary_pseudogene"                
| gene_biotype == "unprocessed_pseudogene"    
| is.na(activation_merged_ann$gene_biotype))
```
```{r}
activation_lncRNAs_pvalues <-as.data.frame(cbind(activation_lncRNAs$external_gene_name,                                           
activation_lncRNAs$`Borett E14 DE cortex vs GE RPs FWER`,                        
activation_lncRNAs$`Borett P6/P7 DE cortical vs GE dNSC FWER`,            
activation_lncRNAs$`Delgado comp3 PValue`,                                 
activation_lncRNAs$`Delgado comp5 PValue`,                                
activation_lncRNAs$`Delgado comp8 PValue`,                                  
activation_lncRNAs$`Mizrak SCOPE_aNSC FDR`,                                
activation_lncRNAs$`Mizrak GCERT2_aNSC FDR`,                               
activation_lncRNAs$`Mizrak NESFLPO_aNSC FDR`,                              
activation_lncRNAs$`Zywitza Pval`,                                        
activation_lncRNAs$`LB padj`))

colnames(activation_lncRNAs_pvalues) <- c(
"external_gene_name",
"Borett E14 DE cortex vs GE RPs FWER",                        
"Borett P6/P7 DE cortical vs GE dNSC FWER",            
"Delgado comp3 PValue",                                 
"Delgado comp5 PValue",                                
"Delgado comp8 PValue",                                  
"Mizrak SCOPE_aNSC FDR",                                
"Mizrak GCERT2_aNSC FDR",                               
"Mizrak NESFLPO_aNSC FDR",                              
"Zywitza Pval",                                        
"LB padj")
```

Calculate mean p-value across columns
```{r}
activation_lncRNAs_pvalues$`Borett E14 DE cortex vs GE RPs FWER` <- as.numeric(activation_lncRNAs_pvalues$`Borett E14 DE cortex vs GE RPs FWER`)                    
activation_lncRNAs_pvalues$`Borett P6/P7 DE cortical vs GE dNSC FWER`<- as.numeric(activation_lncRNAs$`Borett P6/P7 DE cortical vs GE dNSC FWER`)      
activation_lncRNAs_pvalues$`Delgado comp3 PValue`<- as.numeric(activation_lncRNAs_pvalues$`Delgado comp3 PValue`)
activation_lncRNAs_pvalues$`Delgado comp5 PValue`<- as.numeric(activation_lncRNAs_pvalues$`Delgado comp5 PValue`)
activation_lncRNAs_pvalues$`Delgado comp8 PValue`<- as.numeric(activation_lncRNAs_pvalues$`Delgado comp8 PValue`)
activation_lncRNAs_pvalues$`Mizrak SCOPE_aNSC FDR`<- as.numeric(activation_lncRNAs_pvalues$`Mizrak SCOPE_aNSC FDR`)
activation_lncRNAs_pvalues$`Mizrak GCERT2_aNSC FDR`<- as.numeric(activation_lncRNAs_pvalues$`Mizrak GCERT2_aNSC FDR`)
activation_lncRNAs_pvalues$`Mizrak NESFLPO_aNSC FDR`<- as.numeric(activation_lncRNAs_pvalues$`Mizrak NESFLPO_aNSC FDR`)
activation_lncRNAs_pvalues$`Zywitza Pval`<- as.numeric(activation_lncRNAs_pvalues$`Zywitza Pval`)
activation_lncRNAs_pvalues$`LB padj` <- as.numeric(activation_lncRNAs_pvalues$`LB padj`)

activation_lncRNAs_pvalues2 <- activation_lncRNAs_pvalues[,-c(1)]
activation_lncRNAs_pvalues_mean <- rowMeans(activation_lncRNAs_pvalues2, na.rm=TRUE)
activation_lncRNAs_means <- cbind(activation_lncRNAs_pvalues$external_gene_name,activation_lncRNAs_pvalues_mean ,activation_lncRNAs_pvalues2)
```

calculate summed ranks across columns
```{r}
rank_calc <- activation_lncRNAs_means
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

activation_lncRNAs_means <- cbind(activation_lncRNAs_means[,1],summed_ranks$`rowSums(rank_calc)`,activation_lncRNAs_means[,-c(1)])
colnames(activation_lncRNAs_means)[c(1,2,3)] <- c("external_gene_name","summed_rank","mean_pvalue")
activation_lncRNAs_means$summed_rank <- as.numeric(activation_lncRNAs_means$summed_rank)
activation_lncRNAs_means$mean_pvalue <- as.numeric(activation_lncRNAs_means$mean_pvalue)
log_pvalues <- -log10(activation_lncRNAs_means$mean_pvalue)
multiplied <- activation_lncRNAs_means$summed_rank*log_pvalues
activation_lncRNA_ranks <- cbind(activation_lncRNAs_means$external_gene_name,multiplied,log_pvalues,activation_lncRNAs_means[,-c(1)])
colnames(activation_lncRNA_ranks)[2:5] <- c("Activation multiplied","Activation -log10(mean pvalue)","Activation summed ranks","Activation mean pvalues")
write.csv(activation_lncRNA_ranks,"activation_lncRNA_ranks_24-02-22.csv")
write.csv(activation_lncRNAs,"activation_lncRNAs_24-02-22.csv")
```
