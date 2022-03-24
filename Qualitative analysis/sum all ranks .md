---
title: "All together"
author: "Jemima Becker"
date: "24/03/2022"
output: html_document
---

```{r setup, include=FALSE}
names_sums_celltype <- celltype_summed_ranks[,c(1,11)]
names_sums_increasing <- increasing_summed_ranks[,c(1,2)]
names_sums_decreasing <- decreasing_summed_ranks[,c(1,2)]
names_sums_activation <- activation_summed_ranks[,c(1,8)]

colnames(names_sums_celltype)[2] <- "Rank sum cell type"
colnames(names_sums_increasing)[2] <- "Rank sum increasing"
colnames(names_sums_decreasing)[2] <- "Rank sum decreasing"
colnames(names_sums_activation)[2] <- "Rank sum activation"

all_sums <- merge(names_sums_celltype,names_sums_increasing,
                  by.x="external_gene_name",
                  by.y="external_gene_name",
                  all.x=TRUE,
                  all.y=TRUE)
all_sums <- merge(all_sums,names_sums_decreasing,
                  by.x="external_gene_name",
                  by.y="external_gene_name",
                  all.x=TRUE,
                  all.y=TRUE)
all_sums <- merge(all_sums,names_sums_activation,
                  by.x="external_gene_name",
                  by.y="external_gene_name",
                  all.x=TRUE,
                  all.y=TRUE)
```

```{r}
names <- all_sums[,1]
ranks <- all_sums[,-c(1)]
ranks[is.na(ranks)] <- 0
all_summed_ranks <- cbind(names,
                      rowSums(ranks),
                      ranks)
colnames(all_summed_ranks)[2] <- "cumulative sums"
all_summed_ranks <- merge(all_summed_ranks, ensEMBL2id,
                          by.x="names",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=FALSE)
write.csv(all_summed_ranks,"OUTPUT_all_summed_ranks.csv")
```

```{r}
info <- Sheet_1_healthy_rank_info_8_11_21[,c(2,6:13)]
all_info <- merge(all_summed_ranks, info,
by.x="names",
by.y="Gene",
all.x=TRUE,
all.y=FALSE)
write.csv(all_info,"OUTPUT_allinfo.csv")
```
