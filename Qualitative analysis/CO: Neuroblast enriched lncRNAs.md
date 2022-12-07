```{r}
NB_lncRNAS <- all_lncRNAs_merged[,c(
"external_gene_name",
"Del Aguila 2022 NB",
"Kriska 2021 NB1",
"Kriska 2021 NB2" ,                              
"Mizrak 2019 NB",
"Mizrak 2020 NESFPLO NB",
"Mizrak 2020 GCERT2 NB",
"Mizrak 2020 SCOPE aNSC",
"Mizrak 2020 SCOPE NB",
"Zamboni 2020 NB",
"Zywitza 2020 NB"  )]

NB_lncRNAS$`Zamboni 2020 NB` <- as.character(NB_lncRNAS$`Zamboni 2020 NB`)

names4 <- NB_lncRNAS[,1]
ranks4 <- NB_lncRNAS[,-c(1)]
ranks4[!is.na(ranks4)] <- as.numeric(1)
ranks4[is.na(ranks4)] <- as.numeric(0)


ranks5 <- data.frame(apply(ranks4, 2, function(x) as.numeric(as.character(x))))

sapply(ranks5, class)

NB_all_summed_ranks <- cbind(names,
                      rowSums(ranks5),
                      ranks5)

colnames(NB_all_summed_ranks)[2] <- "cumulative sums NB"
NB_all_summed_ranks <- merge(NB_all_summed_ranks, ensEMBL2id,
                          by.x="names",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=FALSE)
write.csv(NB_all_summed_ranks,"NB_all_summed_ranks.csv")
```
