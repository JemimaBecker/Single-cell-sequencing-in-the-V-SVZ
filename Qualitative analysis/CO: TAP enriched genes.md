```{r}
TAP_lncRNAS <- all_lncRNAs_merged[,c(
"external_gene_name",
"Del Aguila 2022 aNSCs TAPs"      , 
"Cebrian-Silla 2021 A-cell cluster",
"Kriska 2021 TAP",
"Mizrak 2019 TAP",
"Mizrak 2020 NESFPLO TAC",
"Mizrak 2020 GCERT2 TAC",
"Mizrak 2020 SCOPE TAC",
"Zamboni 2020 TAP",
"Zywitza 2020 TAP"    )]    

TAP_lncRNAS$`Zamboni 2020 TAP` <- as.character(TAP_lncRNAS$`Zamboni 2020 TAP`)


names6 <- TAP_lncRNAS[,1]
ranks6 <- TAP_lncRNAS[,-c(1)]
ranks6[!is.na(ranks6)] <- as.numeric(1)
ranks6[is.na(ranks6)] <- as.numeric(0)


ranks7 <- data.frame(apply(ranks6, 2, function(x) as.numeric(as.character(x))))

sapply(ranks7, class)

TAP_all_summed_ranks <- cbind(names,
                      rowSums(ranks7),
                      ranks7)

colnames(TAP_all_summed_ranks)[2] <- "cumulative sums TAP"
TAP_all_summed_ranks <- merge(TAP_all_summed_ranks, ensEMBL2id,
                          by.x="names",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=FALSE)
write.csv(TAP_all_summed_ranks,"TAP_all_summed_ranks.csv")
```
