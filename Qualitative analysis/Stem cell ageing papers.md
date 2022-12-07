```{r}
up_in_quiescent_with_age <- aging_lncRNAs[,c("external_gene_name","Xie 2020 up 12M qNSCs", "Kalamakis_qNSC1_sig_pos","Kalamakis_qNSC2_sig_pos")]
down_in_active_with_age <- aging_lncRNAs[,c("external_gene_name","Kalamakis_aNSC1_sig_neg","Kalamakis_aNSC2_sig_neg","Xie 2020 up 2M aNSCs" )]
down_in_quiescent_with_age <- aging_lncRNAs[,c("external_gene_name","Kalamakis_qNSC1_sig_neg","Kalamakis_qNSC2_sig_neg" , "Xie 2020 up 2M qNSCs")]
```

```{r}
#up_in_quiescent_with_age

names <- up_in_quiescent_with_age[,1]
ranks <- up_in_quiescent_with_age[,-c(1)]
ranks[!is.na(ranks)] <- as.numeric(1)
ranks[is.na(ranks)] <- as.numeric(0)


ranks2 <- data.frame(apply(ranks, 2, function(x) as.numeric(as.character(x))))

sapply(ranks2, class)

up_in_quiescent_with_age_ranksum <- cbind(names,
                      rowSums(ranks2),
                      ranks2)

colnames(up_in_quiescent_with_age_ranksum)[2] <- "cumulative sums up_in_quiescent_with_age"
up_in_quiescent_with_age_ranksum <- merge(up_in_quiescent_with_age_ranksum, ensEMBL2id,
                          by.x="names",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=FALSE)
write.csv(up_in_quiescent_with_age_ranksum,"up_in_quiescent_with_age_ranksum.csv")

```

```{r}
#down_in_active_with_age

names <- down_in_active_with_age[,1]
ranks <- down_in_active_with_age[,-c(1)]
ranks[!is.na(ranks)] <- as.numeric(1)
ranks[is.na(ranks)] <- as.numeric(0)


ranks2 <- data.frame(apply(ranks, 2, function(x) as.numeric(as.character(x))))

sapply(ranks2, class)

down_in_active_with_age_ranksum <- cbind(names,
                      rowSums(ranks2),
                      ranks2)

colnames(down_in_active_with_age_ranksum)[2] <- "cumulative sums down_in_active_with_age"
down_in_active_with_age_ranksum <- merge(down_in_active_with_age_ranksum, ensEMBL2id,
                          by.x="names",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=FALSE)
write.csv(down_in_active_with_age_ranksum,"down_in_active_with_age_ranksum.csv")

```

```{r}
#down_in_quiescent_with_age

names <- down_in_quiescent_with_age[,1]
ranks <- down_in_quiescent_with_age[,-c(1)]
ranks[!is.na(ranks)] <- as.numeric(1)
ranks[is.na(ranks)] <- as.numeric(0)


ranks2 <- data.frame(apply(ranks, 2, function(x) as.numeric(as.character(x))))

sapply(ranks2, class)

down_in_quiescent_with_age_ranksum <- cbind(names,
                      rowSums(ranks2),
                      ranks2)

colnames(down_in_quiescent_with_age_ranksum)[2] <- "cumulative sums down_in_quiescent_with_age"
down_in_quiescent_with_age_ranksum <- merge(down_in_quiescent_with_age_ranksum, ensEMBL2id,
                          by.x="names",
                          by.y="external_gene_name",
                          all.x=TRUE,
                          all.y=FALSE)
write.csv(down_in_quiescent_with_age_ranksum,"down_in_quiescent_with_age_ranksum.csv")

```
