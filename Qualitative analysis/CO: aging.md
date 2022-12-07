# Calculation of cumulative occurences for lncRNAs that accumulate with aging


```{r}
up_in_quiescent_with_age <- aging_lncRNAs[,c("external_gene_name","Xie 2020 up 12M qNSCs", "Kalamakis_qNSC1_sig_pos","Kalamakis_qNSC2_sig_pos")]

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
## Data references
- Kalamakis, G., Brüne, D., Ravichandran, S., Bolz, J., Fan, W., Ziebell, F., Stiehl, T., Catalá-Martinez, F., Kupke, J., Zhao, S., et al. (2019). Quiescence Modulates Stem Cell Maintenance and Regenerative Capacity in the Aging Brain. Cell 176, 1407-1419.e14. https://doi.org/10.1016/j.cell.2019.01.040.
- Xie, X.P., Laks, D.R., Sun, D., Poran, A., Laughney, A.M., Wang, Z., Sam, J., Belenguer, G., Fariñas, I., Elemento, O., et al. (2020). High-resolution mouse subventricular zone stem-cell niche transcriptome reveals features of lineage, anatomy, and aging. PNAS 117, 31448–31458. https://doi.org/10.1073/pnas.2014389117.
