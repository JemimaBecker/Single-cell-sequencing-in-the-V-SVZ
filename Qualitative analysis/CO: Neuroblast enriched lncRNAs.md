# Calculation of cumulative occurences for neuroblast-enriched lncRNAs

Input data is results from FindAllMarkers listing lncRNAs that are significantly enriched in neuroblast clusters relative to other cell types.


## Code 

1: Make a single dataframe with the neuroblast enriched genes from each study
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
```
2: Calculate cumulative occurences
```{r}
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

## Data references

del Águila, Á., Adam, M., Ullom, K., Shaw, N., Qin, S., Ehrman, J., Nardini, D., Salomone, J., Gebelein, B., Lu, Q.R., et al. (2022). Olig2 defines a subset of neural stem cells that produce specific olfactory bulb interneuron subtypes in the subventricular zone of adult mice. Development 149, dev200028. https://doi.org/10.1242/dev.200028.
Kriska, J., Janeckova, L., Kirdajova, D., Honsa, P., Knotek, T., Dzamba, D., Kolenicova, D., Butenko, O., Vojtechova, M., Capek, M., et al. (2021). Wnt/β-Catenin Signaling Promotes Differentiation of Ischemia-Activated Adult Neural Stem/Progenitor Cells to Neuronal Precursors. Frontiers in Neuroscience 15. .
Mizrak, D., Bayin, N.S., Yuan, J., Liu, Z., Suciu, R., Niphakis, M.J., Ngo, N., Lum, K.M., Cravatt, B.F., Joyner, A.L., et al. (2019a). Single-cell profiling and SCOPE-seq reveal the lineage dynamics of adult neurogenesis and NOTUM as a key V-SVZ regulator. BioRxiv 770610. https://doi.org/10.1101/770610.
Mizrak, D., Levitin, H.M., Delgado, A.C., Crotet, V., Yuan, J., Chaker, Z., Silva-Vargas, V., Sims, P.A., and Doetsch, F. (2019b). Single-Cell Analysis of Regional Differences in Adult V-SVZ Neural Stem Cell Lineages. Cell Reports 26, 394-406.e5. https://doi.org/10.1016/j.celrep.2018.12.044.
Zywitza, V., Misios, A., Bunatyan, L., Willnow, T.E., and Rajewsky, N. (2018). Single-Cell Transcriptomics Characterizes Cell Types in the Subventricular Zone and Uncovers Molecular Defects Impairing Adult Neurogenesis. Cell Reports 25, 2457-2469.e8. https://doi.org/10.1016/j.celrep.2018.11.003.
