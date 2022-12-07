# Calculation of cumulative occurences for transit amplifying progenitor-enriched lncRNAs

Input data is results from FindAllMarkers listing lncRNAs that are significantly enriched in TAP/TAC clusters relative to other cell types.

## Code

1: Load data
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
```
2: Caluclate cumulative occurences
```{r}
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

## Data references

- del Águila, Á., Adam, M., Ullom, K., Shaw, N., Qin, S., Ehrman, J., Nardini, D., Salomone, J., Gebelein, B., Lu, Q.R., et al. (2022). Olig2 defines a subset of neural stem cells that produce specific olfactory bulb interneuron subtypes in the subventricular zone of adult mice. Development 149, dev200028. https://doi.org/10.1242/dev.200028.
- Cebrian-Silla, A., Nascimento, M.A., Redmond, S.A., Mansky, B., Wu, D., Obernier, K., Romero Rodriguez, R., Gonzalez-Granero, S., García-Verdugo, J.M., Lim, D.A., et al. (2021). Single-cell analysis of the ventricular-subventricular zone reveals signatures of dorsal and ventral adult neurogenesis. ELife 10, e67436. https://doi.org/10.7554/eLife.67436.
- Kriska, J., Janeckova, L., Kirdajova, D., Honsa, P., Knotek, T., Dzamba, D., Kolenicova, D., Butenko, O., Vojtechova, M., Capek, M., et al. (2021). Wnt/β-Catenin Signaling Promotes Differentiation of Ischemia-Activated Adult Neural Stem/Progenitor Cells to Neuronal Precursors. Frontiers in Neuroscience 15. .
- Mizrak, D., Bayin, N.S., Yuan, J., Liu, Z., Suciu, R., Niphakis, M.J., Ngo, N., Lum, K.M., Cravatt, B.F., Joyner, A.L., et al. (2019a). Single-cell profiling and SCOPE-seq reveal the lineage dynamics of adult neurogenesis and NOTUM as a key V-SVZ regulator. BioRxiv 770610. https://doi.org/10.1101/770610.
- Mizrak, D., Levitin, H.M., Delgado, A.C., Crotet, V., Yuan, J., Chaker, Z., Silva-Vargas, V., Sims, P.A., and Doetsch, F. (2019b). Single-Cell Analysis of Regional Differences in Adult V-SVZ Neural Stem Cell Lineages. Cell Reports 26, 394-406.e5. https://doi.org/10.1016/j.celrep.2018.12.044.
- Zamboni, M., Llorens-Bobadilla, E., Magnusson, J.P., and Frisén, J. (2020). A Widespread Neurogenic Potential of Neocortical Astrocytes Is Induced by Injury. Cell Stem Cell 27, 605-617.e5. https://doi.org/10.1016/j.stem.2020.07.006.
- Zywitza, V., Misios, A., Bunatyan, L., Willnow, T.E., and Rajewsky, N. (2018). Single-Cell Transcriptomics Characterizes Cell Types in the Subventricular Zone and Uncovers Molecular Defects Impairing Adult Neurogenesis. Cell Reports 25, 2457-2469.e8. https://doi.org/10.1016/j.celrep.2018.11.003.
