# Data summary for scRNAseq meta-analysis 

<img width="639" alt="Screenshot 2022-04-19 at 13 57 50" src="https://user-images.githubusercontent.com/67189202/164009216-30ee741b-264d-4d9c-8446-53f16fcbf2ef.png">

> Figure 3: Data summary for scRNAseq meta-analysis. 

> A: Total number of studies in which each lncRNA is differentially expressed across all categories, data is expressed as log2(number of lncRNAs) vs the number of studies that each gene was differentially expressed in. 

> B: Overlap between lncRNAs differentially expressed with stem cell activation (n=750) and those that are downregulated with ageing (n=28). 

> C: EnrichR gene ontology results for transcription factors that bind the promoters of lncRNAs downregulated with age. EnrichR was run on PSCAN outputs to find ontologies of transcription factors that bind each promoter, the number of lncRNA promoters that were associated with each GO term is shown.

## 1: Conduct PSCAN

PCAN was conducted with via the web server (http://159.149.160.88/pscan/) to examine which transcription factor binding motifs are enriched in the promoters of lncRNAs that are downregulated with age.
This was conducted in a window -450bp upstream and +50bp downstream of the transcriptional start site for each gene.
A transcription factor binding site was treated as "enriched" with a PSCAN p<0.05.

## 2: Conduct GO via Enrichr

Enrichr was run individually for the transcritption factors that bound each promoter. 

## 3: Plot
