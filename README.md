# GO_cluster
R package for clustering gene ontology (GO) terms with shared genes to simplify reporting GO term analysis. 

GO terms 

![Bt_Lr_Dendrogram_CC_2025-02-06](https://github.com/user-attachments/assets/e8ecedd2-f65f-4ca7-a3e6-6dba2e04c16f)


This package is ment to be used to parse and further process the results produced by "clusterProfiler" [1]

GO_cluster stands out for its ability to identify overlapping GO terms across multiple GO analyses, wether they are different species or treatment conditions. Additionally, it can determine representative GO terms for shared clusters, providing valuable insights into functional similarities. 

The GO term selected to represent a GO cluster is the term with the largest average propotion of genes change between GO analyses. 


**Single Species Transcriptome Analysis**
Inputs:
1. Gene feature files for both species
2. Expression quantification data from Salmon output files for both species.



**Two Species Transcriptome Analysis**
Inputs:
1. Gene feature files for both species
2. Expression quantification data from Salmon output files for both species.

References:

[1] Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012 May;16(5):284-7. doi: 10.1089/omi.2011.0118. Epub 2012 Mar 28. PMID: 22455463; PMCID: PMC3339379.
