# GO_cluster

R package provides a set of tools for analyzing and visualizing gene ontology (GO) term relationships within and across species. Users can group GO terms based on shared genes within a single species and compare GO term similarities between species. The package further allows for the identification of representative GO clusters, either within a single species or shared between two species, providing insights into conserved or species-specific biological functions. Additionally, the package includes visualization features that plot the distribution of log fold changes (lfcs) of genes within representative GO terms, enabling users to better understand gene expression dynamics. By streamlining complex GO term analyses, GO_cluster offers an invaluable resource for genomic research and comparative functional studies

This package is ment to be used to parse and further process the results produced by "clusterProfiler" [1]

![Bt_Lr_Dendrogram_CC_2025-02-06](https://github.com/user-attachments/assets/e8ecedd2-f65f-4ca7-a3e6-6dba2e04c16f)

GO_cluster stands out for its ability to identify overlapping GO terms across multiple GO analyses, wether they are different species or treatment conditions. Additionally, it can determine representative GO terms for shared clusters, providing valuable insights into functional similarities. 

*Multispecies Analysis*

Corrylation between GO term dissimilarity between species can be assessed. Average shared GO term dissimilarity between species is used for multispecies clustering.  
The GO term selected to represent a GO cluster is the term with the largest average propotion of genes change between GO analyses. 


***Single Species Transcriptome Analysis***

Inputs:

1. Deseq output (see example "P_aeruginosa_results.csv" in Data folder)
2. Genome annotation database (org.Xy.eg.db)

Steps:
1. use "GO.parse" function to format the enrichGO output into a dataframe
2. run "combine.GOs" to cluster GO terms 
3. visulize with "plot.shared.GOs"

*#RUN "Single_species_tutorial.25.04.08.R" for in depth exploration of single species analysis!*


***Two Species Transcriptome Analysis***

Inputs:
1. Deseq outputs for both species (see example "P_aeruginosa_results.csv" in Data folder)
2. Genome annotation databasesfor both species (org.Xy.eg.db) 

Steps:
1. use "GO.parse" function to format each enrichGO output into a dataframe
2. run "combine.GOs" to cluster GO terms from multiple species. 
3. visulize with "plot.shared.GOs"
#*RUN "Two_species_tutorial.25.04.08.R" for in depth exploration of single species analysis!*

References:

[1] Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012 May;16(5):284-7. doi: 10.1089/omi.2011.0118. Epub 2012 Mar 28. PMID: 22455463; PMCID: PMC3339379.
