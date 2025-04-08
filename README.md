# GO_cluster


GO_cluster is used to identify overlapping GO terms across multiple GO analyses, wether they are different species or treatment conditions. Additionally, it can determine representative GO terms for shared clusters, providing valuable insights into functional similarities. 

This package is ment to be used to parse and further process the results produced by "clusterProfiler" [1]

Corrylation between GO term dissimilarity between species can be assessed. Average shared GO term dissimilarity between species is used for multispecies clustering.  
The GO term selected to represent a GO cluster is the term with the largest average propotion of genes change between GO analyses. 

## Go cluster binning:  
![Ec_Pa__Dendrogram_MF_2025-02-19](https://github.com/user-attachments/assets/b6510d71-af9f-49d5-877f-e3b2ac9d4716)

## Go cluster differential expression:
![Ec_Pa_MF_Shared_Only_2025-02-19](https://github.com/user-attachments/assets/717224de-5243-49ee-b115-16659d0d15c6)
Y-axis: GO terms, each with the number of significantly differentially expressed genes vs. background in parentheses
X-axis: Log2 fold change in gene expression associated with each GO term.


## ***Single Species Transcriptome Analysis***

Inputs:

1. Deseq output (see example "P_aeruginosa_results.csv" in Data folder)
2. Genome annotation database (org.Xy.eg.db)

Steps:
1. import Deseq2 results
2. Filter Significant Genes
3. use "GO.parse" function to format the enrichGO output into a dataframe
4. run "combine.GOs" to cluster GO terms 
5. visulize with "plot.shared.GOs"
##### *RUN "Single_species_tutorial.25.04.08.R" for in depth exploration of single species analysis!*


## ***Two Species Transcriptome Analysis***

Inputs:
1. Deseq outputs for both species (see example "P_aeruginosa_results.csv" in Data folder)
2. Genome annotation databasesfor both species (org.Xy.eg.db) 

Steps:
1. import Deseq2 results
2. Filter Significant Genes
3. use "GO.parse" function to format the enrichGO output into a dataframe
4. run "combine.GOs" to cluster GO terms 
5. visulize with "plot.shared.GOs"
##### *RUN "Two_species_tutorial.25.04.08.R" for in depth exploration of two species analysis!*

References:

[1] Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012 May;16(5):284-7. doi: 10.1089/omi.2011.0118. Epub 2012 Mar 28. PMID: 22455463; PMCID: PMC3339379.
