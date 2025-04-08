library(AnnotationDbi)
library(GenomicFeatures)
library(rtracklayer)
library(tximport)
library(clusterProfiler)

library(org.Bthetaiotaomicron.eg.db)

source("Functions.25.02.19.R")

# =========================
# Step 1: Import DESeq2 Results
# =========================

# Load the DESeq2 results CSV into a dataframe
Btheta.glu_results <- read.csv("Data/Btheta_glucose_results.csv",
                               header = T, row.names = 1)

# remove missing values.
Btheta.glu_results <- Btheta.glu_results[!is.na(Btheta.glu_results$padj) , ]


# =========================
# Step 2: Filter Significant Genes
# =========================

# Keep only rows where adjusted p-value (padj) is ≤ 0.05 (commonly used threshold for statistical significance).
d1.pick <- Btheta.glu_results[Btheta.glu_results$padj <= 0.05 , ]

# Set log2 fold change cutoff to 0 (i.e., any positive or negative change counts).
lcf.cut_off <- 0

# Identify genes that are upregulated (log2FoldChange > 0).
genes1.up <- rownames(d1.pick)[d1.pick$log2FoldChange > lcf.cut_off ]

# Identify genes that are downregulated (log2FoldChange < 0).
genes1.down <- rownames(d1.pick)[d1.pick$log2FoldChange < lcf.cut_off ]


# =========================
# Step 3: GO Enrichment Analysis
# =========================

# Define the GO categories: Biological Process (BP), Cellular Component (CC), and Molecular Function (MF).
types <- c("BP", "CC", "MF")

# Loop through each GO category for enrichment analysis.
for(t in 1 : 3){
  
  # Perform GO enrichment on down regulated genes
  go.1.down <- enrichGO(gene = genes1.down, 
                        OrgDb = org.Bthetaiotaomicron.eg.db, 
                        keyType = "REFSEQ",
                        ont =  types[t],
                        universe = rownames(Btheta.glu_results), 
                        readable = TRUE)
  
  # Perform GO enrichment on up regulated genes
  go.1.up <- enrichGO(gene = genes1.up, 
                      OrgDb = org.Bthetaiotaomicron.eg.db, 
                      keyType = "REFSEQ",
                      ont =  types[t],
                      universe = rownames(Btheta.glu_results), 
                      readable = TRUE)
  
  # Parse the GO results using `GO.parse` and label direction
  go.1.down.parsed <- GO.parse(go.result = go.1.down, 
                               db = org.Bthetaiotaomicron.eg.db)
  go.1.down.parsed$direction <- "down"
  
  go.1.up.parsed <- GO.parse(go.1.up, db = org.Bthetaiotaomicron.eg.db)
  go.1.up.parsed$direction <- "up"
  
  # Combine parsed GO results from up- and down-regulated sets
  go.1 <- rbind2(go.1.down.parsed, 
                 go.1.up.parsed)
  
  par(mar = c(5.1, 5.1, 2, 15))
  
  # Combine and format GO terms using a custom function `combine.GOs`
  # Note: Passing the same inputs twice to make plot for only one species instead of two (tricking it)
  GOs.pick <- combine.GOs(go.1.down.parsed, 
                          go.1.up.parsed, 
                          go.1.down.parsed,
                          go.1.up.parsed,
                          dseq1.pick = d1.pick,
                          dseq2.pick = d1.pick)# this will be fixed for future versions 
  
  GOs.pick <- GOs.pick[seq(1, nrow(GOs.pick), 2), ] #even rows are duplicated (finish tricking it))
  
  par(mar = c(4, 5, 4, 15))
  
  plot.clustered.GO(GOs.pick, 
                    col =  c("goldenrod2"))
  title(types[t])
}




