
##############
##############Step 1: load software:
library(DESeq2)
library(AnnotationDbi)
library(GenomicFeatures)
library(rtracklayer)
library(tximport)

library(clusterProfiler)
library(org.EcK12.eg.db)
library(org.PaeruginosaPAO1.eg.db)

source("Functions.25.02.19.R")

# =========================
# Step 1: Import DESeq2 Results
# =========================

# Import Data: Deseq results dataframe saved as a csv
E.coli <- read.csv("Data/E_coli_results.csv", header = T, row.names = 1)
P.aeruginosa <- read.csv( "Data/P_aeruginosa_results.csv", header = T, row.names = 1)


# remove missing values.
E.coli <- E.coli[!is.na(E.coli$padj) , ]
P.aeruginosa <- P.aeruginosa[!is.na(P.aeruginosa$padj) , ]

# =========================
# Step 2: Filter Significant Genes
# =========================

# set adjusted p value cutoff:
d1.pick <- E.coli[E.coli$padj <= 0.05 , ]
d2.pick <- P.aeruginosa[P.aeruginosa$padj <= 0.05 , ]

# set log fold change cutoff:
lcf.cut_off <- 0
genes1.up <- rownames(d1.pick)[d1.pick$log2FoldChange > lcf.cut_off ]
genes1.down <- rownames(d1.pick)[d1.pick$log2FoldChange < lcf.cut_off ]
genes2.up <- rownames(d2.pick)[d2.pick$log2FoldChange > lcf.cut_off ]
genes2.down <- rownames(d2.pick)[d2.pick$log2FoldChange < lcf.cut_off ]



# Define the GO categories: Biological Process (BP), Cellular Component (CC), and Molecular Function (MF).
types <- c("BP", "CC", "MF")

# Loop through each GO category for enrichment analysis.
for(t in 1 : 3){
  
  # =========================
  # Step 3: GO Enrichment Analysis 
  # =========================
  
  go.1.down <- enrichGO(gene = genes1.down, 
                        OrgDb = org.EcK12.eg.db, 
                        keyType = "REFSEQ",
                        ont =  types[t],
                        universe = rownames(E.coli), 
                        readable = TRUE)
  
  go.1.up <- enrichGO(gene = genes1.up, 
                      OrgDb = org.EcK12.eg.db, 
                      keyType = "REFSEQ",
                      ont =  types[t],
                      universe = rownames(E.coli), 
                      readable = TRUE)
  
  go.2.down <- enrichGO(gene = genes2.down, 
                        OrgDb = org.PaeruginosaPAO1.eg.db, 
                        keyType = "REFSEQ",
                        ont =  types[t],
                        universe = rownames(P.aeruginosa), 
                        readable = TRUE)
  
  go.2.up <- enrichGO(gene = genes2.up, 
                      OrgDb = org.PaeruginosaPAO1.eg.db, 
                      keyType = "REFSEQ",
                      ont = types[t],
                      universe = rownames(P.aeruginosa), 
                      readable = TRUE)
  
  # Combine parsed GO results from up- and down-regulated sets
  go.1.down.parsed <- GO.parse(go.1.down, db = org.EcK12.eg.db)
  go.1.down.parsed$direction <- "down"
  
  go.1.up.parsed <- GO.parse(go.1.up, db = org.EcK12.eg.db)
  go.1.up.parsed$direction <- "up"
  
  go.2.down.parsed <- GO.parse(go.2.down, db = org.PaeruginosaPAO1.eg.db)
  go.2.down.parsed$direction <- "down"
  
  go.2.up.parsed <- GO.parse(go.2.up, db = org.PaeruginosaPAO1.eg.db)
  go.2.up.parsed$direction <- "up"
  
  #####################
  go.1 <- rbind2(go.1.down.parsed, 
                 go.1.up.parsed)
  go.2 <- rbind2(go.2.down.parsed,
                 go.2.up.parsed)
  
  par(mar = c(5.1, 5.1, 2, 15))
  
  # =========================
  # Step 4: Cluster and combine GO terms
  # =========================
  
  # Combine and format GO terms
  GOs.pick <- combine.GOs(go.1.down.parsed, 
                          go.1.up.parsed, 
                          go.2.down.parsed, 
                          go.2.up.parsed, only.shared = F)
  # only.shared set to TRUE to show only the GO clusters with overlapping terms between the two species

  par(mar = c(4, 15, 4, 2))
  
  # =========================
  # Step 5: Visualize representative clusters 
  # =========================
  
  plot.shared.GOs(GOs.pick, 
                  cols =  c("goldenrod2", "darkorchid4"),
                  labels = c("B. thetaiotaomicron", "L. rhamnosus"))
  
}


