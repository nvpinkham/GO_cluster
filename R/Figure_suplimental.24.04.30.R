################################################################################
#read in data:
library(AnnotationDbi)
# BiocManager::install("org.EcK12.eg.db")
library(org.EcK12.eg.db)
library(org.PaeruginosaPAO1.eg.db)
library(org.Lrhamnosus.eg.db)
library(org.Bthetaiotaomicron.eg.db)
library(DESeq2)
library(GenomicFeatures)
library(rtracklayer)
library(tximport)
library(clusterProfiler)
source("Functions.24.03.20.R")
################################################################################
{gr <- import("gffs/Ecoli.gff")
txdb <- makeTxDbFromGRanges(gr)
k <- keys(txdb, keytype = "CDSNAME")
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "CDSNAME", 
                                 columns = "GENEID")
tx2gene["REFSEQ"] <- tx2gene[,1]
samples <- list.files(path = "Salmon output files/Ecoli", full.names = T,
                      pattern = "quant.sf", recursive = T)
txi <- tximport(samples, type = "salmon", txOut = TRUE,
                ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)
cnt <- txi$counts
colnames(cnt) <- c("CCday12A", "CCday12B", "CCday8A", "CCday8B", "EcMonoA",
                   "EcMonoB")
metadata <- read.csv("Salmon output files/Ec_samples.csv")
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ condition)
dseq2_results <- DESeq(dds)
results_clean <- results(dseq2_results, 
                         contrast = c("condition", "Co", "mono"))
E.coli <- results_clean[complete.cases(results_clean),]

###
###
###
gr <- import("gffs/Paeruginosa.gff")
txdb <- makeTxDbFromGRanges(gr)
k <- keys(txdb, keytype = "CDSNAME")
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "CDSNAME", 
                                 columns = "CDSNAME")
tx2gene["REFSEQ"] <- tx2gene[,1]
samples <- list.files(path = "Salmon output files/Paeruginosa", full.names = T,
                      pattern = "quant.sf", recursive = T)
txi <- tximport(samples, type = "salmon", tx2gene = tx2gene, 
                ignoreAfterBar = TRUE)
cnt <- txi$counts
colnames(cnt) <- c("day12A", "day12B", "day8A", "day8B", "MonoA", "MonoB")
metadata <- read.csv("Salmon output files/Pa_samples.csv")
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ condition)
dseq2_results <- DESeq(dds)
results_clean <- results(dseq2_results, 
                         contrast = c("condition", "Co", "mono"))
P.aeruginosa <- results_clean[complete.cases(results_clean),]
###########
#######
###

gr <- import("gffs/Btheta.gff")
txdb <- makeTxDbFromGRanges(gr)
k <- keys(txdb, keytype = "CDSNAME")
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "CDSNAME", 
                                 columns = "CDSNAME")
#now we'll try a silly hack.

tx2gene["REFSEQ"] <- tx2gene[,1]

#Now we're importing the salmon results from a series of nested folders.

samples_full <- list.files(path = "Salmon output files/Btheta", full.names = T,
                           pattern = "quant.sf", recursive = T)

samples <- samples_full[c(1:2, 7:10)]

txi <- tximport(samples, type = "salmon", tx2gene = tx2gene, 
                ignoreAfterBar = TRUE)

cnt <- txi$counts
colnames(cnt) <- c("48hBtA", "48hBtB",
                   "p4CCA", "p4CCB", "p6CCA", "p6CCB")

metadata <- read.csv("Salmon output files/Bt_no72_samples.csv")
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ condition)
dseq2_results <- DESeq(dds)
results_clean <- results(dseq2_results, 
                         contrast = c("condition", "Co", "Mono"))
B.theta <- results_clean[complete.cases(results_clean),]
#######
gr <- import("gffs/Lrham.gff")
txdb <- makeTxDbFromGRanges(gr)
k <- keys(txdb, keytype = "CDSNAME")
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "CDSNAME", 
                                 columns = "CDSNAME")

tx2gene["REFSEQ"] <- tx2gene[,1]

samples_full <- list.files(path = "Salmon output files/Lrham", full.names = T,
                           pattern = "quant.sf", recursive = T)

samples <- samples_full[c(1:2, 7:10)]

txi <- tximport(samples, type = "salmon", tx2gene = tx2gene, 
                ignoreAfterBar = TRUE)

cnt <- txi$counts

colnames(cnt) <- c("48hLrA", "48hLrB", 
                   "p4CCA", "p4CCB", "p6CCA", "p6CCB")


metadata <- read.csv("Salmon output files/Lr_no72_samples.csv")
metadata <- read.csv("csv files/metadata_LR.csv")

dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ condition)

dseq2_results <- DESeq(dds)

results_clean <- results(dseq2_results, 
                         contrast = c("condition", "Co", "Mono"))
L.rham <- results_clean[complete.cases(results_clean),]
}

################################################################################
types <- c("BP", "CC", "MF")
org1 <- c("genes1.up", "genes1.down")
org2 <- c("genes2.up", "genes2.down")

E.coli <- E.coli[complete.cases(E.coli) , ]
P.aeruginosa <- P.aeruginosa[complete.cases(P.aeruginosa) , ]

d1.pick <- E.coli[E.coli$padj <= 0.05 , ]
d2.pick <- P.aeruginosa[P.aeruginosa$padj <= 0.05 , ]

genes1.up <- rownames(d1.pick)[d1.pick$log2FoldChange > 0 ]
genes1.down <- rownames(d1.pick)[d1.pick$log2FoldChange < 0 ]
genes2.up <- rownames(d2.pick)[d2.pick$log2FoldChange > 0 ]
genes2.down <- rownames(d2.pick)[d2.pick$log2FoldChange < 0 ]

dir.create("Supplimental")

for(i in 1 : 2){
  
  print(i)
  
  for(j in 1 : 2){
    
    print(j)
    
    
    go.1 <- enrichGO(gene = get(org1[i]), 
                     OrgDb = org.EcK12.eg.db, 
                     keyType = "REFSEQ",
                     ont = types[1],
                     universe = rownames(E.coli), 
                     readable = TRUE)
    
    go.1.parsed <- parse.go(go.result = go.1, 
                            db = org.EcK12.eg.db, 
                            use.p.adjust = T)
    
    go.2 <- enrichGO(gene = get(org2[j]), 
                     OrgDb = org.PaeruginosaPAO1.eg.db, 
                     keyType = "REFSEQ",
                     ont = types[1],
                     universe = rownames(P.aeruginosa), 
                     readable = TRUE)
    
    go.2.parsed <- parse.go(go.result = go.2, 
                            db = org.PaeruginosaPAO1.eg.db, 
                            use.p.adjust = T)
    
    h <-  sum(max.e0(go.1.parsed$go.cluster), 
              max.e0(go.2.parsed$go.cluster))
    
    tiff(paste0("Supplimental/EP_", types[1], "_", i, j, "_", Sys.Date(),".tiff"), 
         units="in", 
         width=8, 
         height = max(c(h / 3), 6) , 
         res=600)
    
    par(mar = c(5, 25, 5, 1))
    
    barplot.2gos(go.1 = go.1.parsed,
                 go.2 = go.2.parsed, 
                 dseq_results1 = E.coli,
                 dseq_results2 = P.aeruginosa)
    title(types[k], font.main = 4)
    dev.off()
    
  }
}


go.2.parsed <- parse.go(go.result = go.2, 
                        db = org.PaeruginosaPAO1.eg.db)


################################################################################
B.theta <- B.theta[complete.cases(B.theta) , ]
L.rham <- L.rham[complete.cases(L.rham) , ]

d1.pick <- B.theta[B.theta$padj <= 0.05 , ]
d2.pick <- L.rham[L.rham$padj <= 0.05 , ]

genes1.up <- rownames(d1.pick)[d1.pick$log2FoldChange > 0 ]
genes1.down <- rownames(d1.pick)[d1.pick$log2FoldChange < 0 ]
genes2.up <- rownames(d2.pick)[d2.pick$log2FoldChange > 0 ]
genes2.down <- rownames(d2.pick)[d2.pick$log2FoldChange < 0 ]

for(i in 1 : 2){
  
  print(i)
  
  for(j in 1 : 2){
    
    print(j)
    
    
    
    go.1 <- enrichGO(gene = get(org1[i]), 
                     OrgDb = org.Bthetaiotaomicron.eg.db, 
                     keyType = "REFSEQ",
                     ont = types[1],
                     universe = rownames(B.theta), 
                     readable = TRUE)
    
    go.1.parsed <- parse.go(go.result = go.1, 
                            db = org.Bthetaiotaomicron.eg.db, 
                            similarity.threshold = .3)
    
    go.2 <- enrichGO(gene = get(org2[j]), 
                     OrgDb = org.Lrhamnosus.eg.db, 
                     keyType = "REFSEQ",
                     ont = types[1],
                     universe = rownames(L.rham), 
                     readable = TRUE)
    
    go.2.parsed <- parse.go(go.result = go.2, 
                            db = org.Lrhamnosus.eg.db, 
                            similarity.threshold = .3)
    
    h  <-  sum(c(max.e0(go.1.parsed$go.cluster), 
                 max.e0(go.2.parsed$go.cluster)))
    
    tiff(paste0("Supplimental/BL_", types[1], "_", i, j, "_", Sys.Date(),".tiff"), 
         units="in", 
         width=8, 
         height = max(c(h / 3), 6) , 
         res=600)
    
    par(mar = c(5, 25, 5, 1))
    
    barplot.2gos(go.1 = go.1.parsed,
                 go.2 = go.2.parsed, 
                 dseq_results1 = B.theta,
                 dseq_results2 = L.rham, 
                 colors = c("goldenrod2", "darkorchid4"))
    try(title(types[k], font.main = 4))
    dev.off()
  }
}






################################################################################
dir.create("RESULTS_padj")

{
  
  d1.pick <- E.coli[E.coli$padj <= 0.05 , ]
  d2.pick <- P.aeruginosa[P.aeruginosa$padj <= 0.05 , ]
  
  genes1.up <- rownames(d1.pick)[d1.pick$log2FoldChange > 0 ]
  genes1.down <- rownames(d1.pick)[d1.pick$log2FoldChange < 0 ]
  genes2.up <- rownames(d2.pick)[d2.pick$log2FoldChange > 0 ]
  genes2.down <- rownames(d2.pick)[d2.pick$log2FoldChange < 0 ]
  
  for(i in 1 : 2){
    
    print(i)
    
    for(j in 1 : 2){
      
      print(j)
      
      for(k in 1 : 3){
        
        print(k)
        
        go.1 <- enrichGO(gene = get(org1[i]), 
                         OrgDb = org.EcK12.eg.db, 
                         keyType = "REFSEQ",
                         ont = types[k],
                         universe = rownames(E.coli), 
                         readable = TRUE)
        
        go.1.parsed <- parse.go(go.result = go.1, 
                                db = org.EcK12.eg.db,
                                use.p.adjust = T)
        
        go.2 <- enrichGO(gene = get(org2[j]), 
                         OrgDb = org.PaeruginosaPAO1.eg.db, 
                         keyType = "REFSEQ",
                         ont = types[k],
                         universe = rownames(P.aeruginosa), 
                         readable = TRUE)
        
        go.2.parsed <- parse.go(go.result = go.2, 
                                db = org.PaeruginosaPAO1.eg.db,
                                use.p.adjust = T)
        
        h <-  sum(max.e0(go.1.parsed$go.cluster), 
                  max.e0(go.2.parsed$go.cluster))
        
        tiff(paste0("RESULTS_padj/EP_", types[k], "_", i, j, "_", Sys.Date(),".tiff"), 
             units="in", 
             width=8, 
             height = max(c(h / 3), 6) , 
             res=600)
        
        par(mar = c(5, 25, 5, 1))
        
        barplot.2gos(go.1 = go.1.parsed,
                     go.2 = go.2.parsed, 
                     dseq_results1 = E.coli,
                     dseq_results2 = P.aeruginosa)
        title(types[k], font.main = 4)
        dev.off()
      }
    }
  }
  ################################################################################
  
  d1.pick <- B.theta[B.theta$padj <= 0.05 , ]
  d2.pick <- L.rham[L.rham$padj <= 0.05 , ]
  
  genes1.up <- rownames(d1.pick)[d1.pick$log2FoldChange > 0 ]
  genes1.down <- rownames(d1.pick)[d1.pick$log2FoldChange < 0 ]
  genes2.up <- rownames(d2.pick)[d2.pick$log2FoldChange > 0 ]
  genes2.down <- rownames(d2.pick)[d2.pick$log2FoldChange < 0 ]
  
  for(i in 1 : 2){
    
    print(i)
    
    for(j in 1 : 2){
      
      print(j)
      
      for(k in 1 : 3){
        
        print(k)
        
        go.1 <- enrichGO(gene = get(org1[i]), 
                         OrgDb = org.Bthetaiotaomicron.eg.db, 
                         keyType = "REFSEQ",
                         ont = types[k],
                         universe = rownames(B.theta), 
                         readable = TRUE)
        
        go.1.parsed <- parse.go(go.result = go.1, 
                                db = org.Bthetaiotaomicron.eg.db,
                                use.p.adjust = T)
        
        go.2 <- enrichGO(gene = get(org2[j]), 
                         OrgDb = org.Lrhamnosus.eg.db, 
                         keyType = "REFSEQ",
                         ont = types[k],
                         universe = rownames(L.rham), 
                         readable = TRUE)
        
        go.2.parsed <- parse.go(go.result = go.2, 
                                db = org.Lrhamnosus.eg.db,
                                use.p.adjust = T)
        
        h  <-  sum(c(max.e0(go.1.parsed$go.cluster), 
                     max.e0(go.2.parsed$go.cluster)))
        
        tiff(paste0("RESULTS_padj/BL_", types[k], "_", i, j, "_", Sys.Date(),".tiff"), 
             units="in", 
             width=8, 
             height = max(c(h / 3), 6) , 
             res=600)
        
        par(mar = c(5, 25, 5, 1))
        
        barplot.2gos(go.1 = go.1.parsed,
                     go.2 = go.2.parsed, 
                     dseq_results1 = B.theta,
                     dseq_results2 = L.rham, 
                     colors = c("goldenrod2", "darkorchid4"))
        try(title(types[k], font.main = 4))
        dev.off()
      }
    }
  }
  ############
  ############
  ############
}




