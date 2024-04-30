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
source("Functions.24.04.25.R")
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
# only present the "red" shared go clusters
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

################################################################################
################################################################################

res <- NULL
group <- 0

for(i in 1 : 2){
  
  print(i)
  
  for(j in 1 : 2){
    
    print(j)
    
    go.1 <- enrichGO(gene = get(org1[i]), 
                     OrgDb = org.EcK12.eg.db, 
                     keyType = "REFSEQ",
                     ont = types[types == "BP"],
                     universe = rownames(E.coli), 
                     readable = TRUE)
    
    go.1.parsed <- parse.go(go.result = go.1, 
                            db = org.EcK12.eg.db)
    
    go.2 <- enrichGO(gene = get(org2[j]), 
                     OrgDb = org.PaeruginosaPAO1.eg.db, 
                     keyType = "REFSEQ",
                     ont = types[types == "BP"],
                     universe = rownames(P.aeruginosa), 
                     readable = TRUE)
    
    go.2.parsed <- parse.go(go.result = go.2, 
                            db = org.PaeruginosaPAO1.eg.db)
    
    res.i.j <- get.2gos_onlyshared(go.1 = go.1.parsed,
                                   go.2 = go.2.parsed, 
                                   dseq_results1 = E.coli,
                                   dseq_results2 = P.aeruginosa)
    
    group  <-  group + 1
    
    if(!is.null(res.i.j)){
      res.i.j$group <- group 
    }
    res <- rbind(res, res.i.j)
  }
}

agg <- aggregate(res$group, list(res$group), length)

################################################################################
tiff(paste0("RESULTS/EP_BP_", Sys.Date(),".tiff"),  units="in",  width=6,  height = 8 , res=600)
par(mar = c(4, 15, 4, 2))

bp <- boxplot(res$lfc,
              xlim = c(-1, nrow(res)),
              yaxt = "n",
              horizontal = T, 
              col = res$col, 
              cex.axis=0.4, 
              xlab = "log2 fold change")
abline(h = seq(2 , nrow(res), 2) + .5, col = "grey80")
abline(h = cumsum(agg$x) + .5, col = "grey30", lwd = 2)
abline(v = 0, col = "red", lwd = 2)

legend("bottomright",
       bty = "n", cex = .5, 
       fill = c("cyan4", "firebrick3"), 
       legend = c("P. aeruginosa", "E. coli"), 
       text.font = 4)

axis(2,
     hadj = 0.5,
     font = 3,
     las = 1,
     at=1 : nrow(res),
     labels=paste0(" (n=", res$Count, "/", res$genes.in.term,")"),
     tck = -.01,
     col.ticks ="grey60", 
     col.axis ="grey60", 
     cex.axis=0.4)

res$Description2 <- sapply(res$Description, add_carriage_return)

axis(2, 
     hadj = 1,
     font = 3,
     las = 1,
     at = seq(2 , nrow(res), 2) - .5,
     #labels = paste(res$Description[seq(2 , nrow(res), 2)]),
     labels = paste(res$Description2[seq(2 , nrow(res), 2)], "     "),
     
     tck = -0.01,
     col.ticks ="grey20", 
     col.axis="grey20",
     cex.axis=0.6)

dev.off()
################################################################################
types <- c("BP", "CC", "MF")
org1 <- c("genes1.up", "genes1.down")
org2 <- c("genes2.up", "genes2.down")

B.theta <- B.theta[complete.cases(B.theta) , ]
L.rham <- L.rham[complete.cases(L.rham) , ]

d1.pick <- B.theta[B.theta$padj <= 0.05 , ]
d2.pick <- L.rham[L.rham$padj <= 0.05 , ]

genes1.up <- rownames(d1.pick)[d1.pick$log2FoldChange > 0 ]
genes1.down <- rownames(d1.pick)[d1.pick$log2FoldChange < 0 ]
genes2.up <- rownames(d2.pick)[d2.pick$log2FoldChange > 0 ]
genes2.down <- rownames(d2.pick)[d2.pick$log2FoldChange < 0 ]
################################################################################
res <- NULL
group <- 0

for(i in 1 : 2){
  
  print(i)
  
  for(j in 1 : 2){
    
    print(j)
    
    go.1 <- enrichGO(gene = get(org1[i]), 
                     OrgDb = org.Bthetaiotaomicron.eg.db, 
                     keyType = "REFSEQ",
                     ont = types[types == "BP"],
                     universe = rownames(B.theta), 
                     readable = TRUE)
    
    go.1.parsed <- parse.go(go.result = go.1, 
                            db = org.Bthetaiotaomicron.eg.db)
    
    go.2 <- enrichGO(gene = get(org2[j]), 
                     OrgDb = org.Lrhamnosus.eg.db, 
                     keyType = "REFSEQ",
                     ont = types[types == "BP"],
                     universe = rownames(L.rham), 
                     readable = TRUE)
    
    go.2.parsed <- parse.go(go.result = go.2, 
                            db = org.Lrhamnosus.eg.db)
    
    res.i.j <- get.2gos_onlyshared(go.1 = go.1.parsed,
                                   go.2 = go.2.parsed, 
                                   dseq_results1 = B.theta,
                                   dseq_results2 = L.rham, 
                                   colors = c("goldenrod2", "darkorchid4"))
    
    group  <-  group + 1
    
    if(!is.null(res.i.j)){
      res.i.j$group <- group 
    }
    res <- rbind(res, res.i.j)
  }
}


agg <- aggregate(res$group, list(res$group), length)

################################################################################
tiff(paste0("RESULTS/BL_BP_", Sys.Date(),".tiff"),  units="in",  width=6,  height = 3, res=600)
par(mar = c(4, 15, 4, 2))

bp <- boxplot(res$lfc,
              xlim = c(-.5, nrow(res) + .5),
              yaxt = "n",
              horizontal = T, 
              col = res$col, 
              cex.axis=0.4, 
              xlab = "log2 fold change")
abline(h = seq(2 , nrow(res), 2) + .5, col = "grey80")



legend("bottomright", bty = "n",
       cex = .5, 
       fill =  c("darkorchid4", "goldenrod2"), 
       legend = c("L. rhamnosus", "B. thetaiotaomicron"),
       text.font = 4)

axis(2,
     hadj = 0.5,
     font = 3,
     las = 1,
     at=1 : nrow(res),
     labels=paste0(" (n=", res$Count, "/", res$genes.in.term,")"),
     tck = -.01,
     col.ticks ="grey60", 
     col.axis ="grey60", 
     cex.axis=0.4)

res$Description2 <- sapply(res$Description, add_carriage_return)

axis(2,
     hadj = 0.5,
     font = 3,
     las = 1,
     at=1 : nrow(res),
     labels=paste0(" (n=", res$Count, "/", res$genes.in.term,")"),
     tck = -.01,
     col.ticks ="grey60", 
     col.axis ="grey60", 
     cex.axis=0.4)

axis(2, 
     hadj = 1,
     font = 3,
     las = 1,
     at = seq(2 , nrow(res), 2) - .5,
     labels = paste(res$Description[seq(2 , nrow(res), 2)], "     "),
     tck = -0.01,
     col.ticks ="grey20", 
     col.axis="grey20",
     cex.axis=0.6)

dev.off()

################################################################################



