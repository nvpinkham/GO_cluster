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
source("Functions.24.05.25.R")

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
E.coli <- E.coli[complete.cases(E.coli) , ]
P.aeruginosa <- P.aeruginosa[complete.cases(P.aeruginosa) , ]

d1.pick <- E.coli[E.coli$padj <= 0.05 , ]
d2.pick <- P.aeruginosa[P.aeruginosa$padj <= 0.05 , ]

genes1.up <- rownames(d1.pick)[d1.pick$log2FoldChange > 0 ]
genes1.down <- rownames(d1.pick)[d1.pick$log2FoldChange < 0 ]
genes2.up <- rownames(d2.pick)[d2.pick$log2FoldChange > 0 ]
genes2.down <- rownames(d2.pick)[d2.pick$log2FoldChange < 0 ]

types <- c("BP", "CC", "MF")

go.1.down <- enrichGO(gene = genes1.down, 
                      OrgDb = org.EcK12.eg.db, 
                      keyType = "REFSEQ",
                      ont = types[types == "BP"],
                      universe = rownames(E.coli), 
                      readable = TRUE)

go.1.up <- enrichGO(gene = genes1.up, 
                    OrgDb = org.EcK12.eg.db, 
                    keyType = "REFSEQ",
                    ont = types[types == "BP"],
                    universe = rownames(E.coli), 
                    readable = TRUE)

go.2.down <- enrichGO(gene = genes2.down, 
                      OrgDb = org.PaeruginosaPAO1.eg.db, 
                      keyType = "REFSEQ",
                      ont = types[types == "BP"],
                      universe = rownames(P.aeruginosa), 
                      readable = TRUE)

go.2.up <- enrichGO(gene = genes2.up, 
                    OrgDb = org.PaeruginosaPAO1.eg.db, 
                    keyType = "REFSEQ",
                    ont = types[types == "BP"],
                    universe = rownames(P.aeruginosa), 
                    readable = TRUE)
##################################################

go.1.down.parsed <- GO.parse(go.1.down, db = org.EcK12.eg.db)
go.1.down.parsed$direction <- "down"

go.1.up.parsed <- GO.parse(go.1.up, db = org.EcK12.eg.db)
go.1.up.parsed$direction <- "up"

go.2.down.parsed <- GO.parse(go.2.down, db = org.PaeruginosaPAO1.eg.db)
go.2.down.parsed$direction <- "down"

go.2.up.parsed <- GO.parse(go.2.up, db = org.PaeruginosaPAO1.eg.db)
go.2.up.parsed$direction <- "up"

#####################

go.1 <- rbind(go.1.down.parsed, 
              go.1.up.parsed)

go.2 <- rbind(go.2.down.parsed,
              go.2.up.parsed)

go.1$col <- "firebrick3"
go.2$col <- "cyan4"

go.1.shared <- go.1[go.1$ID %in% go.2$ID , ]
go.2.shared <- go.2[go.2$ID %in% go.1$ID , ]

go.2.shared <- go.2.shared[match(go.1.shared$ID, go.2.shared$ID) , ]

all(go.1.shared$ID == go.2.shared$ID)

go.1.shared.dist <- GO.dist(parsed.GO = go.1.shared)
go.2.shared.dist <- GO.dist(parsed.GO = go.2.shared)

tiff(paste0("RESULTS/Ec_Pa_Dendrogram_", Sys.Date(),".tiff"),  units="in",  width=6,  height = 8 , res=600)
par(mar = c(15, 4, 4, 2))

GO.clusters <- GO.cluster(go.1.shared.dist, go.2.shared.dist)

dev.off() 

go.1.shared <- go.1.shared[match(names(GO.clusters), go.1.shared$Description) , ]
go.2.shared <- go.2.shared[match(names(GO.clusters), go.2.shared$Description) , ]

all(names(GO.clusters) == go.1.shared$Description)
all(names(GO.clusters) == go.2.shared$Description)

go.1.shared$cluster <- GO.clusters
go.2.shared$cluster <- GO.clusters

go.1.shared$lfc <- get.lfc(go.1.shared, d1.pick)
go.2.shared$lfc <- get.lfc(go.2.shared, d2.pick)

porp.means <- rowMeans(cbind(go.1.shared$proportion.changed,  
                             go.2.shared$proportion.changed))

go.1.shared$representative.go.term <- F
go.2.shared$representative.go.term <- F

res <- NULL
discriptions.long <- NULL
for(i in 1 : length(unique(GO.clusters))){
  
  p <- GO.clusters == i
  
  res[i] <- which(porp.means == max(porp.means[p]) & p)[1]
  
  discriptions.long[i] <- paste(names(GO.clusters[which(porp.means == max(porp.means[p]) & p)]), collapse = " OR ")
  
}


go.1.shared$representative.go.term[res] <- T
go.2.shared$representative.go.term[res] <- T

go.1.shared$Description.cluster[res] <- discriptions.long
go.2.shared$Description.cluster[res] <- discriptions.long

o <- order(paste(go.1.shared$sample, go.2.shared$sample))
go.1.shared <- go.1.shared[o,]
go.2.shared <- go.2.shared[o,]

go.1.shared$order <- 1 : nrow(go.1.shared)
go.2.shared$order <- 1 : nrow(go.2.shared)

GOs.pick <- rbind(go.1.shared,
                  go.2.shared)

GOs.pick <- GOs.pick[order(GOs.pick$order), ]
GOs.pick <- GOs.pick[GOs.pick$representative.go.term, ]

g1 <- seq(1, nrow(GOs.pick), 2)
g2 <- seq(2, nrow(GOs.pick), 2)

GOs.pick$group <- NA
GOs.pick$group[g1] <- paste(GOs.pick$direction[g1], 
                            GOs.pick$direction[g2])

GOs.pick$group[g2] <- paste(GOs.pick$direction[g1], 
                            GOs.pick$direction[g2])

agg <- aggregate(GOs.pick$group, list(GOs.pick$group), length)

tiff(paste0("RESULTS/Ec_Pa_", Sys.Date(),".tiff"),  units="in",  width=6,  height = 8 , res=600)
par(mar = c(4, 15, 4, 2))

bp <- boxplot(GOs.pick$lfc,
              xlim = c(0, nrow(GOs.pick)),
              yaxt = "n",
              horizontal = T, 
              col = GOs.pick$col, 
              bg.pt = GOs.pick$col, 
              cex.axis=0.4, 
              xlab = "log2 fold change")
abline(h = seq(2 , nrow(GOs.pick), 2) + .5, col = "grey80")

abline(h = cumsum(agg$x)[-nrow(agg)] + .5, col = "grey30", lwd = 2)
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
     at=1 : nrow(GOs.pick),
     labels=paste0(" (n=", GOs.pick$Count, "/", GOs.pick$genes.in.term,")"),
     tck = -.01,
     col.ticks ="grey60", 
     col.axis ="grey60", 
     cex.axis=0.4)

GOs.pick$Description2 <- sapply(GOs.pick$Description.cluster, add_carriage_return)

axis(2, 
     hadj = 1,
     font = 3,
     las = 1,
     at = seq(2 , nrow(GOs.pick), 2) - .5,
     labels = paste(GOs.pick$Description2[seq(2 , nrow(GOs.pick), 2)], "     "),
     tck = -0.01,
     col.ticks ="grey20", 
     col.axis="grey20",
     cex.axis=0.6)
dev.off()
################################################################################
################################################################################
################################################################################
types <- c("BP", "CC", "MF")
B.theta <- B.theta[complete.cases(B.theta) , ]
L.rham <- L.rham[complete.cases(L.rham) , ]

d1.pick <- B.theta[B.theta$padj <= 0.05 , ]
d2.pick <- L.rham[L.rham$padj <= 0.05 , ]

genes1.up <- rownames(d1.pick)[d1.pick$log2FoldChange > 0 ]
genes1.down <- rownames(d1.pick)[d1.pick$log2FoldChange < 0 ]
genes2.up <- rownames(d2.pick)[d2.pick$log2FoldChange > 0 ]
genes2.down <- rownames(d2.pick)[d2.pick$log2FoldChange < 0 ]

go.1.down <- enrichGO(gene = genes1.down, 
                      OrgDb = org.Bthetaiotaomicron.eg.db, 
                      keyType = "REFSEQ",
                      ont = types[types == "BP"],
                      universe = rownames(B.theta), 
                      readable = TRUE)

go.1.up <- enrichGO(gene = genes1.up, 
                    OrgDb = org.Bthetaiotaomicron.eg.db, 
                    keyType = "REFSEQ",
                    ont = types[types == "BP"],
                    universe = rownames(B.theta), 
                    readable = TRUE)

go.2.down <- enrichGO(gene = genes2.down, 
                      OrgDb = org.Lrhamnosus.eg.db, 
                      keyType = "REFSEQ",
                      ont = types[types == "BP"],
                      universe = rownames(L.rham), 
                      readable = TRUE)

go.2.up <- enrichGO(gene = genes2.up, 
                    OrgDb = org.Lrhamnosus.eg.db, 
                    keyType = "REFSEQ",
                    ont = types[types == "BP"],
                    universe = rownames(L.rham), 
                    readable = TRUE)
##################################################

go.1.down.parsed <- GO.parse(go.1.down, db = org.Bthetaiotaomicron.eg.db)
go.1.down.parsed$direction <- "down"

go.1.up.parsed <- GO.parse(go.1.up, db = org.Bthetaiotaomicron.eg.db)
go.1.up.parsed$direction <- "up"

go.2.down.parsed <- GO.parse(go.2.down, db = org.Lrhamnosus.eg.db)
go.2.down.parsed$direction <- "down"

go.2.up.parsed <- GO.parse(go.2.up, db = org.Lrhamnosus.eg.db)
go.2.up.parsed$direction <- "up"

#####################

go.1 <- rbind(go.1.down.parsed, 
              go.1.up.parsed)

go.2 <- rbind(go.2.down.parsed,
              go.2.up.parsed)

go.1$col <- "goldenrod2"
go.2$col <- "darkorchid4"

go.1.shared <- go.1[go.1$ID %in% go.2$ID , ]
go.2.shared <- go.2[go.2$ID %in% go.1$ID , ]

go.2.shared <- go.2.shared[match(go.1.shared$ID, go.2.shared$ID) , ]

all(go.1.shared$ID == go.2.shared$ID)

go.1.shared.dist <- GO.dist(parsed.GO = go.1.shared)
go.2.shared.dist <- GO.dist(parsed.GO = go.2.shared)

tiff(paste0("RESULTS/Bt_Lr_Dendrogram_", Sys.Date(),".tiff"),  units="in",  width=6,  height = 8 , res=600)

GO.clusters <- GO.cluster(go.1.shared.dist, go.2.shared.dist)
unique(GO.clusters)

dev.off()

go.1.shared <- go.1.shared[match(names(GO.clusters), go.1.shared$Description) , ]
go.2.shared <- go.2.shared[match(names(GO.clusters), go.2.shared$Description) , ]

all(names(GO.clusters) == go.1.shared$Description)
all(names(GO.clusters) == go.2.shared$Description)

go.1.shared$cluster <- GO.clusters
go.2.shared$cluster <- GO.clusters

go.1.shared$lfc <- get.lfc(go.1.shared, d1.pick)
go.2.shared$lfc <- get.lfc(go.2.shared, d2.pick)

porp.means <- rowMeans(cbind(go.1.shared$proportion.changed,  
                             go.2.shared$proportion.changed))

go.1.shared$representative.go.term <- F
go.2.shared$representative.go.term <- F

res <- NULL
discriptions.long <- NULL
for(i in 1 : length(unique(GO.clusters))){
  
  p <- GO.clusters == i
  
  res[i] <- which(porp.means == max(porp.means[p]) & p)[1]
  
  discriptions.long[i] <- paste(names(GO.clusters[which(porp.means == max(porp.means[p]) & p)]), collapse = " OR ")
  
}


go.1.shared$representative.go.term[res] <- T
go.2.shared$representative.go.term[res] <- T

go.1.shared$Description.cluster[res] <- discriptions.long
go.2.shared$Description.cluster[res] <- discriptions.long

o <- order(paste(go.1.shared$sample, go.2.shared$sample))
go.1.shared <- go.1.shared[o,]
go.2.shared <- go.2.shared[o,]

go.1.shared$order <- 1 : nrow(go.1.shared)
go.2.shared$order <- 1 : nrow(go.2.shared)

GOs.pick <- rbind(go.1.shared,
                  go.2.shared)

GOs.pick <- GOs.pick[order(GOs.pick$order), ]
GOs.pick <- GOs.pick[GOs.pick$representative.go.term, ]

g1 <- seq(1, nrow(GOs.pick), 2)
g2 <- seq(2, nrow(GOs.pick), 2)

GOs.pick$group <- NA
GOs.pick$group[g1] <- paste(GOs.pick$direction[g1], 
                            GOs.pick$direction[g2])

GOs.pick$group[g2] <- paste(GOs.pick$direction[g1], 
                            GOs.pick$direction[g2])

agg <- aggregate(GOs.pick$group, list(GOs.pick$group), length)


tiff(paste0("RESULTS/Bt_Lr_", Sys.Date(),".tiff"),  units="in",  width=6,  height = 3 , res=600)
par(mar = c(4, 15, 4, 2))
bp <- boxplot(GOs.pick$lfc,
              xlim = c(0, nrow(GOs.pick) + 1),
              yaxt = "n",
              horizontal = T, 
              col = GOs.pick$col, 
              cex.axis=0.4, 
              xlab = "log2 fold change")
abline(h = seq(2 , nrow(GOs.pick), 2) + .5, col = "grey80")

abline(h = cumsum(agg$x)[-nrow(agg)] + .5, col = "grey30", lwd = 2)
abline(v = 0, col = "red", lwd = 2)

legend("bottomright", bty = "n",
       cex = .5, 
       fill =  c("goldenrod2", "darkorchid4"), 
       legend = c("B. thetaiotaomicron", "L. rhamnosus"),
       text.font = 4)

axis(2,
     hadj = 0.5,
     font = 3,
     las = 1,
     at=1 : nrow(GOs.pick),
     labels=paste0(" (n=", GOs.pick$Count, "/", GOs.pick$genes.in.term,")"),
     tck = -.01,
     col.ticks ="grey60", 
     col.axis ="grey60", 
     cex.axis=0.4)

GOs.pick$Description2 <- sapply(GOs.pick$Description.cluster, add_carriage_return)

axis(2, 
     hadj = 1,
     font = 3,
     las = 1,
     at = seq(2 , nrow(GOs.pick), 2) - .5,
     labels = paste(GOs.pick$Description2[seq(2 , nrow(GOs.pick), 2)], "     "),
     tck = -0.01,
     col.ticks ="grey20", 
     col.axis="grey20",
     cex.axis=0.6)
dev.off()
################################################################################





