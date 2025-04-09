get.lfc <- function(go, dseq_results){
  #extracts log2 fold-change (log2FC) values for genes associated with specific Gene Ontology (GO) terms,

  go.lfc <- NULL

  if( nrow( go) > 0){
    for(i in 1 : nrow( go)){

      these <-  strsplit(go$symbols[i], "/")[[1]]
      dseq_results.i <- dseq_results[row.names(dseq_results) %in% these , ]
      go.lfc[[i]] <- dseq_results.i$log2FoldChange
    }

    #names(go.lfc) <- paste("\n(n=", comp.pick$Count, "out of", comp.pick$genes.in.term, ")")
    names(go.lfc) <- paste0("(n=", go$Count,")")

    return(go.lfc)
  }else{
    return()
  }
}

#' @export
max.e0 <- function(x){
  tryCatch({
    return(max(x))
  }, error = function(e) {
    return(0)  # Return 0 when an error occurs
  })
}

which.min <- function(x){
  which(x == min(x))[1]
}

add_carriage_return <- function(input_string, n.words = 5) {
  # Split the input string into words
  words <- strsplit(input_string, " ")[[1]]

  # Initialize a vector to store the modified words
  modified_words <- c()

  # Iterate through the words
  for (i in seq_along(words)) {
    # Add the current word to the list of modified words
    modified_words <- c(modified_words, words[i])

    # If the current index is divisible by 4 and is not the last index,
    # add a carriage return after the current word
    if (i %% n.words == 0 && i != length(words)) {
      modified_words <- c(modified_words, "\n")
    }
  }

  # Combine the modified words back into a single string
  modified_string <- paste(modified_words, collapse = " ")
  modified_string <- gsub("\n ", "\n", modified_string)

  return(modified_string)
}

add_newlines <- function(text, width = 35) {
  # Split the text into words
  words <- unlist(strsplit(text, "\\s+"))

  # Initialize result and character counter
  result <- ""
  line_length <- 0

  for (word in words) {
    word_length <- nchar(word)

    # Check if adding the word exceeds the width
    if (line_length + word_length > width) {
      result <- paste0(result, "\n") # Insert a newline
      line_length <- 0
    }

    # Append the word and update the line length
    result <- paste(result, word)
    line_length <- line_length + word_length + 1
  }

  # Trim leading spaces from the final result
  trimws(result)
}

GO.parse <- function(go.result, db, use.p.adjust = F){

  if(is.null(go.result)){
    return(NULL)
  }

  # Extract the sample name from the go.result expression.
  sample <- deparse(substitute(go.result))
  db.used <- deparse(substitute(db))

  # Extract the GO analysis result data frame from go.result.
  go <- go.result@result

  go$symbols <- NA

  for(i in 1 : nrow(go)){
    # Extract gene symbols from gene2Symbol mapping based on gene IDs in the current row.
    symbols.i <- names(go.result@gene2Symbol)[go.result@gene2Symbol %in% strsplit(go$geneID[i], "/")[[1]]]
    # Combine the extracted symbols into a single string and assign it to the 'symbols' column.
    go$symbols[i] <- paste0(symbols.i, collapse = "/")
  }

  go$genes.in.term <- NA

  if(use.p.adjust){
    #
    go.pick <- go[go$p.adjust <= 0.05 , ]
  }else{
    go.pick <- go[go$pvalue <= 0.05 , ]
  }
  # Check if any significant GO terms are found.
  if(nrow(go.pick) > 0){
    # Iterate through each row in the filtered 'go.pick' data frame.
    for(i in 1 : nrow(go.pick)){

      # Get total number of genes in term
      genes.in.term.i <- clusterProfiler::bitr(go.pick$ID[i], fromType="GOALL", toType="SYMBOL", OrgDb=db)

      # Assign the number of genes mapped to the current GO term to the 'genes.in.term' column.
      go.pick$genes.in.term[i] <- nrow(genes.in.term.i)
    }

    # Calculate the proportion of changed genes for each significant GO term.
    go.pick$proportion.changed <- go.pick$Count / go.pick$genes.in.term


    go.pick$sample <- sample
    go.pick$db.used <- db.used

    # Return the annotated and modified 'go.pick' data frame.
    return(go.pick)

  } else {
    # If no significant GO terms are found, print a message.
    print("no significant GO terms found")
  }
}

GO.dist <- function(parsed.GO = parsed.GO){
  # Initialize a matrix 'unique.dist' to store pairwise similarity between GO terms.
  go.dist <- matrix(nrow = nrow(parsed.GO), ncol = nrow(parsed.GO))
  row.names(go.dist) <- parsed.GO$Description
  colnames(go.dist) <- parsed.GO$Description

  # Calculate pairwise similarity between significant GO terms based on shared genes.
  for(i in 1 : nrow(parsed.GO)){
    genes.i <- strsplit(parsed.GO$geneID[i], "/")[[1]]

    for(j in 1 : i){
      shared.j1 <- !genes.i %in% strsplit(parsed.GO$geneID[j], "/")[[1]]
      shared.j2 <- !strsplit(parsed.GO$geneID[j], "/")[[1]] %in% genes.i

      # Calculate the portion shared of the smallest set of genes between two GO terms.
      go.dist[i, j] <- min(c(length(shared.j1[shared.j1]), length(shared.j2[shared.j2]))) /
        min(c(parsed.GO$Count[i], parsed.GO$Count[j]))
    }
  }
  # Convert the matrix to a distance object.
  go.dist <- as.dist(go.dist)
  return(go.dist)
}

GO.cluster <- function(GO.1.dist, GO.2.dist, similarity.threshold = 0.25, plot = T){

  if (missing(GO.2.dist)) {
    GO.2.dist <- GO.1.dist
  }

  print(vegan::mantel(GO.1.dist, GO.2.dist))

  if(plot){

    lm <- lm(as.vector(GO.1.dist) ~  as.vector(GO.2.dist))

    plot(as.vector(GO.1.dist) ~
           as.vector(GO.2.dist), pch = 21, bg = "grey", main = paste("R^2 = ", summary(lm)$r.squared))

    try(abline(lm(as.vector(GO.1.dist) ~  as.vector(GO.2.dist)), col = "red", lty = 2))
  }

  go.shared.dist <- analogue::fuse(GO.1.dist, GO.2.dist)
  #Combines dissimilarities from two or more dissimilarity objects into a single dissimilarity object so that both original dissimilarities contribute equally

  shared.tre <- hclust(go.shared.dist, method = "complete")

  # Cut the hierarchical clustering tree to form clusters based on the similarity threshold.
  clusts <- cutree(shared.tre, h = similarity.threshold)

  if(plot){
    # If the maximum similarity in the distance matrix is greater than the threshold, plot the clustering tree.
    if(max(go.shared.dist) > similarity.threshold & min(go.shared.dist) < 1){

      shared.tre$labels <- sapply( shared.tre$labels, function(x) add_newlines(x))

      s <- .8
      if(length(   shared.tre$labels) > 20){
        s <- .2
      }

      nodePar <- list(lab.cex = .75, pch = c(NA, 19),
                      cex = s, col = "blue")

      shared.tre <- as.dendrogram(shared.tre)

      plot(shared.tre,
           xlab = "GO term dissimilarity",
           nodePar = nodePar,
           horiz = TRUE,
           ylab= "hierarchical clustering\n(complete linkage)")

      abline(v = similarity.threshold, col = 2)

    }
  }

  return(clusts)
}

rbind2 <- function(go.down.parsed, go.up.parsed){

  if(is.null(nrow(go.down.parsed))){
    return(go.up.parsed)
  }
  if(is.null(nrow(go.up.parsed))){
    return(go.down.parsed)
  }
  return(rbind(go.down.parsed, go.up.parsed))
}

#' @export
plot.GO.clusters <- function(clusts){


  agg <- aggregate(GOs.pick$group, list(GOs.pick$group), length)

  bp <- boxplot(GOs.pick$lfc,
                #  xlim = c(1.5, nrow(GOs.pick)),
                yaxt = "n",
                xaxt = "n",
                horizontal = T,
                col = GOs.pick$col,
                bg.pt = GOs.pick$col,
                cex.axis=0.5,
                xlab = "log2 fold change")

  abline(h = seq(4 , nrow(GOs.pick), 2) - 1.5, col = 1, lty = 3)
  abline(h = cumsum(agg$x)[-nrow(agg)] + .5, col = "grey30", lwd = 2)
  abline(v = 0, col = "red", lwd = 2)

  legend("bottomright",
         bty = "n",
         cex = .5,
         fill = cols,
         legend = labels,
         xpd = TRUE,
         # inset = c(-0.05, -0.15),
         text.font = 4)

  axis(2,
       line = 0,
       hadj = 0.5,
       font = 3,
       las = 1,
       at=1 : nrow(GOs.pick),
       labels=paste0(" (n=", GOs.pick$Count, "/", GOs.pick$genes.in.term,")"),
       tck = -.01,
       col.ticks ="grey60",
       col.axis ="grey60",
       cex.axis=0.4)

  axis(1,
       mgp = c(0, 0, 0)  ,
       line = 0,
       font = 3,
       las = 1,
       at=-10:10,
       labels=-10:10,
       tck = -.01,
       cex.axis=0.4)

  GOs.pick$Description2 <- sapply(GOs.pick$Description.cluster,
                                  function(x) add_newlines(x))
  axis(2,
       mgp = c(3, 2, 0)  ,
       font = 3,
       las = 1,
       at = seq(2 , nrow(GOs.pick), 2) - .5,
       labels = GOs.pick$Description2[seq(2 , nrow(GOs.pick), 2)],
       tck = -0.1,
       col.ticks ="grey20",
       col.axis="grey20",
       cex.axis=0.8)
}

combine.GOs <- function(go.1.down.parsed,
                        go.1.up.parsed,
                        go.2.down.parsed,
                        go.2.up.parsed,
                        dseq1.pick = d1.pick,
                        dseq2.pick = d2.pick,
                        only.shared = T){

  go.1 <- rbind2(go.1.down.parsed,
                 go.1.up.parsed)

  go.2 <- rbind2(go.2.down.parsed,
                 go.2.up.parsed)

  go.1.shared <- go.1[go.1$ID %in% go.2$ID , ]

  if(nrow(go.1.shared) == 0){
    return("no shared")
  }

  go.2.shared <- go.2[go.2$ID %in% go.1$ID , ]

  go.2.shared <- go.2.shared[match(go.1.shared$ID, go.2.shared$ID) , ]

  go.1.shared.dist <- GO.dist(parsed.GO = go.1.shared)
  go.2.shared.dist <- GO.dist(parsed.GO = go.2.shared)

  GO.clusters <- GO.cluster(GO.1.dist = go.1.shared.dist,
                            GO.2.dist = go.2.shared.dist)

  go.1.shared <- go.1.shared[match(names(GO.clusters), go.1.shared$Description) , ]
  go.2.shared <- go.2.shared[match(names(GO.clusters), go.2.shared$Description) , ]

  #all(names(GO.clusters) == go.1.shared$Description)
  #all(names(GO.clusters) == go.2.shared$Description)

  go.1.shared$cluster <- GO.clusters
  go.2.shared$cluster <- GO.clusters

  go.1.shared$lfc <- get.lfc(go.1.shared, dseq1.pick)
  go.2.shared$lfc <- get.lfc(go.2.shared, dseq2.pick)

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

  go.1.shared$Description.cluster <- NA
  go.2.shared$Description.cluster <- NA

  go.1.shared$Description.cluster[res] <- discriptions.long
  go.2.shared$Description.cluster[res] <- discriptions.long

  o <- order(paste(go.1.shared$sample, go.2.shared$sample))
  go.1.shared <- go.1.shared[o,]
  go.2.shared <- go.2.shared[o,]

  go.1.shared$order <- 1 : nrow(go.1.shared)
  go.2.shared$order <- 1 : nrow(go.2.shared)

  GOs.shared <- rbind(go.1.shared,
                      go.2.shared)

  GOs.shared <- GOs.shared[order(GOs.shared$order), ]
  GOs.pick <- GOs.shared[GOs.shared$representative.go.term, ]

  g1 <- seq(1, nrow(GOs.pick), 2)
  g2 <- seq(2, nrow(GOs.pick), 2)

  GOs.pick$group <- NA
  GOs.pick$group[g1] <- paste(GOs.pick$direction[g1],
                              GOs.pick$direction[g2])

  GOs.pick$group[g2] <- paste(GOs.pick$direction[g1],
                              GOs.pick$direction[g2])

  if(only.shared){

    GOs.pick$only.shared <- only.shared

    return(GOs.pick)

  }else{

    GOs.pick$order <- NULL
    GOs.pick$cluster <- GOs.pick$cluster + 20000

    go.1.dist <- GO.dist(parsed.GO = go.1)
    go.1.clusts <- GO.cluster(go.1.dist, plot = T)
   # match(names(go.1.clusts), go.1$Description)# looks good
    go.1$cluster <- go.1.clusts
    shared.clusters.1 <- go.1$cluster[go.1$ID %in% GOs.pick$ID]
    go.1.pick <- go.1[ !go.1$cluster  %in% shared.clusters.1 , ]
    go.1.pick$lfc <- get.lfc(go.1.pick, dseq1.pick)

    go.1.pick$representative.go.term <- F

    res <- NULL
    discriptions.long <- NULL

    clusters <- unique(go.1.pick$cluster)
    for(i in 1 : length(clusters)){

      p <- go.1.pick$cluster == clusters[i]

      this <- which(go.1.pick$proportion.changed == max(go.1.pick$proportion.changed[p]) & p)
      res[i] <- this[1]

      discriptions.long[i] <- paste(go.1.pick$Description[this], collapse = " OR ")
    }

    go.1.pick$representative.go.term[res] <- T
    go.1.pick$Description.cluster <- NA
    go.1.pick$Description.cluster[res] <- discriptions.long

    go.1.pick <- go.1.pick[order(go.1.pick$direction, decreasing = F) , ]
    go.1.pick$cluster <- go.1.pick$cluster + 10000
    go.1.pick$group <- go.1.pick$direction

    go.2.dist <- GO.dist(parsed.GO = go.2)
    go.2.clusts <- GO.cluster(go.2.dist, plot = T)
    # match(names(go.2.clusts), go.2$Description)# looks good
    go.2$cluster <- go.2.clusts
    shared.clusters.2 <- go.2$clust[go.2$ID %in% GOs.pick$ID]
    go.2.pick <- go.2[ !go.2$cluster  %in% shared.clusters.2 , ]
    go.2.pick$lfc <- get.lfc(go.2.pick, dseq2.pick)

    go.2.pick$representative.go.term <- F

    res <- NULL
    discriptions.long <- NULL

    clusters <- unique(go.2.pick$cluster)
    for(i in 1 : length(clusters)){

      p <- go.2.pick$cluster == clusters[i]

      this <- which(go.2.pick$proportion.changed == max(go.2.pick$proportion.changed[p]) & p)
      res[i] <- this[1]

      discriptions.long[i] <- paste(go.2.pick$Description[this], collapse = " OR ")
    }

    go.2.pick$representative.go.term[res] <- T
    go.2.pick$Description.cluster <- NA
    go.2.pick$Description.cluster[res] <- discriptions.long

    go.2.pick <- go.2.pick[order(go.2.pick$direction, decreasing = T) , ]
    go.2.pick$cluster <- go.2.pick$cluster + 30000
    go.2.pick$group <- go.2.pick$direction

    colnames(go.1.pick) == colnames(GOs.pick)

    GOs <- rbind(go.1.pick, GOs.pick, go.2.pick)

    GOs <- GOs[GOs$representative.go.term , ]

    return(GOs)
  }
}

combine.GOs.1 <- function(go.1.down.parsed,
                          go.1.up.parsed,
                          dseq1.pick = d1.pick){

  go.1 <- rbind2(go.1.down.parsed,
                 go.1.up.parsed)

  go.1.dist <- GO.dist(parsed.GO = go.1)

  GO.clusters <- GO.cluster(GO.1.dist = go.1.dist,
                            GO.2.dist = go.1.dist)


  go.1$cluster <- GO.clusters

  go.1$lfc <- get.lfc(go.1, dseq1.pick)

  porp.means <- go.1$proportion.changed


  go.1$representative.go.term <- F

  res <- NULL
  discriptions.long <- NULL

  for(i in 1 : length(unique(GO.clusters))){

    p <- GO.clusters == i

    res[i] <- which(porp.means == max(porp.means[p]) & p)[1]

    discriptions.long[i] <- paste(names(GO.clusters[which(porp.means == max(porp.means[p]) & p)]), collapse = " OR ")
  }

  go.1$representative.go.term[res] <- T

  go.1$Description.cluster <- NA
  go.1$Description.cluster[res] <- discriptions.long

  GOs.pick <- go.1[go.1$representative.go.term, ]

  g1 <- seq(1, nrow(GOs.pick), 2)
  g2 <- seq(2, nrow(GOs.pick), 2)

  GOs.pick$group <- NA
  GOs.pick$group[g1] <- paste(GOs.pick$direction[g1],
                              GOs.pick$direction[g2])

  GOs.pick$group[g2] <- paste(GOs.pick$direction[g1],
                              GOs.pick$direction[g2])

  return(GOs.pick)
}

#' Plot shared GOs
#'
#' Custom plot method for shared.GOs objects.
#'
#' @param x An object of class 'shared.GOs'.
#' @param ... Additional arguments passed to plotting functions.
#' @exportS3Method plot shared.GOs
plot.shared.GOs <- function(GOs.pick,
                            cols = c("cyan4", "firebrick3"),
                            labels = c("P. aeruginosa", "E. coli")){
  #makes plot
  dbs <- unique(GOs.pick$db.used)
  GOs.pick$col[GOs.pick$db.used == dbs[1]] <- cols[1]
  GOs.pick$col[GOs.pick$db.used == dbs[2]] <- cols[2]

  agg <- aggregate(GOs.pick$group, list(GOs.pick$group), length)

  bp <- boxplot(GOs.pick$lfc,
               # xlim = c(0, nrow(GOs.pick)),
                yaxt = "n",
                xaxt = "n",
                horizontal = T,
                col = GOs.pick$col,
                bg.pt = GOs.pick$col,
                cex.axis=0.5,
                xlab = "log2 fold change")

  abline(h = seq(4 , nrow(GOs.pick), 2) - 1.5, col = 1, lty = 3)
  abline(h = cumsum(agg$x)[-nrow(agg)] + .5, col = "grey30", lwd = 2)
  abline(v = 0, col = "red", lwd = 2)

  axis(2,
       line = 0,
       hadj = 0.5,
       font = 3,
       las = 1,
       at=1 : nrow(GOs.pick),
       labels=paste0(" (n=", GOs.pick$Count, "/", GOs.pick$genes.in.term,")"),
       tck = -.01,
       col.ticks ="grey60",
       col.axis ="grey60",
       cex.axis=0.4)

  axis(1,
       mgp = c(0, 0, 0)  ,
       line = 0,
       font = 3,
       las = 1,
       at=-10:10,
       labels=-10:10,
       tck = -.01,
       cex.axis=0.4)

  GOs.pick$Description2 <- sapply(GOs.pick$Description.cluster,
                                  function(x) add_newlines(x))
  axis(2,
       mgp = c(3, 2, 0)  ,
       font = 3,
       las = 1,
       at = seq(2 , nrow(GOs.pick), 2) - .5,
       labels = GOs.pick$Description2[seq(2 , nrow(GOs.pick), 2)],
       tck = -0.1,
       col.ticks ="grey20",
       col.axis="grey20",
       cex.axis=0.8)

  # Get plot boundaries
  x_limits <- par("usr")[1:2]  # x-axis limits
  y_limits <- par("usr")[3:4]  # y-axis limits

  # Define coordinates for the legend
  x_pos <- x_limits[1]  # Far right of the x-axis
  y_pos <- y_limits[1] - min(c(0.1 * diff(y_limits), .7))  # Below the x-axis

  legend("bottomright",
         bty = "n",
         cex = .5,
         fill = cols,
         legend = labels,
         xpd = TRUE,
         x = x_pos,
         y = y_pos,
         text.font = 4)

}

#' @export
plot.clustered.GO <- function(GOs.pick,
                              col =  "goldenrod2"){

  dbs <- unique(GOs.pick$db.used)
  GOs.pick$col <- col

  agg <- aggregate(GOs.pick$group, list(GOs.pick$group), length)

  bp <- boxplot(GOs.pick$lfc,
                # xlim = c(0, nrow(GOs.pick)),
                yaxt = "n",
                xaxt = "n",
                horizontal = T,
                col = GOs.pick$col,
                bg.pt = GOs.pick$col,
                cex.axis=0.5,
                xlab = "log2 fold change")

  abline(h = seq(2 , nrow(GOs.pick), 1) - .5, col = 1, lty = 3)
  abline(h = cumsum(agg$x)[-nrow(agg)] + .5, col = "grey30", lwd = 2)
  abline(v = 0, col = "red", lwd = 2)

  axis(2,
       line = 0,
       hadj = 0.5,
       font = 3,
       las = 1,
       at=1 : nrow(GOs.pick),
       labels=paste0(" (n=", GOs.pick$Count, "/", GOs.pick$genes.in.term,")"),
       tck = -.01,
       col.ticks ="grey60",
       col.axis ="grey60",
       cex.axis=0.4)

  axis(1,
       mgp = c(0, 0, 0)  ,
       line = 0,
       font = 3,
       las = 1,
       at=-10:10,
       labels=-10:10,
       tck = -.01,
       cex.axis=0.4)

  GOs.pick$Description2 <- sapply(GOs.pick$Description.cluster,
                                  function(x) add_newlines(x))
  axis(4,
       mgp = c(3, 2, 0)  ,
       font = 3,
       las = 1,
       at=1 : nrow(GOs.pick),
       labels = GOs.pick$Description2,
       tck = -0.1,
       col.ticks ="grey20",
       col.axis="grey20",
       cex.axis=0.5)
}

#' @export
plot.GOs <- function(GOs.pick,
                            cols = c("cyan4", "firebrick3"),
                            labels = c("P. aeruginosa", "E. coli")){

  dbs <- unique(GOs.pick$db.used)
  GOs.pick$col[GOs.pick$db.used == dbs[1]] <- cols[1]
  GOs.pick$col[GOs.pick$db.used == dbs[2]] <- cols[2]

  agg <- aggregate(GOs.pick$group, list(GOs.pick$group), length)

  bp <- boxplot(GOs.pick$lfc,
                yaxt = "n",
                xaxt = "n",
                horizontal = T,
                col = GOs.pick$col,
                bg.pt = GOs.pick$col,
                cex.axis=0.5,
                xaxs = "i",
                bty = "l",
                xlab = "log2 fold change")

  here <- seq(1 , nrow(GOs.pick), 1)
  a <- grep(" ", GOs.pick$group)
  b <- a[seq(1, length(a), 2)]
  here <- here[!here %in% b]
  abline(h = here +.5, col = 1, lty = 3)

  abline(v = 0, col = "red", lwd = 2)

  axis(2,
       line = 0,
       hadj = 0.5,
       font = 3,
       las = 1,
       at=1 : nrow(GOs.pick),
       labels=paste0(" (n=", GOs.pick$Count, "/", GOs.pick$genes.in.term,")"),
       tck = -.01,
       col.ticks ="grey60",
       col.axis ="grey60",
       cex.axis=0.4)

  axis(1,
       mgp = c(0, 0, 0)  ,
       line = 0,
       font = 3,
       las = 1,
       at=-10:10,
       labels=-10:10,
       tck = -.01,
       cex.axis=0.4)

  GOs.pick$Description2 <- sapply(GOs.pick$Description.cluster,
                                  function(x) add_newlines(x, width = 40))



  p <- !duplicated(GOs.pick$Description2) & grepl(" ", GOs.pick$group)
  p <- !p

  axis(2,
       line = 0,
       mgp = c(3, 2, 0)  ,
       font = 3,
       las = 1,
       at= (1 : nrow(GOs.pick))[p],
       labels= GOs.pick$Description2[p],
       tck = -0,
       col.ticks ="grey20",
       col.axis="grey20",
       cex.axis=0.8)

  # Get plot boundaries
  x_limits <- par("usr")[1:2]  # x-axis limits
  y_limits <- par("usr")[3:4]  # y-axis limits

  # Define coordinates for the legend
  x_pos <- x_limits[1]  # Far right of the x-axis
  y_pos <- y_limits[1] - min(c(0.1 * diff(y_limits), .7))  # Below the x-axis

  legend("bottomright",
         bty = "n",
         cex = .5,
         fill = cols,
         legend = labels,
         xpd = TRUE,
         x = x_pos,
         y = y_pos,
         text.font = 4)

}

