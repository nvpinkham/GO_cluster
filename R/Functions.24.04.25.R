parse.go <- function(go.result, db, similarity.threshold = 0.3, use.p.adjust = F){
  
  # Extract the sample name from the go.result expression.
  sample <- deparse(substitute(go.result))
  
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
      # Print the current iteration index.
      print(i)
      
      # Map GO term IDs to gene symbols using the 'bitr' function from clusterProfiler package.
      genes.in.term.i <- clusterProfiler::bitr(go.pick$ID[i], fromType="GOALL", toType="SYMBOL", OrgDb=db) 
      
      # Assign the number of genes mapped to the current GO term to the 'genes.in.term' column.
      go.pick$genes.in.term[i] <- nrow(genes.in.term.i)
    }
    
    # Calculate the proportion of changed genes for each significant GO term.
    go.pick$proportion.changed <- go.pick$Count / go.pick$genes.in.term
    
    if(nrow(go.pick) > 1){
      # Initialize a matrix 'unique.dist' to store pairwise similarity between GO terms.
      unique.dist <- matrix(nrow = nrow(go.pick), ncol = nrow(go.pick))
      row.names(unique.dist) <- go.pick$Description
      colnames(unique.dist) <- go.pick$Description
      
      # Calculate pairwise similarity between significant GO terms based on shared genes.
      for(i in 1 : nrow(go.pick)){
        genes.i <- strsplit(go.pick$geneID[i], "/")[[1]]
        
        for(j in 1 : i){
          shared.j1 <- !genes.i %in% strsplit(go.pick$geneID[j], "/")[[1]]
          shared.j2 <- !strsplit(go.pick$geneID[j], "/")[[1]] %in% genes.i
          
          # Calculate the portion shared of the smallest set of genes between two GO terms.
          unique.dist[i, j] <- min(c(length(shared.j1[shared.j1]), length(shared.j2[shared.j2]))) / 
            min(c(go.pick$Count[i], go.pick$Count[j]))
        }
      }
      
      # Convert the matrix to a distance object.
      unique.dist <- as.dist(unique.dist)
      
      # Perform hierarchical clustering on the distance matrix.
      unique.tre <- hclust(unique.dist, method = "complete")
      
      # Cut the hierarchical clustering tree to form clusters based on the similarity threshold.
      clusts <- cutree(unique.tre, h = similarity.threshold)
      
      # If the maximum similarity in the distance matrix is greater than the threshold, plot the clustering tree.
      if(max(unique.dist) > similarity.threshold & min(unique.dist) < 1){
        
        plot(unique.tre, hang = -1, cex = .5)
        abline(h = similarity.threshold, col = 2)
      }
      
      # Assign the cluster labels to the 'go.cluster' column.
      go.pick$go.cluster <- clusts
    }else{
      go.pick$go.cluster <- 1
    }
    # Extract unique clusters from the 'go.cluster' column.
    clusters <- unique(go.pick$go.cluster)
    
    # Initialize a vector 'keep' to store representative GO terms for each cluster.
    keep <- NULL
    
    # Iterate through each cluster.
    for(i in clusters){
      # Extract GO terms from the current cluster.
      go.pick.i <- go.pick[go.pick$go.cluster == i , ]
      
      # Select the GO term with the maximum proportion of changed genes as the representative.
      keep[i] <- go.pick.i$ID[which(go.pick.i$proportion.changed == max(go.pick.i$proportion.changed))[1]]
    }
    
    # Mark the representative GO terms in the 'representative.go.term' column.
    go.pick$representative.go.term <- go.pick$ID %in% keep
    # Assign the sample name to the 'sample' column.
    go.pick$sample <- sample
    
    # Return the annotated and modified 'go.pick' data frame.
    return(go.pick)
    
  } else {
    # If no significant GO terms are found, print a message.
    print("no significant GO terms found")
  }
}

get.lfc <- function(go, dseq_results){
  
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

max.e0 <- function(x){
  tryCatch({
    return(max(x))
  }, error = function(e) {
    return(0)  # Return 0 when an error occurs
  })
}

barplot.1go <- function(go.1, dseq_results1, color = "firebrick3"){
  
  id.1 <- deparse(substitute(dseq_results1))
  
  go.1$lfc <- get.lfc(go.1, dseq_results1)
  
  go.1 <- go.1[go.1$representative.go.term , ]
  
  bp <- boxplot(go.1$lfc,
                xlim = c(-2, nrow(go.1) + 1),
                yaxt = "n",
                horizontal = T, 
                col = color, 
                cex.axis=0.4, 
                xlab = "log2 fold change")
  
  axis(2,
       hadj = 0.5,
       font = 3,
       las = 1,
       at=1 : nrow(go.1),
       labels=paste0(" (n=", go.1$Count, "/", go.1$genes.in.term,")"),
       tck = -.01,
       col.ticks ="grey60", 
       col.axis ="grey60", 
       cex.axis=0.4)
  
  axis(2, 
       hadj = 1,
       font = 3,
       las = 1,
       at=1:nrow(go.1),
       labels=paste(go.1$Description, "     "),
       tck = -0.01,
       col.ticks ="grey20", 
       col.axis="grey20",
       cex.axis=0.6)
}

barplot.2gos <- function(go.1 = go.1.dn.p,
                         go.2 = go.2.up.p, 
                         dseq_results1,
                         dseq_results2, 
                         colors = c("firebrick3", "cyan4")){
  
  # parse go 1st to culster terms!!
  id.1 <- deparse(substitute(dseq_results1))
  id.2 <- deparse(substitute(dseq_results2))
  
  if(identical(go.1, go.2) & identical(go.1, "no significant GO terms found")){
    return("no significant GO terms found")
  }
  if(identical(go.1, "no significant GO terms found")){
    barplot.1go(go.1 = go.2, 
                dseq_results1 = dseq_results2, 
                color = colors[2])
    legend("bottomright", 
           bty = "n",
           fill = colors[2], 
           legend =  id.2)
    return("only 1 go")
  }
  if(identical(go.2, "no significant GO terms found")){
    barplot.1go(go.1 = go.1, 
                dseq_results1 = dseq_results1)
    legend("bottomright",
           bty = "n",
           fill = colors[1], 
           legend =  id.1)
    return("only 1 go")
  }
  
  go.1$lfc <- get.lfc(go.1, dseq_results1)
  go.2$lfc <- get.lfc(go.2, dseq_results2)
  
  # get go clusters shared between go results
  shared.1 <- unique(go.1$go.cluster[go.1$ID %in% go.2$ID])
  shared.2 <- unique(go.2$go.cluster[go.2$ID %in% go.1$ID])
  
  go.1.unique <- go.1[!go.1$go.cluster %in% shared.1,]
  go.1.unique <- go.1.unique[go.1.unique$representative.go.term,]
  
  go.2.unique <- go.2[!go.2$go.cluster %in% shared.2,]
  go.2.unique <- go.2.unique[go.2.unique$representative.go.term,]
  
  go.1.shared <- go.1[go.1$go.cluster %in% shared.1,]
  go.2.shared <- go.2[go.2$go.cluster %in% shared.2,]
  
  shared.gos <- c(go.1.shared$ID, go.2.shared$ID)
  
  if(length(shared.gos) > 0){
    
    key <- matrix(nrow = length(shared.gos), ncol = 2)
    row.names(key) <- shared.gos
    # Iterate through each shared GO term
    for(i in 1 : length(shared.gos)){
      
      try(key[i, 1] <- go.1$go.cluster[go.1$ID == shared.gos[i]], silent = T)
      try(key[i, 2] <- go.2$go.cluster[go.2$ID == shared.gos[i]], silent = T)
      
    }
    
    key <- key[order(key[,1]),]
    key <- unique(key)
    key <- key[!is.na(rowSums(key)),, drop = FALSE]
    
    shared.clusters <- list()
    
    for(i in 1 : nrow(key)){
      
      a <- key[i,1]
      matches <- key[which(key[,1] == a), 2]
      shared.clusters[[i]] <- subset(key, key[,2] %in% matches)
    }
    
    shared.clusters <- unique(shared.clusters)
    
    res <- NULL
    for(i in 1 : length(shared.clusters)){
      
      go.1.pick <- go.1.shared[go.1.shared$go.cluster %in% shared.clusters[[i]][,1],]
      go.2.pick <- go.2.shared[go.2.shared$go.cluster %in% shared.clusters[[i]][,2],]
      
      res.i <-  rbind(go.1.pick[which(go.1.pick$proportion.changed == max(go.1.pick$proportion.changed))[1] , ] ,
                      go.2.pick[which(go.2.pick$proportion.changed == max(go.2.pick$proportion.changed))[1] , ])
      
      # use the descripotion from the best (highest proportion.changed) to represent the SHARED GO cluster
      res.i$Description <- res.i$Description[which(res.i$proportion.changed == max(res.i$proportion.changed))[1]]
      
      res <- rbind(res,
                   res.i)
      
    }
    
    res <- unique(res)
    
    dividers <- cumsum(c(nrow(go.1.unique), 
                         nrow(res),
                         nrow(go.2.unique)))
    
  }else{
    res <- NULL
    
    dividers <- cumsum(c(nrow(go.1.unique), 
                         0,
                         nrow(go.2.unique)))
  }
  
  aa <- rbind(go.1.unique, 
              res,
              go.2.unique)
  
  group.there <- c(nrow(go.1.unique) > 0,
                   !is.null(res) ,
                   nrow(go.2.unique) > 0)
  
  
  aa$col <- CAF::color.groups(aa$sample, colors)
  
  bp <- boxplot(aa$lfc,
                xlim = c(-3, nrow(aa) + 1),
                yaxt = "n",
                horizontal = T, 
                col = aa$col, 
                cex.axis=0.4, 
                xlab = "log2 fold change")
  
  if(group.there[1]){
    abline(h = 1 : (dividers[1]) + .5, col = "grey80")
  }
  if(group.there[2]){
    abline(h = seq((dividers[1]),  dividers[2], 2) + .5, col = "grey80")
  }
  if(group.there[3]){
    abline(h = dividers[2]:dividers[3] + .5, col = "grey80")
  }
  
  abline(h = dividers + .5, col = "grey30")
  abline(v = 0, col = "red")
  
  legend("bottomright", bty = "n",
         fill = colors, 
         legend = c(id.1, id.2))
  
  axis(2,
       hadj = 0.5,
       font = 3,
       las = 1,
       at=1 : nrow(aa),
       labels=paste0(" (n=", aa$Count, "/", aa$genes.in.term,")"),
       tck = -.01,
       col.ticks ="grey60", 
       col.axis ="grey60", 
       cex.axis=0.4)
  
  if(group.there[1]){
    axis(2, 
         hadj = 1,
         font = 3,
         las = 1,
         at=1:dividers[1],
         labels=paste(go.1.unique$Description, "     "),
         tck = -0.01,
         col.ticks ="grey20", 
         col.axis="grey20",
         cex.axis=0.6)
  }
  if(group.there[2]){
    axis(2, 
         hadj = 1,
         font = 3,
         las = 1,
         at= seq((dividers[1] + 1),  dividers[2], 2) + .5,
         labels= paste0(res$Description[seq(1, nrow(res), 2)], "      "),
         col.axis="red",
         cex.axis=0.6, 
         tck = -0.05, 
         col.ticks = "red")
  }
  if(group.there[3]){
    axis(2, 
         hadj = 1,
         font = 3,
         las = 1,
         at= (dividers[2] + 1) : dividers[3],
         labels=paste0(go.2.unique$Description, "      "),
         tck = -0.01,
         col.ticks ="grey20", 
         col.axis="grey20",
         cex.axis=0.6)
  }
} 

get.2gos_onlyshared <- function(go.1 = go.1.dn.p,
                                go.2 = go.2.up.p, 
                                dseq_results1,
                                dseq_results2, 
                                colors = c("firebrick3", "cyan4")){
  
  if(all(grepl("no significant GO terms found", go.1)) | all(grepl("no significant GO terms found", go.2))){
    return(NULL)
  }
  
  # parse go 1st to culster terms!!
  id.1 <- deparse(substitute(dseq_results1))
  id.2 <- deparse(substitute(dseq_results2))
  
  go.1$lfc <- get.lfc(go.1, dseq_results1)
  go.2$lfc <- get.lfc(go.2, dseq_results2)
  
  # get go clusters shared between go results
  shared.1 <- unique(go.1$go.cluster[go.1$ID %in% go.2$ID])
  shared.2 <- unique(go.2$go.cluster[go.2$ID %in% go.1$ID])
  
  go.1.shared <- go.1[go.1$go.cluster %in% shared.1,]
  go.2.shared <- go.2[go.2$go.cluster %in% shared.2,]
  
  shared.gos <- c(go.1.shared$ID, go.2.shared$ID)
  
  if(length(shared.gos) > 0){
    
    key <- matrix(nrow = length(shared.gos), ncol = 2)
    row.names(key) <- shared.gos
    # Iterate through each shared GO term
    for(i in 1 : length(shared.gos)){
      
      try(key[i, 1] <- go.1$go.cluster[go.1$ID == shared.gos[i]], silent = T)
      try(key[i, 2] <- go.2$go.cluster[go.2$ID == shared.gos[i]], silent = T)
      
    }
    
    key <- key[order(key[,1]),]
    key <- unique(key)
    key <- key[!is.na(rowSums(key)),, drop = FALSE]
    
    shared.clusters <- list()
    
    for(i in 1 : nrow(key)){
      
      a <- key[i,1]
      matches <- key[which(key[,1] == a), 2]
      shared.clusters[[i]] <- subset(key, key[,2] %in% matches)
    }
    
    shared.clusters <- unique(shared.clusters)
    
    res <- NULL
    for(i in 1 : length(shared.clusters)){
      
      go.1.pick <- go.1.shared[go.1.shared$go.cluster %in% shared.clusters[[i]][,1],]
      go.2.pick <- go.2.shared[go.2.shared$go.cluster %in% shared.clusters[[i]][,2],]
      
      res.i <-  rbind(go.1.pick[which(go.1.pick$proportion.changed == max(go.1.pick$proportion.changed))[1] , ] ,
                      go.2.pick[which(go.2.pick$proportion.changed == max(go.2.pick$proportion.changed))[1] , ])
      
      # use the descripotion from the best (highest proportion.changed) to represent the SHARED GO cluster
      res.i$Description <- res.i$Description[which(res.i$proportion.changed == max(res.i$proportion.changed))[1]]
      
      res <- rbind(res,
                   res.i)
      
    }
    
    Descriptions <- unique(res$Description)
    samples <- unique(res$sample)
    res2 <- NULL
    
    for(i in 1 : length(Descriptions)){
      
      res.i1 <- res[res$Description == Descriptions[i] & res$sample == samples[1], ]
      res.i2 <- res[res$Description == Descriptions[i] & res$sample == samples[2], ]
      
      res2 <- rbind(res2, 
                    res.i1[which.min(res.i1$pvalue) , ],
                    res.i2[which.min(res.i2$pvalue) , ])
      
    }    
    
    res2$col <- CAF::color.groups(res2$sample, colors)
    
  }else{
    res2 <- NULL
  }
  
  return(res2)
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
  modified_string <- gsub(" \n ", "     \n", modified_string)
  
  return(modified_string)
}


