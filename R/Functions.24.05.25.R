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

GO.parse <- function(go.result, db, use.p.adjust = F){
  
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
      # Print the current iteration index.
      print(i)
      
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

GO.cluster <- function(GO.1.dist, GO.2.dist, similarity.threshold = 0.25){
  
  print(vegan::mantel(go.1.shared.dist, go.2.shared.dist))
  
  plot(as.vector(go.1.shared.dist) ~ 
         as.vector(go.2.shared.dist), pch = 21, bg = "grey")
  
  abline(lm(as.vector(go.1.shared.dist) ~ 
              as.vector(go.2.shared.dist)), col = "red")
  
  go.shared.dist <- analogue::fuse(GO.1.dist, GO.2.dist)
  
  shared.tre <- hclust(go.shared.dist, method = "complete")
  
  # Cut the hierarchical clustering tree to form clusters based on the similarity threshold.
  clusts <- cutree(shared.tre, h = similarity.threshold)
  
  # If the maximum similarity in the distance matrix is greater than the threshold, plot the clustering tree.
  if(max(go.shared.dist) > similarity.threshold & min(go.shared.dist) < 1){
    
    plot(shared.tre, hang = -1, cex = .5, xlab="")
    abline(h = similarity.threshold, col = 2)
  }
  
  return(clusts)
}
