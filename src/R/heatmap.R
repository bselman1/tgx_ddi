blue_yellow_col_vec = function (nrgcols = 12) {
  k <- trunc(nrgcols/2)
  if (2 * k == nrgcols) {
    r <- c(rep(0, k), seq(0, 1, length = k))
    g <- r
    b <- c(rev(seq(0, 1, length = k)), rep(0, k))
    #change
    #colvec <- rgb(r, g, b, rep(0, 2 * k))
    colvec <- rgb(r, g, b)
  }
  else {
    r <- c(rep(0, k + 1), seq(1/(2 * k - 1), 1, length = k))
    g <- r
    b <- c(rev(seq(1/(2 * k - 1), 1, length = k)), rep(0, 
                                                       k + 1))
    #colvec <- rgb(r, g, b, rep(0, 2 * k + 1))
    colvec <- rgb(r, g, b)
    
  }
  colvec
}

#' @title Build heatmap comparing a single test material against the training materials.
#' @description 
#' This function builds a heatmap comparing a single test material against the training materials.
#' The test material should be the last column of the fc_matrix.
plot_heatmap_against_training = function(fc_matrix, classification_df, gene_clust_dendo, plot_title, timestamp) {
  
  # Insert a blank column between the training materials and the test material which is in the last column
  train_chunk = fc_matrix[, seq(1, ncol(fc_matrix) - 1)]
  blank_col = rep(NA, nrow(fc_matrix))
  test_chunk = fc_matrix[, ncol(fc_matrix), drop=FALSE]
  heat_df = cbind(
    train_chunk, 
    blank_col, 
    test_chunk
  )
  # Remove the column name from the blank column
  colnames(heat_df)[colnames(heat_df) == "blank_col"] = ""
  
  # Capture the classification of the test material (last row)
  classification = tail(classification_df, 1)[["classification"]]
  ddi_label = "NDDI"
  if (classification == "Genotoxic") {
    ddi_label = "DDI"
  }
  
  # Add spacer 
  
  heatmap3::heatmap3(
    x = heat_df,
    Rowv = gene_clust_dendo,
    Colv = NA,
    scale = "row",
    margins = c(15, 10),
    cexRow=0.65,
    cexCol=1,
    col = blue_yellow_col_vec(ncol(heat_df)),
    ColSideColors = rep("white", ncol(heat_df)),
    ColSideLabs = ""
  )
  mtext(plot_title, side=3, line=-3, at=c(0.45), cex=1)
  mtext(
    text = c("DNA Damage Inducing (DDI)", "Non DNA Damage Inducing (NDDI)"),
    side=3,
    line=-4.5,
    at=c(0.24, 0.46),
    adj=0,
    cex=0.9
  )
  mtext(timestamp, side=1, line=3, cex=0.8)
  mtext(ddi_label, side=3, line=-5, at=c(0.74), adj=0, cex=0.7)
  
  #left-v
  segments(0.229, 0.84, 0.229, 0.86, col= 'black', lwd=1)
  #mid-v
  segments(0.452, 0.84, 0.452, 0.86, col= 'black', lwd=1)
  #right-v
  segments(0.727, 0.84, 0.727, 0.86, col= 'black', lwd=1)
  #horizental
  segments(0.23, 0.85, 0.726, 0.85, col= 'black', lwd=1)
}


plot_heatmap = function(fc_matrix, classification_df, gene_clust_dendo, plot_title, timestamp) {

  # Build a df that we will use to display colored labels on top of the heatmap
  class_colors = data.frame(
    classification = c("Genotoxic", "Non-Genotoxic", "Unclassified"),
    color = c("red", "blue", "grey")
  )
  classification_df = merge(classification_df, class_colors, by = "classification")
  
  plot_legend = function() {
    heatmap3::showLegend(legend = class_colors$classification, col = class_colors$color)
  }
  
  heatmap(
    x = fc_matrix[, classification_df[["chem_id"]]],
    Rowv = gene_clust_dendo,
    Colv = NA,
    scale = "row",
    margins = c(10, 5),
    cexRow=0.65,
    cexCol=0.65,
    col = blue_yellow_col_vec(ncol(fc_matrix)),
    ColSideColors = classification_df[["color"]]
  )
}


plot_cluster <- function(fc_matrix, plot_title, timestamp) {
  # Calculate the mean and sd for the training set per gene symbol
  # Every column but the last is the training set
  train_df = fc_matrix[, -ncol(fc_matrix)]
  train_mean = rowMeans(train_df)
  train_sd = apply(train_df, 1, sd)
  
  # Calculate z-score using the training mean and sd from above
  z = (fc_matrix - train_mean) / train_sd
  
  # Keep a reference to just the training portion of the z-scores
  z_train = z[, -ncol(fc_matrix)]
  
  # Do principal component analysis on the training set
  pca = prcomp(
    x = t(z_train),
    retx = TRUE,
    center = FALSE,
    scale. = FALSE
  )
  pc1 = as.vector(pca$rotation[,1]) %*% z
  pc2 = as.vector(pca$rotation[,2]) %*% z

  
  # Draw the PCA plot
  par(mar = c(5,4,4,2), oma=c(1,1,2,1), ps=12)
  layout(matrix(c(1,2), byrow=T, ncol=2))
  small_cex <- 0.5
  large_cex <- small_cex * 3
  
  plot(pc1, pc2, type = "n", main = "", xlab = "PC 1", ylab = "PC 2")
  points(pc1[1:13], pc2[1:13], pch = 20, cex=small_cex, col = "red")
  points(pc1[14:29], pc2[14:29], pch = 20, cex=small_cex, col = "blue")
  points(pc1[30:ncol(z)], pc2[30:ncol(z)], pch = 17, cex=large_cex, col = "dark green")
  
  
  chem_labels = colnames(fc_matrix)
  text(pc1[1:13], pc2[1:13], labels = chem_labels[1:13], col = "red", cex = 0.7, pos = 1)
  text(pc1[14:29], pc2[14:29], labels = chem_labels[14:29], col = "blue", cex = 0.7, pos = 1)
  text(pc1[30:ncol(z)], pc2[30:ncol(z)], labels = chem_labels[30:ncol(z)], col = "dark green", cex = 0.7, pos = 1)
  
  abline(v = 0, col = "red", lwd = 2)
  
  # Draw the dendogram plot
  chem_dist <- dist(t(z))
  dendro <- hclust(chem_dist, method = "average")
  dendro <- as.dendrogram(dendro)
  
  # Assigning the labels of dendrogram object with new colors:
  dendextend::labels_colors(dendro) <- c(rep("red", 13), rep("blue", 16), rep("dark green", 5))[order.dendrogram(dendro)]
  
  par(mar = c(3,1,1,10), cex = 0.7)
  plot(dendro, horiz=T)
  
  # Add labels
  mtext(plot_title, side=3, line=-1, outer=T, adj=0.5, cex=1.5)
  mtext(timestamp, side=1, outer=T, line=-0.5, adj=0.49, cex=0.8)
  mtext("Euclidean distance", side=1, outer=T, line=-4.5, at=(0.7), cex=0.8)
  
  return(chem_dist)
}


save_heatmap_pdf = function(fc_matrix, classification_df, gene_clust_dendo, plot_title, timestamp, out_dir, outfile_prefix) {
  filename = glue::glue('{outfile_prefix}_heatmap.pdf')
  filepath = file.path(out_dir, filename)

  Cairo::CairoPDF(
    file = filepath, 
    width = 11, 
    height = 8.5,
    pointsize = 12, 
    family = "Courier"
  )
  plot_heatmap(
    fc_matrix = fc_matrix,
    classification_df = classification_df,
    gene_clust_dendo = gene_clust_dendo, 
    plot_title = plot_title, 
    timestamp = timestamp
  )
  
  dev.off()
  invisible()
}


save_cluster_pdf = function(fc_matrix, plot_title, timestamp, out_dir, outfile_prefix) {
  filename = glue::glue('{outfile_prefix}_cluster.pdf')
  filepath = file.path(out_dir, filename)

  Cairo::CairoPDF(
    file = filepath, 
    width = 11, 
    height = 8.5,
    pointsize = 12, 
    family = "Courier"
  )
  plot_cluster(
    fc_matrix = fc_matrix,
    plot_title = plot_title, 
    timestamp = timestamp
  )
  
  dev.off()
  invisible()
}


save_heatmap_png = function(fc_matrix, classification_df, gene_clust_dendo, plot_title, timestamp, out_dir, outfile_prefix) {
  filename = glue::glue('{outfile_prefix}_heatmap.png')
  filepath = file.path(out_dir, filename)
  
  res2X = 144
  png(
    filename=filepath, 
    width=11 * res2X, 
    height=8.5 * res2X, 
    res=res2X, 
    type="cairo"
  )
  plot_heatmap(
    fc_matrix = fc_matrix,
    classification_df = classification_df,
    gene_clust_dendo = gene_clust_dendo, 
    plot_title = plot_title, 
    timestamp = timestamp
  )
  dev.off()
  
  invisible()
}


save_heatmap_with_cluster_training_pdf = function(fc_matrix, classification_df, gene_clust_dendo, plot_title, timestamp, out_dir, outfile_prefix) {
  filename = glue::glue('{outfile_prefix}_heatmap.pdf')
  filepath = file.path(out_dir, filename)

  Cairo::CairoPDF(
    file = filepath, 
    width = 11, 
    height = 8.5,
    pointsize = 12, 
    family = "Courier"
  )
  plot_heatmap_against_training(
    fc_matrix = fc_matrix,
    classification_df = classification_df,
    gene_clust_dendo = gene_clust_dendo, 
    plot_title = plot_title, 
    timestamp = timestamp
  )
  plot_cluster(
    fc_matrix = fc_matrix,
    plot_title = plot_title, 
    timestamp = timestamp
  )
  
  dev.off()
  invisible()
}


save_heatmap_training_png <- function(fc_matrix, classification_df, gene_clust_dendo, plot_title, timestamp, out_dir, outfile_prefix) {
  filename = glue::glue('{outfile_prefix}_heatmap.png')
  filepath = file.path(out_dir, filename)
  
  res2X = 144
  png(
    filename=filepath, 
    width=11 * res2X, 
    height=8.5 * res2X, 
    res=res2X, 
    type="cairo"
  )
  plot_heatmap_against_training(
    fc_matrix = fc_matrix,
    classification_df = classification_df,
    gene_clust_dendo = gene_clust_dendo, 
    plot_title = plot_title, 
    timestamp = timestamp
  )
  dev.off()
  
  invisible()
}
