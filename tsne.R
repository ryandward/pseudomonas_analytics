
# tsne_out <- Rtsne(data_grid_matrix, check_duplicates = FALSE, perplexity = 10)
# 
# tsne_out <- data.table(spacer = rownames(data_grid_matrix), Dimension1 = tsne_out$Y[,1], Dimension2 = tsne_out$Y[,2])
# 
# plot  (Dimension2 ~ Dimension1, tsne_out, col = alpha("black"), cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[like(spacer, "Ctrl")], col = alpha("green", 0.5), pch = 20, cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[like(spacer, "pur")], col = alpha("blue", 0.5), pch = 20, cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[like(spacer, "orf")], col = alpha("red", 0.5), pch = 20, cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[spacer %in% island$spacer], col = alpha("cyan", 0.5), pch = 20, cex = 2)

# points(Dimension2 ~ Dimension1, tsne_out[like(strain, "Carsonella")], col = alpha("magenta", 0.5), pch = 20, cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[like(strain, "Wiggles")], col = alpha("magenta", 0.5), pch = 20, cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[like(strain, "Blochmannia")], col = alpha("magenta", 0.5), pch = 20, cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[like(strain, "Buchnera")], col = alpha("magenta", 0.5), pch = 20, cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[like(strain, "Sodalis")], col = alpha("magenta", 0.5), pch = 20, cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[like(strain, "Yersinia")], col = alpha("saddlebrown", 0.5), pch = 20, cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[like(strain, "Klebsiella")], col = alpha("cyan", 0.5), pch = 20, cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[like(strain, "Acinetobacter")], col = alpha("green", 0.5), pch = 20, cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[like(strain, "Pseudomonas")], col = alpha("red", 0.5), pch = 20, cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[like(strain, "Enterobacter")], col = alpha("blue", 0.5), pch = 20, cex = 2)
# points(Dimension2 ~ Dimension1, tsne_out[like(strain, "Escherichia")], col = alpha("black", 0.5), pch = 20, cex = 2)
# 
