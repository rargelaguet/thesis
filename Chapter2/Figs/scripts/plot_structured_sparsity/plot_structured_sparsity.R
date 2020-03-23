library(MOFA2)

library(RColorBrewer)
library(gplots)
library(pheatmap)

model <- load_model("/Users/ricard/data/mofa2_simulations/test/hdf5/test.hdf5", remove_inactive_factors = F)

model <- subset_factors(model, factors=1:11)
plot_variance_explained(model)




outdir <- "/Users/ricard/thesis/ricard_thesis/Chapter2/Figs/scripts/plot_structured_sparsity"

# Factors
factors_palette <- colorRampPalette(c("white", "darkgrey"))(n=10)
Z <- get_factors(model, scale = T)[[1]] %>% abs
Z[Z>0.55] <- 0.55; Z[Z<(-0.55)] <- -(0.55)
pdf(paste0(outdir,"/Z.pdf"), height=3, width=6)
pheatmap(Z, cluster_rows=F, cluster_cols=F, show_colnames=F, show_rownames=F, color=factors_palette, legend = F)
dev.off()


# Weights
weights_palette <- colorRampPalette(c("white", "red"))(n=5)

W <- get_weights(model, scale = T, abs=T)
for (m in 1:length(W)) {
  
  Wm <- W[[m]]
  Wm[Wm>0.5] <- 0.5
  
  pdf(sprintf("%s/W_%d.pdf",outdir,m), height=3, width=6)
  pheatmap(Wm, cluster_rows=F, cluster_cols=F, show_colnames=F, show_rownames=F, color=weights_palette, border_color="black", legend=F)
  dev.off()
}
