# visualize_latent_space.R
# Loads latent embeddings (CSV produced by `build_and_train_autoencoder`) and
# creates PCA, t-SNE, and UMAP visualizations (PNG + interactive HTML).

visualize_latent_space <- function(embeddings_path = "./Data/depmap_rppa_embeddings.csv",
                                   output_prefix = "./Data/depmap_embeddings",
                                   tsne_perplexity = 30,
                                   random_seed = 42) {
  if (!file.exists(embeddings_path)) stop("Embeddings file not found: ", embeddings_path)

  # Load packages (warn if missing)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
  if (!requireNamespace("Rtsne", quietly = TRUE)) stop("Please install Rtsne")
  if (!requireNamespace("uwot", quietly = TRUE)) stop("Please install uwot")
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Please install plotly")
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) stop("Please install htmlwidgets")

  library(ggplot2)
  library(Rtsne)
  library(uwot)
  library(plotly)
  library(htmlwidgets)

  df <- read.csv(embeddings_path, stringsAsFactors = FALSE)
  # Expect first column to be Sample or similar
  if (ncol(df) < 2) stop("Embeddings file must contain at least one latent dimension column")

  # If there's a 'Sample' column, use it as rownames
  if (tolower(names(df)[1]) %in% c("sample", "samples", "id") ) {
    rownames(df) <- df[[1]]
    df <- df[ , -1, drop = FALSE]
  }

  mat <- as.matrix(df)
  mode(mat) <- "numeric"

  # PCA
  pca <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
  pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Sample = rownames(mat), stringsAsFactors = FALSE)

  p_pca <- ggplot(pca_df, aes(PC1, PC2, label = Sample)) +
    geom_point(alpha = 0.8) + theme_minimal() + ggtitle("PCA of latent embeddings")

  ggsave(paste0(output_prefix, "_pca.png"), p_pca, width = 7, height = 6, dpi = 150)

  # t-SNE (use perplexity smaller of 30 or nrow/3)
  n <- nrow(mat)
  perp <- min(tsne_perplexity, floor((n - 1) / 3))
  set.seed(random_seed)
  tsne_out <- Rtsne::Rtsne(mat, perplexity = max(2, perp), verbose = TRUE, check_duplicates = FALSE)
  tsne_df <- data.frame(TSNE1 = tsne_out$Y[,1], TSNE2 = tsne_out$Y[,2], Sample = rownames(mat), stringsAsFactors = FALSE)

  p_tsne <- ggplot(tsne_df, aes(TSNE1, TSNE2, label = Sample)) +
    geom_point(alpha = 0.8) + theme_minimal() + ggtitle("t-SNE of latent embeddings")

  ggsave(paste0(output_prefix, "_tsne.png"), p_tsne, width = 7, height = 6, dpi = 150)

  # UMAP
  set.seed(random_seed)
  umap_out <- uwot::umap(mat, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
  umap_df <- data.frame(UMAP1 = umap_out[,1], UMAP2 = umap_out[,2], Sample = rownames(mat), stringsAsFactors = FALSE)

  p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, label = Sample)) +
    geom_point(alpha = 0.8) + theme_minimal() + ggtitle("UMAP of latent embeddings")

  ggsave(paste0(output_prefix, "_umap.png"), p_umap, width = 7, height = 6, dpi = 150)

  # Interactive plotly (combine PCA, t-SNE, UMAP into tabs via simple HTML)
  # Create interactive versions for each
  p_pca_int <- plotly::ggplotly(p_pca)
  p_tsne_int <- plotly::ggplotly(p_tsne)
  p_umap_int <- plotly::ggplotly(p_umap)

  # Create a small HTML with the three plots stacked
  html_file <- paste0(output_prefix, "_interactive.html")
  htmlwidgets::saveWidget(htmltools::tagList(p_pca_int, p_tsne_int, p_umap_int), file = html_file, selfcontained = TRUE)

  message("Saved: ", paste0(output_prefix, "_pca.png"), ", ", paste0(output_prefix, "_tsne.png"), ", ", paste0(output_prefix, "_umap.png"))
  message("Interactive HTML saved: ", html_file)

  invisible(list(pca = p_pca, tsne = p_tsne, umap = p_umap,
                 pca_int = p_pca_int, tsne_int = p_tsne_int, umap_int = p_umap_int,
                 files = list(pca = paste0(output_prefix, "_pca.png"),
                              tsne = paste0(output_prefix, "_tsne.png"),
                              umap = paste0(output_prefix, "_umap.png"),
                              interactive = html_file)))
}

# If the file is sourced interactively, do not auto-run; call `visualize_latent_space()` manually.

```