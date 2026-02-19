# KAGGLE _ aiMATCH project ####

# LIBRARIES ####
# Install and load biolasso from GitHub
# install.packages("devtools")
# devtools::install_github("dmhenke/BioPrimeLASSO")

# DATA ####
# Read in data from the Data folder
train_data <- read.csv("./Data/train.csv") # 385 1011
# test_data <- read.csv("./Data/test.csv")
# External DepMAP RPPA data
depmap_rppa <- read.csv("./Data/harmonized_RPPA_CCLE.csv") # 899 145
depmap_id <- read.csv("./Data/uniprot_hugo_entrez_id_mapping.csv") # 16398     3
depmatp <- read.csv("./Data/harmonized_MS_CCLE_Gygi.csv") #  375 12559


# FUNCTIONS ####
source("./Data/aiMATCH_functions.R")

# Identify columns with no variability (zero variance)
# Check both numeric (zero variance) and character/factor (only one unique value)
zero_var_numeric <- names(train_data)[sapply(train_data, function(x) {
  if (is.numeric(x)) {
    var(x, na.rm = TRUE) == 0 | all(is.na(x))
  } else {
    FALSE
  }
})]

zero_var_character <- names(train_data)[sapply(train_data, function(x) {
  if (!is.numeric(x)) {
    length(unique(na.omit(x))) <= 1 | all(is.na(x))
  } else {
    FALSE
  }
})]

all_zero_var <- c(zero_var_numeric, zero_var_character)

# Remove zero variance columns from the training data
train_data <- train_data[, !(names(train_data) %in% all_zero_var)]


# Extract RPPA columns: identify block from 'ABL1' to 'ZAP70' (as present in header)
if (any(names(train_data) == "ABL1") && any(names(train_data) == "ZAP70")) {
  rppa_start <- which(names(train_data) == "ABL1")[1]
  rppa_end <- which(names(train_data) == "ZAP70")[1]
  train_rppa <- train_data[, rppa_start:rppa_end] # 385 446
  message("Extracted RPPA columns: ", rppa_start, " to ", rppa_end, " (", ncol(train_rppa), " columns)")
  # print(head(train_rppa))
} else {
  stop("Could not locate RPPA column markers 'ABL1' and/or 'ZAP70' in train_data column names.")
}

# CLEAN data ####
depmap_rppa_clean <- rppa_clean(depmap_rppa);lapply(depmap_rppa_clean,dim)
depmatp_clean <- rppa_clean(depmatp);lapply(depmatp_clean,dim)
train_rppa_clean <- rppa_clean(train_rppa);lapply(train_rppa_clean,dim)

# lapply(c(depmap_rppa_clean,depmatp_clean,train_rppa_clean),function(x)colnames(x)[grep("MD",colnames(x))])



depmap_rppa_cleaner <- clean_ColNA(depmap_rppa_clean$rppa);dim(depmap_rppa_cleaner) # 899 89
depmatp_cleaner <- clean_ColNA(depmatp_clean$rppa);dim(depmatp_cleaner) # 375 269
train_rppa_cleaner <- clean_ColNA(train_rppa_clean$rppa);dim(train_rppa_cleaner) # 385 364


# Intersect colnames of DepMap and Train





# Train autoencoder on aligned DepMap RPPA if available
if (exists("depmap_rppa_aligned")) {
  message("\n=== Starting Autoencoder Training on DepMap RPPA (aligned) ===")
  message("Samples: ", nrow(depmap_rppa_aligned), ", Features: ", ncol(depmap_rppa_aligned))

  ae_res <- tryCatch(
    build_and_train_autoencoder(depmap_rppa_aligned, encoding_dim = NULL, epochs = 50, batch_size = 32),
    error = function(e) {
      message("Error: ", e$message)
      e
    }
  )

  if (inherits(ae_res, "error")) {
    warning("Autoencoder training failed: ", ae_res$message)
  } else {
    message("\n=== Training Successful ===")
    message("Embeddings shape: ", nrow(ae_res$embeddings), " x ", ncol(ae_res$embeddings))

    # Save embeddings to CSV
    embed_df <- as.data.frame(ae_res$embeddings)
    colnames(embed_df) <- paste0("latent_", seq_len(ncol(embed_df)))
    embed_df <- cbind(Sample = rownames(embed_df), embed_df)
    write.csv(embed_df, "./Data/depmap_rppa_embeddings.csv", row.names = FALSE)
    message("Saved embeddings to ./Data/depmap_rppa_embeddings.csv")

    # Optionally save models (uncomment to use)
    # keras3::save_model(ae_res$autoencoder, "./Data/depmap_autoencoder_model.keras")
    # keras3::save_model(ae_res$encoder, "./Data/depmap_encoder_model.keras")
  }
} else {
  message("depmap_rppa_aligned not found; skipping autoencoder step.")
}

#train_rppa from train_data
tmp <- build_and_train_autoencoder(train_rppa, encoding_dim = NULL, epochs = 50, batch_size = 32)
library(umap)
library(ggplot2)
umap_out<- umap(tmp$embeddings, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
  umap_df <- data.frame(UMAP1 = umap_out$layout[,1], UMAP2 = umap_out$layout[,2], Sample = rownames(tmp$embeddings), stringsAsFactors = FALSE)

  p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = train_data$Diagnosis.1)) +
    geom_point(alpha = 0.6) + theme_minimal() + ggtitle("UMAP of latent embeddings")

  # ggsave(paste0(output_prefix, "_umap.png"), p_umap, width = 7, height = 6, dpi = 150)

  # Interactive HTML plots using plotly
  p_pca_plotly <- plotly::ggplotly(p_pca)
  p_tsne_plotly <- plotly::ggplotly(p_tsne)
  p_umap_plotly <- plotly::ggplotly(p_umap)

  htmlwidgets::saveWidget(p_pca_plotly, paste0(output_prefix, "_pca.html"))
  htmlwidgets::saveWidget(p_tsne_plotly, paste0(output_prefix, "_tsne.html"))
  htmlwidgets::saveWidget(p_umap_plotly, paste0(output_prefix, "_umap.html"))

  message("Visualizations saved with prefix: ", output_prefix)
}