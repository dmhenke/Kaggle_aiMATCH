rppa_clean <- function(dp){
  # If the first column contains row identifiers (cell line names), set them as rownames
  if (ncol(dp) > 0 && (names(dp)[1] == "X" || names(dp)[1] == "")) {
    rownames(dp) <- dp[[1]]
    dp <- dp[ , -1, drop = FALSE]
  }
  
  # Remove any "MD." prefix from column names (e.g., "MD.P53" -> "P53")
  rm_MD <- grep("MD.",colnames(dp),fixed=T,value=T)
  if(length(rm_MD) > 0) {
    message("Removing 'MD.' columns from training RPPA columns: ", paste(rm_MD, collapse = ", "))
    # dp_MD <- dp[,grep("^MD\\.", colnames(dp))]
    dp_MD <- dp[,rm_MD]
    dp <- dp[,!colnames(dp)%in%rm_MD]
  } else dp_MD <- NULL
  
  # Clean RPPA columns: Remove any ".{xxxx}" columns (active form of protein)
  rm_active <- grep(".", colnames(dp), value = TRUE,fixed=T)
  if (length(rm_active) > 0) {
    message("Removing active form columns from training RPPA columns: ", paste(rm_active, collapse = ", "))
    dp_Active <- dp[,rm_active]
    dp <- dp[,!colnames(dp)%in%rm_active]
  } else dp_Active <- NULL
  
  list(rppa=dp,rppa_MD=dp_MD,rppa_active = dp_Active)
}

map_ColNams <- function(align_df,map=depmap_id){
  # Build mapping from UniprotID -> Symbol
  # map <- depmap_id
  # Ensure characters
  map$UniprotID <- as.character(map$UniprotID)
  map$Symbol <- as.character(map$Symbol)
  
  # Current DepMap column names are Uniprot IDs; map them to gene symbols where possible
  dep_cols <- names(align_df)
  sym_map <- setNames(map$Symbol, map$UniprotID)
  mapped_names <- sym_map[dep_cols]
  # If any Uniprot IDs lack a mapping, keep the original Uniprot ID as the column name
  mapped_names[is.na(mapped_names)] <- dep_cols[is.na(mapped_names)]
  names(align_df) <- mapped_names
  
  # Now align columns to the training RPPA block: keep order of `train_rppa` columns
  train_cols <- names(train_rppa)
  common_cols <- intersect(train_cols, names(align_df))
  
  # Create aligned DepMap matrix with same column order as train_rppa; fill missing with NA
  dp_aligned <- align_df[ , common_cols, drop = FALSE]
  missing_cols <- setdiff(train_cols, names(dp_aligned))
  if (length(missing_cols) > 0) {
    na_mat <- matrix(NA, nrow = nrow(align_df), ncol = length(missing_cols))
    colnames(na_mat) <- missing_cols
    rownames(na_mat) <- rownames(align_df)
    dp_aligned <- cbind(dp_aligned, na_mat)
  }
  
  # Reorder to exactly match training RPPA column order
  dp_aligned <- dp_aligned[ , train_cols, drop = FALSE]
  
  message("DepMap RPPA columns mapped: total DepMap cols = ", length(dep_cols),
          "; mapped to symbols = ", sum(!(dep_cols %in% mapped_names)), " (kept Uniprot where unmapped)")
  message("Common RPPA columns with training data: ", length(common_cols), " / ", length(train_cols))
  
  # Export aligned object into environment for downstream use
  return(dp_aligned)
  # Optionally write to file (uncomment to save)
  # write.csv(depmap_rppa_aligned, "./Data/depmap_rppa_aligned.csv", row.names = TRUE)
} 

clean_ColNA <- function(x){
  xdf <- map_ColNams(x,map=depmap_id)
  xdf <- xdf[,!sapply(xdf, function(x) all(is.na(x)))]
  xdf
}

# match_ColNams <- function(x, template){
#   
# }

# --- Autoencoder with keras3: impute, scale, build/train, extract embeddings ---
build_and_train_autoencoder <- function(mat, encoding_dim = NULL, epochs = 50, batch_size = 32, validation_split = 0.2, seed = 42) {
  if (!requireNamespace("keras3", quietly = TRUE)) {
    stop("The 'keras3' package is required but not installed. Install with: install.packages('keras3')")
  }
  library(keras3)
  
  set.seed(seed)
  #keras3::op_reset_seed(seed)
  
  # mat: rownames = samples, cols = features
  if (is.data.frame(mat)) mat <- as.data.frame(mat)
  # Convert to numeric matrix and preserve rownames
  rn <- rownames(mat)
  mat_num <- as.matrix(sapply(mat, as.numeric))
  rownames(mat_num) <- rn
  
  # Impute NAs by column medians
  col_medians <- apply(mat_num, 2, function(x) median(x, na.rm = TRUE))
  for (j in seq_len(ncol(mat_num))) {
    nas <- is.na(mat_num[, j])
    if (any(nas)) mat_num[nas, j] <- col_medians[j]
  }
  
  # Scale (center & scale) by column
  col_means <- colMeans(mat_num)
  col_sds <- apply(mat_num, 2, sd)
  col_sds[col_sds == 0] <- 1
  mat_scaled <- scale(mat_num, center = col_means, scale = col_sds)
  
  input_dim <- ncol(mat_scaled)
  if (is.null(encoding_dim)) {
    encoding_dim <- max(2, floor(input_dim / 8))
  }
  
  message("Building autoencoder: input_dim=", input_dim, ", encoding_dim=", encoding_dim)
  
  # Build autoencoder (simple symmetric dense network) with keras3
  inputs <- keras3::layer_input(shape = input_dim)
  encoded <- inputs %>%
    keras3::layer_dense(units = min(256, input_dim), activation = "relu") %>%
    keras3::layer_dense(units = max(encoding_dim, 2), activation = "relu")
  decoded <- encoded %>%
    keras3::layer_dense(units = min(256, input_dim), activation = "relu") %>%
    keras3::layer_dense(units = input_dim, activation = "linear")
  
  autoencoder <- keras3::keras_model(inputs = inputs, outputs = decoded)
  encoder <- keras3::keras_model(inputs = inputs, outputs = encoded)
  
  autoencoder %>% keras3::compile(optimizer = "adam", loss = "mse")
  
  # Early stopping callback
  callbacks_list <- list(
    keras3::callback_early_stopping(monitor = "val_loss", patience = 10, restore_best_weights = TRUE)
  )
  
  message("Starting training...")
  history <- autoencoder %>% keras3::fit(
    x = mat_scaled, y = mat_scaled,
    epochs = epochs,
    batch_size = batch_size,
    validation_split = validation_split,
    callbacks = callbacks_list,
    verbose = 1
  )
  
  # Get latent embeddings using the encoder
  embeddings <- predict(encoder, mat_scaled)
  rownames(embeddings) <- rownames(mat_scaled)
  
  message("Training complete. Final loss: ", round(tail(history$metrics$loss, 1), 4))
  
  list(
    autoencoder = autoencoder,
    encoder = encoder,
    history = history,
    embeddings = embeddings,
    scaling = list(means = col_means, sds = col_sds)
  )
}

