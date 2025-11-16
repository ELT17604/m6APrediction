#' One-hot encode DNA sequences
#'
#' @description An internal helper function to convert a vector of DNA strings
#' (like 5-mers) into a one-hot encoded data frame suitable for ML models.
#'
#' @param dna_strings A character vector of DNA sequences. All must be the same length.
#'
#' @return A data.frame where each column represents a nucleotide position
#' (e.g., nt_pos1) and values are factors (A, T, C, G).
#'
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Predict m6A sites for multiple samples
#'
#' @description Predicts m6A probability and status for a data frame of input features
#' using a pre-trained random forest model.
#'
#' @param ml_fit A trained random forest model object (e.g., from `readRDS`).
#' @param feature_df A data.frame containing the required input features:
#' "gc_content", "RNA_type", "RNA_region", "exon_length",
#' "distance_to_junction", "evolutionary_conservation", "DNA_5mer".
#' @param positive_threshold Numeric. The probability threshold to call a "Positive" status.
#'
#' @return The original `feature_df` with two new columns appended:
#' `predicted_m6A_prob` and `predicted_m6A_status`.
#'
#' @examples
#' # 1. Load the model included with the package
#' model_path <- system.file("extdata", "rf_fit.rds", package = "m6APrediction")
#'
#' # Check if model file exists before proceeding
#' if (file.exists(model_path)) {
#'   ml_fit <- readRDS(model_path)
#'
#'   # 2. Load the example data included with the package
#'   data_path <- system.file("extdata", "m6A_input_example.csv", package = "m6APrediction")
#'   if (file.exists(data_path)) {
#'     example_df <- read.csv(data_path)
#'
#'     # 3. Run prediction
#'     predictions <- prediction_multiple(ml_fit, example_df)
#'     print(head(predictions))
#'   }
#' }
#' @importFrom stats predict
#' @import randomForest
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))

  original_df <- feature_df
  feature_df$RNA_type   <- factor(feature_df$RNA_type,   levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  feature_df_encoded <- cbind(feature_df, dna_encoding(feature_df$DNA_5mer))
  pred_probs <- predict(ml_fit, newdata = feature_df_encoded, type = "prob")

  original_df$predicted_m6A_prob <- pred_probs[, "Positive"]
  original_df$predicted_m6A_status <- ifelse(original_df$predicted_m6A_prob > positive_threshold, "Positive", "Negative")
  return(original_df)
}


#' Predict m6A site for a single sample
#'
#' @description Predicts m6A probability and status for a single set of feature values.
#' This is a wrapper for `prediction_multiple`.
#'
#' @param ml_fit A trained random forest model object.
#' @param gc_content Numeric value for GC content.
#' @param RNA_type Character string (e.g., "mRNA", "lincRNA").
#' @param RNA_region Character string (e.g., "CDS", "3'UTR").
#' @param exon_length Numeric value for exon length.
#' @param distance_to_junction Numeric value for distance.
#' @param evolutionary_conservation Numeric value for conservation.
#' @param DNA_5mer Character string for the 5-mer sequence (e.g., "GGACA").
#' @param positive_threshold Numeric. The probability threshold to call a "Positive" status.
#'
#' @return A named character vector with two elements: "predicted_m6A_prob"
#' and "predicted_m6A_status".
#'
#' @examples
#' # 1. Load the model included with the package
#' model_path <- system.file("extdata", "rf_fit.rds", package = "m6APrediction")
#'
#' # Check if model file exists before proceeding
#' if (file.exists(model_path)) {
#'   ml_fit <- readRDS(model_path)
#'
#'   # 2. Run prediction with single values
#'   single_pred <- prediction_single(
#'     ml_fit,
#'     gc_content = 0.5,
#'     RNA_type = "mRNA",
#'     RNA_region = "CDS",
#'     exon_length = 10,
#'     distance_to_junction = 8,
#'     evolutionary_conservation = 0.5,
#'     DNA_5mer = "GGACA"
#'   )
#'   print(single_pred)
#' }
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){

  input_df <- data.frame(
    gc_content = as.numeric(gc_content),
    RNA_type = as.character(RNA_type),
    RNA_region = as.character(RNA_region),
    exon_length = as.numeric(exon_length),
    distance_to_junction = as.numeric(distance_to_junction),
    evolutionary_conservation = as.numeric(evolutionary_conservation),
    DNA_5mer = as.character(DNA_5mer),
    stringsAsFactors = FALSE
  )

  result_df <- prediction_multiple(ml_fit, input_df, positive_threshold)

  prob_value <- result_df$predicted_m6A_prob[1]
  status_value <- result_df$predicted_m6A_status[1]

  returned_vector <- c(
    "predicted_m6A_prob" = prob_value,
    "predicted_m6A_status" = status_value
  )

  return(returned_vector)
}
