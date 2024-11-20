#' @title Verify the provided classifier version is valid.
#' @param classifier_version
#' Version of the trained DDI classifier to use. Possible values are 1 or 2.
#' @export
check_classifier_version = function(classifier_version = c(1,2)) {
  if (classifier_version %in% c(1,2)) {
    return(invisible(NULL))
  }
  stop("Unknown version.")
}


#' @title Load the trained classifier version into a dataframe.
#' @inheritParams check_classifier_version
#' @export
load_classifier = function(classifier_version) {
  check_classifier_version(classifier_version)

  classifier_filepath = get_pkg_filepath(
    "classifier", 
    paste0("v", classifier_version), 
    "general", 
    "Classifier.txt"
  )
  
  log_msg(glue::glue("Loading classifier: {classifier_filepath}"))
  result = read.table(
    file = classifier_filepath, 
    header = TRUE, 
    sep = "\t", 
    quote = "", 
    fill = TRUE, 
    comment.char = ""
  )

  # Keep only the necessary columns
  result = result[c("ID", "Genotoxic.score", "Non.Genotoxic.score", "my.train.sd")]

  # Rename the ID column to Probe
  colnames(result)[colnames(result) == 'ID'] = 'Probe'

  return(result)
}


#' @title Load the log10 fold-change data for the training chemicals.
#' @inheritParams check_classifier_version
#' @export
load_training_fc = function(classifier_version) {
  check_classifier_version(classifier_version)
  
  fc_filepath = get_pkg_filepath(
    "classifier", 
    paste0("v", classifier_version), 
    "general", 
    "Heatmap_Data.txt"
  )
  
  log_msg(glue::glue("Loading training fold-change data: {fc_filepath}"))
  read.table(
    file = fc_filepath, 
    header = TRUE,
    check.names = FALSE,
    sep = "\t", 
    quote = "", 
    fill = TRUE, 
    comment.char = ""
  )
}


#' @title Load merged training data.
#' @description
#' Load the specified trained classifier joined to the training chemicals
#' joined by gene symbol.
load_merged_classifier = function(classifier_version) {
  # Load the classifier keeping only necessary columns
  classifier = load_classifier(classifier_version)
  
  # Load the heatmap data
  training_fc = load_training_fc(classifier_version)
  
  # Merge the data frames
  merge(classifier, training_fc, by.x = "Probe", by.y = "Probe")
}


#' @title Calculate genotox score
#' @param probe_values Vector of Log10 fold-change values for each probe for the
#' test chemical.
#' @param probe_scores Vector of training scores for each probe (either genotoxic
#' or non-genotoxic) provided by the classifier data frame.
#' @param probe_score_sds Vector of training standard deviation values for each
#' probe provided by the classifier data frame.
calculate_score = function(probe_values, probe_scores, probe_score_sds) {
  score_differences = probe_values - probe_scores
  score_variances = score_differences / probe_score_sds
  sum((score_variances ^ 2) - 2 * log(0.5))
}


#' @title Validate the input data frame for classification.
#' @description
#' Validate a data frame that should be a merged combination of the result returned
#' from [load_merged_classifier] with the log10 fold change data from a test 
#' chemical. Verifies that the data frame is in the correct shape.
#' @param input_df
#' A data frame that should have the following columns:
#'   * Probe
#'   * Genotoxic.score
#'   * Non.Genotoxic.score
#'   * my.train.sd
#'   * One or more columns with log10 fold-change data for each chemical instance
#'     we'd like to predict probabilites and classification on.
#' @returns A list with the following members:
#'   * df - The unchanged validated data frame
#'   * probe_colname - The name of the column containing the Probe information.
#'   * genotoxic_score_colname - The name of the column containing positive genotox score weights.
#'   * non_genotoxic_score_colname - The name of the column containing negative genotox score weights.
#'   * sd_colname - The name of the column containing standard deviation info for each probe.
#'   * data_colnames - The names of the columns containing the log10 fold change data.
validate_merged_input_df = function(input_df) {
  # Validate input is the correct shape
  expected_headers = c("Probe", "Genotoxic.score", "Non.Genotoxic.score", "my.train.sd")
  if (ncol(input_df) < length(expected_headers)) {
    stop(glue::glue("Expected at least {length(expected_headers)} columns but got {ncol(input_df)}."))
  }
  
  # Validate the headers are as expected
  actual_headers = colnames(input_df)[seq(1, length(expected_headers))]
  header_matches = expected_headers %in% actual_headers
  if (!all(header_matches)) {
    missing_headers = expected_headers(!header_matches)
    missing_headers_str = paste(missing_headers, collapse = ", ")
    stop(glue::glue("One or more expected headers are missing. Headers: {missing_headers_str}"))
  }
  
  # Data starts on column after the expected headers
  start_column = length(expected_headers) + 1
  data_column_indices = seq(start_column, ncol(input_df))
  
  list(
    df = input_df,
    probe_colname = "Probe",
    genotoxic_score_colname = "Genotoxic.score",
    non_genotoxic_score_colname = "Non.Genotoxic.score",
    sd_colname = "my.train.sd",
    data_colnames = colnames(input_df)[data_column_indices]
  )
}


#' @title Calculate the genotox probability and classification.
#' @param probe_scores Vector of genotox positive training scores for each probe.
#' @param probe_scores Vector of genotox negative training scores for each probe.
#' @inheritParams calculate_score
#' @returns 
#' Named list with the calculated probabilities and classification:
#'   * prob_genotoxic
#'   * prob_non_genotoxic
#'   * classification
calculate_probability = function(probe_values, genotox_scores, non_genotox_scores, probe_score_sds) {
  genotox_score = calculate_score(
    probe_values = probe_values,
    probe_scores = genotox_scores,
    probe_score_sds = probe_score_sds
  )
  non_genotox_score = calculate_score(
    probe_values = probe_values,
    probe_scores = non_genotox_scores,
    probe_score_sds = probe_score_sds
  )
  prob_genotoxic = exp(-genotox_score/2)/(exp(-genotox_score/2) + exp(-non_genotox_score/2))
  prob_non_genotoxic = exp(-non_genotox_score/2)/(exp(-genotox_score/2) + exp(-non_genotox_score/2))
  classification = "Unclassified"
  if (prob_genotoxic > 0.9) {
    classification = "Genotoxic"
  } else if (prob_non_genotoxic > 0.9) {
    classification = "Non-Genotoxic"
  }
  
  list(
    prob_genotoxic = prob_genotoxic,
    prob_non_genotoxic = prob_non_genotoxic,
    classification = classification
  )
}


#' @title Predict genotox probabilites and classifications for an input data frame.
#' @param valid_df A list as returned by [validate_merged_input_df] describing
#' the validated input data
#' @returns
#' A data frame with a row for each tested chemical and the following columns:
#'   * prob_genotoxic
#'   * prob_non_genotoxic
#'   * classification
predict_classification = function(valid_df) {
  # Unpack the valid_df input list
  df = valid_df$df
  data_colnames = valid_df$data_colnames
  genotoxic_score_colname = valid_df$genotoxic_score_colname
  non_genotoxic_score_colname = valid_df$non_genotoxic_score_colname
  sd_colname = valid_df$sd_colname
  
  # Build a result data frame
  n_result = length(data_colnames)
  result = data.frame(
    chem_id = data_colnames,
    prob_genotoxic = numeric(n_result),
    prob_non_genotoxic = numeric(n_result),
    classification = character(n_result)
  )
  
  result_i = 1
  for (data_colname in data_colnames) {
    probability_result = calculate_probability(
      probe_values = df[data_colname],
      genotox_scores = df[[genotoxic_score_colname]],
      non_genotox_scores = df[[non_genotoxic_score_colname]],
      probe_score_sds = df[[sd_colname]]
    )  
    
    result[result_i, "prob_genotoxic"] = probability_result[["prob_genotoxic"]]
    result[result_i, "prob_non_genotoxic"] = probability_result[["prob_non_genotoxic"]]
    result[result_i, "classification"] = probability_result[["classification"]]
    result_i = result_i + 1
  }

  return(result)
}


#' @title
#' Create a matrix of just log10 fold change data from the provided data frame.
#' @inheritParams predict_classification
create_fc_matrix = function(valid_df) {
  # Keep just the data columns
  tmp_df = valid_df$df[valid_df$data_colnames]
  
  # Convert to a matrix
  result = as.matrix(tmp_df)
  rownames(result) = valid_df$df[[valid_df$probe_colname]]
  return(result)
}

#' @title Classify a set of test materials as genotoxic, non-genotoxic, or unclassified.
#' @description
#' Classify a set of test materials as being genotoxic, non-genotoxic, or unclassified. 
#' Optionally include the training chemical set in the classification.
#' @param input_data A two column data frame that has gene symbols in column 1
#' and log10 fold-change data in column 2. Column 1 should be named "Probe" and
#' the name of column 2 will be taken as the test material identifier for labeling
#' purposes.
#' @param include_training If TRUE, include the training chemicals in the classification.
#' @inheritParams check_classifier_version
#' @returns A list with the following members:
#'   * fc_matrix - A matrix of just the log10 fold change data for each chemical. Rows are genes and columns are chemicals.
#'   * classification_df - A data frame with the classification results for each test material.
#'       Columns are:
#'         * chem_id - The name of the test material
#'         * prob_genotoxic - The probability of the test material being genotoxic
#'         * prob_non_genotoxic - The probability of the test material being non-genotoxic
#'         * classification - The classification of the test material (Genotoxic, Non-Genotoxic, or Unclassified)
#'   * gene_dist - The gene distance matrix
#'   * gene_clust_dendo - The dendogram of the gene distance matrix
#' @export
classify = function(input_data, classifier_version, include_training = FALSE) {
  # Load the classifier and optionally the training chemicals heatmap data
  if (include_training) {
    dat = load_merged_classifier(classifier_version)
  } else {
    dat = load_classifier(classifier_version)
  }

  # Add the log10 probe data from the input chemicals
  dat = merge(dat, input_data, by.x = "Probe", by.y = "Probe")
  
  # Verify the data frame is in the correct format
  valid_df = validate_merged_input_df(dat)
  
  # Create a matrix of just the log10 fold change data
  fc_matrix = create_fc_matrix(valid_df)
  
  # Obtain classification predictions
  classification_df = predict_classification(valid_df)

  # Dendogram generation
  log_msg("Calculating gene distance and performing dendogram clustering...")
  gene_dist = dist(fc_matrix)
  gene_clust_dendo = hclust(gene_dist, method = "average")
  gene_clust_dendo = as.dendrogram(gene_clust_dendo)
  
  list(
    fc_matrix = fc_matrix,
    classification_df = classification_df,
    gene_dist = gene_dist,
    gene_clust_dendo = gene_clust_dendo
  )
}

save_foldchange_file <- function(out_dir, outfile_prefix, fc_data) {
  filename = glue::glue('{outfile_prefix}_fold_change_data.txt')
  filepath = file.path(out_dir, filename)
  write.table(
    x = data.frame('Data_Label' = rownames(fc_data), fc_data), 
    file = filepath, 
    append = FALSE, 
    quote = FALSE, 
    row.names=FALSE, 
    sep = "\t"
  )
  invisible()
}

save_gene_distance_file <- function(out_dir, outfile_prefix, gene_distance_data) {
  filename = glue::glue('{outfile_prefix}_gene_cluster_dist.txt')
  filepath = file.path(out_dir, filename)
  
  dist_mt = convert_dist_to_matrix(gene_distance_data)
  write.table(
    x = data.frame('Genes' = rownames(dist_mt), dist_mt), 
    file = filepath, 
    append = FALSE, 
    quote = FALSE, 
    row.names=FALSE, 
    sep = "\t"
  )
  invisible()
}


# Convert dist results to matrix and removed upper-right mirrored data
convert_dist_to_matrix <- function(dist_results) {
  dist_matrix <- as.matrix(dist_results)
  
  for (r in 1:length(dist_matrix[,1])) {
    for (c in 1:length(dist_matrix[1,])) {
      if (c > r) dist_matrix[r,c] <- ""
    }
  }
  
  dist_matrix
}