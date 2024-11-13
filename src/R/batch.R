#' @title Read a batch dataframe
#' @description
#' Read a batch file containing results for multiple chemical instances returning
#' a list of dataframes: 1 dataframe per chemical instance. The log2 fold change
#' data will also be converted to log10 fold change for ease of subsequent analysis.
#' @param batch_df Dataframe containing batch results. The first column should 
#' be the gene symbol followed zero or more metadata columns followed by data 
#' columns where each column represents log2 fold change data for a chemical 
#' instance. Each of these data columns should have a column header that describes
#' the chemical of the format *ChemicalName_Concentration_ConcentrationUnits*.
#' @param data_start_column Index of the column where the log2 fold change data
#' starts. Defaults to 2.
read_batch_df = function(batch_df, data_start_column = 2) {
  
  # Uppercase the gene symbols
  batch_df[, 1] <- toupper(batch_df[, 1]);
  
  #' The incoming dataframe may have metadata columns between the gene symbol column
  #' and the fold change columns. Build a sequence of columns of just the FC columns.
  data_cols_seq = seq(data_start_column, ncol(batch_df))
  n_data_cols = length(data_cols_seq)
  
  log_msg(glue::glue("Found {n_data_cols} chemical instances starting at column: {data_start_column}"))
  
  # Build a list to hold the results of splitting out the batched dataframe
  results = vector(mode = "list", length = n_data_cols)
  names(results) = colnames(batch_df)[data_cols_seq]
  result_i = 1
  
  for (i in data_cols_seq) {
    chemical_instance = colnames(batch_df)[i]
    log_msg(glue::glue("Found {chemical_instance}"))
    
    # Build a 2 column dataframe with the gene symbol column and the data for this instance
    chemical_instance_df = batch_df[, c(1, i)]
    colnames(chemical_instance_df) = c("Probe", chemical_instance)
    
    # Convert from log2 to log10
    chemical_instance_df[, 2] = log10(2 ^ chemical_instance_df[, 2])
    
    # Save result and increment result index
    results[[result_i]] = chemical_instance_df
    result_i = result_i +  1
  }
  
  return(results)
}

#' @title Read a batch TSV file.
#' @description Read a batch file into a dataframe.
#' @param batch_filepath Path to the tab separated input file.
#' @inheritParams read_batch_df
read_batch_file = function(batch_filepath, data_start_column = 2) {
  log_msg(glue::glue("Reading input data file: {batch_filepath}"))
  
  batch_data = read.table(
    file = batch_filepath, 
    header = TRUE, 
    sep = "\t", 
    quote = "", 
    fill = TRUE, 
    comment.char = ""
  )
  
  read_batch_df(batch_data, data_start_column)
}


#' @title Classify a batch of chemical results from a list of dataframes.
#' @param df_list List of dataframes to process.
#' @inheritParams check_classifier_version
#' @inheritParams read_batch_df
#' @inheritParams classify_general
#' @export
classify_batch_df = function(
    df_list,
    out_dir = NULL,
    data_start_column = 2,
    classifier_version = 2,
    plot_title = NULL,
    save_outputs = TRUE
) {
  if (is.null(plot_title)) {
    plot_title = "All Results"
  }
  
  n_result = length(df_list)
  result = data.frame(
    test_material = character(n_result),
    prob_genotoxic = numeric(n_result),
    prob_non_genotoxic = numeric(n_result),
    classification = character(n_result)
  )
  
  for (i in seq_along(df_list)) {
    df = df_list[[i]]
    df_name = names(df_list)[i]
    class_result = classify_general(
      input_data = df,
      classifier_version = classifier_version,
      out_dir = out_dir, 
      outfile_prefix = df_name, 
      plot_title = plot_title,
      save_outputs = save_outputs
    )
    result[i, "test_material"] = df_name
    result[i, "prob_genotoxic"] = class_result[["prob_genotoxic"]]
    result[i, "prob_non_genotoxic"] = class_result[["prob_non_genotoxic"]]
    result[i, "classification"] = class_result[["classification"]]
  }
  return(result)
}


#' @title Classify a batch of chemical results from a file.
#' @param input_file 
#' Path to a file containing a batch of input data to classify as tab delimited
#' text.
#' @inheritParams classify_batch_df
#' @export
classify_batch_file = function(
    batch_filepath,
    out_dir = NULL,
    data_start_column = 2,
    classifier_version = 2,
    plot_title = NULL
) {
  df_list = read_batch_file(batch_filepath, data_start_column)
  classify_batch_df(
    df_list = df_list, 
    out_dir = out_dir, 
    data_start_column = data_start_column,
    classifier_version = classifier_version
  )
}