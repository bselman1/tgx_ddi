#' @title Validate a chemical identifier is in the expected format.
#' @description Verify that the chemical identifier is in the expected format:
#' *ChemicalName_Concentration_Units*.
#' @param chemical_id The chemical identifier to validate.
#' @return A list containing the name, concentration, and units of the chemical.
validate_chemical_id = function(chemical_id) {
  if (!is.character(chemical_id)) {
    stop("Chemical ID must be a string.")
  }
  if (nchar(chemical_id) == 0) {
    stop("Chemical ID must not be empty.")
  }

  split_vals = unlist(strsplit(chemical_id, "_"))
  if (length(split_vals) != 3) {
    stop(glue::glue("Chemical id {chemical_id} does not have the expected format: ChemicalName_Concentration_Units."))
  }

  return(list(
    name = split_vals[1],
    concentration = as.numeric(split_vals[2]),
    units = split_vals[3]
  ))
}

set_valid_chemical_id_attribute = function(x, valid_chem_id) {
  attr(x, "valid_chem_id") = valid_chem_id
  return(x)
}

get_valid_chemical_id_attribute = function(x) {
  return(attr(x, "valid_chem_id"))
}


#' @title Validate a batch dataframe is in the expected format.
#' @description 
#' Validate the incoming batch data frame by:
#' 1. There are at least 2 columns in the dataframe.
#' 2. The first column is the gene symbol.
#' 3. Starting at the provided *data_start_column* (default=2), assume
#'    the remaining columns are log2 fold change data for each test chemical
#'    instance. These columns should have a column header that describe that
#'    uniquely labels the chemical instance of the format 
#'    *ChemicalName_Concentration_ConcentrationUnits*.
#' @param batch_df Dataframe containing batch results.
#' @param data_start_column Index of the column where the log2 fold change data
#' starts. Defaults to 2.
#' @return A validated dataframe with the following form:
#' 1. The gene symbol column is converted to uppercase and named "Probe".
#' 2. The remaining columns with the test chemical log2 fold change data converted
#'    to log10 fold change. The columns are sorted by chemical name then concentration.
validate_batch_df = function(batch_df, data_start_column = 2) {
  if (ncol(batch_df) < 2) {
    stop("Batch dataframe must have at least 2 columns.")
  }
  if (data_start_column < 2) {
    stop("data_start_column must be >= 2.")
  }

  # Uppercase the gene symbols
  result = batch_df
  colnames(result)[1] = "Probe"
  result[["Probe"]] <- toupper(result[["Probe"]])
  
  #' The incoming dataframe may have metadata columns between the gene symbol column
  #' and the fold change columns. Drop these columns.
  data_col_indices = seq(data_start_column, ncol(result))
  n_data_cols = length(data_col_indices)
  result = result[, c(1, data_col_indices)]
  log_msg(glue::glue("Found {n_data_cols} chemical instances starting at column: {data_start_column}"))
  
  # Transform the log2 fold change data to log10 and validate column headers
  for (i in 2:ncol(result)) {
    colname = colnames(result)[i]

    # Check the format of the column header
    valid_chem_id = try(validate_chemical_id(colname))
    if (is.error(valid_chem_id)) {
      error_msg = get_error_msg(valid_chem_id)
      stop(glue::glue("Invalid column header {colname}. Details: {error_msg}."))
    }
    
    col_values = result[[i]]
    
    # Make sure the column contains numeric values
    if (!is.numeric(col_values)) {  
      stop(glue::glue("Column {colname} has non-numeric data."))
    }

    # Convert from log2 to log10
    col_values = log10(2 ^ col_values)

    # Add an attribute to the column with the validated chemical id for easier sorting later
    col_values = set_valid_chemical_id_attribute(col_values, valid_chem_id)
    result[[i]] = col_values
  }

  return(result)
}


#' @title Sort a batch dataframe by chemical name then concentration.
#' @param valid_batch_df Validated batch dataframe as produced by [validate_batch_df].
sort_batch_df = function(valid_batch_df) {
  # We'll keep track of the unique list of chemical names in a list associated to a number of concentrations tested
  chem_conc_map = list()

  # Capture the unique list of chemical names and their tested concentrations
  for (i in 2:ncol(valid_batch_df)) {
    colname = colnames(valid_batch_df)[i]
    col_values = valid_batch_df[[i]]
    valid_chem_id = get_valid_chemical_id_attribute(col_values)

    if (is.null(valid_chem_id)) {
      stop(glue::glue("Column {colname} does not have a valid chemical id attribute."))
    }

    # Check if this is a new chemical name or an additional concentration to a previous chemical
    if (valid_chem_id$name %in% names(chem_conc_map)) {
      # Add the concentration to the list of concentrations for this chemical
      entry = chem_conc_map[[valid_chem_id$name]]
      entry$concentrations = c(entry$concentrations, valid_chem_id$concentration)
      entry$units = c(entry$units, valid_chem_id$units)
      chem_conc_map[[valid_chem_id$name]] = entry
    } else {
      # Add a new entry for this chemical
      chem_conc_map[[valid_chem_id$name]] = list(
        concentrations = valid_chem_id$concentration,
        units = valid_chem_id$units
      )
    }
  }

  # Sort by chemical name then concentration
  sorted_names = vector(mode="character", length=0)
  for (i in order(names(chem_conc_map))) {
    chemical_name = names(chem_conc_map)[i]
    entry = chem_conc_map[[i]]

    # Sort concentrations from low to high
    concentration_order = order(entry$concentrations)
    entry$concentrations = entry$concentrations[concentration_order]
    entry$units = entry$units[concentration_order]

    # Rebuild the chem id as ChemicalName_Concentration_Units and add to the sorted list
    for (conc_i in seq_along(entry$concentrations)) {
      chem_id = glue::glue("{chemical_name}_{entry$concentrations[conc_i]}_{entry$units[conc_i]}")
      sorted_names = c(sorted_names, chem_id)
    }
  }

    # Resort the result data frame to match the sorted chemical names
    result = valid_batch_df[, c("Probe", sorted_names)]
    return(result)
}


#' @title Split a batch dataframe into a list of dataframes.
#' @description
#' Read a batch data frame containing results for multiple chemical instances returning
#' a list of dataframes: 1 dataframe per chemical instance. The log2 fold change
#' data will also be converted to log10 fold change for ease of subsequent analysis.
#' @param valid_batch_df Validated batch dataframe as produced by [validate_batch_df].
#' @return A list of dataframes, one for each chemical instance.
split_validated_batch_df = function(valid_batch_df) {
  # Build a list to hold the results of splitting out the batched dataframe
  data_cols_seq = seq(2, ncol(valid_batch_df))
  results = vector(mode = "list", length = length(data_cols_seq))
  names(results) = colnames(valid_batch_df)[data_cols_seq]
  result_i = 1
  
  for (i in data_cols_seq) {
    chemical_instance = colnames(valid_batch_df)[i]
    
    # Build a 2 column dataframe with the gene symbol column and the data for this instance
    chemical_instance_df = valid_batch_df[, c(1, i)]
    colnames(chemical_instance_df) = c("Probe", chemical_instance)
    
    # Save result and increment result index
    results[[result_i]] = chemical_instance_df
    result_i = result_i +  1
  }
  
  return(results)
}


#' @title Read a batch TSV file.
#' @description Read a batch file into a dataframe.
#' @param batch_filepath Path to the tab separated input file.
read_batch_file = function(batch_filepath) {
  log_msg(glue::glue("Reading batch input data file: {batch_filepath}"))
  
  batch_data = read.table(
    file = batch_filepath, 
    header = TRUE, 
    sep = "\t", 
    quote = "", 
    fill = TRUE, 
    comment.char = ""
  )
}


#' @title Classify a batch of test chemical results from a file.
#' @description Classif the provided test materials from a batch file but do not
#' compare against the training data.
#' @param batch_filepath 
#' Path to a file containing a batch of input data to classify as tab delimited
#' text.
#' @param meta_info_df A dataframe containing metadata information about the test
#' materials in the batch file. Defaults to NULL. If provided, the data frame should
#  have the following columns:
#' * chem_id (Required): The chemical identifier in the format ChemicalName_Concentration_Units.
#' * actual_class (Optional): The actual genotoxicity class of the test material. One of
#'   "Genotoxic", "Non-Genotoxic", or "Unclassified". If provided, a bar will be added to
#'   the heatmap indicating the actual class vs the predicted class.
#' * sort_order (Optional): A numeric value to use to sort the test materials. If not provided,
#'   the test materials will be sorted by chemical name then concentration.
#' @param save_outputs If TRUE, save the resulting figures and tables to the output directory.
#' @param plot_title The title to use in the classification figures. Defaults to NULL.
#' @param out_dir The output directory to place the resulting figures and tables.
#' Defaults to the current working directory
#' @param outfile_prefix A prefix to use in the file names of any saved output from
#' this function. Defaults to "tgx_ddi_results".
#' @inheritParams classify
#' @export
classify_batch_file = function(
  batch_filepath,
  data_start_column = 2,
  classifier_version = 2,
  meta_info_df = NULL,
  save_outputs = TRUE,
  plot_title = NULL,
  out_dir = NULL,
  outfile_prefix = NULL
) {
  # Read and validate the batch file
  batch_df = read_batch_file(batch_filepath)
  valid_batch_df = validate_batch_df(batch_df, data_start_column)

  # If the meta_info_df is provided, validate that the chem_id column is present and matches the column names in the batch data
  if (!is.null(meta_info_df)) {
    if (!"chem_id" %in% colnames(meta_info_df)) {
      stop("meta_info_df must contain a column named 'chem_id'.")
    }
    if (!all(meta_info_df$chem_id %in% colnames(valid_batch_df))) {
      stop("All values in chem_id must match column names in the batch data.")
    }
  } else {
    log_msg("No meta data provided, using defaults.")
    meta_info_df = data.frame(chem_id = tail(colnames(valid_batch_df), -1))
  }

  # Validate actual classes if provided
  has_actual_classes = "actual_class" %in% colnames(meta_info_df)
  if (has_actual_classes) {
    valid_classes = meta_info_df$actual_class %in% c("Genotoxic", "Non-Genotoxic", "Unclassified")
    if (!all(valid_classes)) {
      stop("All values in actual_class must be one of 'Genotoxic', 'Non-Genotoxic', or 'Unclassified'.")
    }
    log_msg("Using actual classes provided from meta data.")
  }

  # Sort the batch data
  sort_order_override = "sort_order" %in% colnames(meta_info_df)
  if (sort_order_override) {
    # Check the names in the sort order override are valid
    log_msg("Using custom sorting from provided meta data.")
    ordered_ids = meta_info_df[order(meta_info_df$sort_order), "chem_id"]
    valid_batch_df = valid_batch_df[, c("Probe", ordered_ids)]
  } else {
    log_msg("Using default sorting: by chemical name then concentration.")
    valid_batch_df = sort_batch_df(valid_batch_df)
    # Add the default sort order to the meta data
    ordered_ids = tail(colnames(valid_batch_df), -1)
    default_sort_order = data.frame(
      chem_id = ordered_ids,
      sort_order = seq_along(ordered_ids)
    )
    meta_info_df = merge(meta_info_df, default_sort_order, by = "chem_id", sort = FALSE)
  }

  # Classify the batch data
  result = classify(
    input_data = valid_batch_df,
    classifier_version = classifier_version,
    include_training = FALSE
  )

  if (save_outputs) {
    # Default to the current working directory if out_dir is not provided
    if (is.null(out_dir)) {
      out_dir = getwd()
    }
    # Make the output directory if it doesn't exist
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    if (is.null(outfile_prefix)) {
      # Use "All" as the default prefix
      outfile_prefix = "tgx_ddi_results"
    }
    # Make sure the prefix doesn't contain any invalid file name characters
    outfile_prefix = fs::path_sanitize(outfile_prefix)

    # Add the meta data to the classification data frame
    result$classification_df = merge(result$classification_df, meta_info_df, by = "chem_id")

    # Create the outputs
    timestamp <- format(Sys.time(), "%d-%b-%Y")
    save_foldchange_file(out_dir, outfile_prefix, result$fc_matrix)
    save_gene_distance_file(out_dir, outfile_prefix, result$gene_dist)
    save_heatmap_pdf(
      fc_matrix = result$fc_matrix,
      classification_df = result$classification_df,
      gene_clust_dendo = result$gene_clust_dendo, 
      plot_title = plot_title, 
      timestamp = timestamp,
      out_dir = out_dir, 
      outfile_prefix = outfile_prefix
    )
    save_cluster_pdf(
      fc_matrix = result$fc_matrix, 
      plot_title = plot_title, 
      timestamp = timestamp, 
      out_dir = out_dir, 
      outfile_prefix = outfile_prefix
    )
    save_heatmap_png(
      fc_matrix = result$fc_matrix,
      classification_df = result$classification_df,
      gene_clust_dendo = result$gene_clust_dendo, 
      plot_title = plot_title, 
      timestamp = timestamp,
      out_dir = out_dir,
      outfile_prefix = outfile_prefix
    )
  }
  return(result)
}


#' @title Classify a batch of chemical results from a file against the training data.
#' @description Classify the provided test materials from a batch file and compare 
#' against the training data.
#' @inheritParams classify_batch_file
#' @export
classify_batch_file_against_training = function(
  batch_filepath,
  out_dir = NULL,
  data_start_column = 2,
  classifier_version = 2,
  plot_title = NULL,
  save_outputs = TRUE
) {
  # Read, validate and split the batch file
  batch_df = read_batch_file(batch_filepath)
  valid_batch_df = validate_batch_df(batch_df, data_start_column)
  df_list = split_validated_batch_df(valid_batch_df)

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
    
    class_result = classify(
      input_data = df,
      classifier_version = classifier_version,
      include_training = TRUE
    )
    
    # Test material will be the last result in the classification data frame
    test_material_result = tail(class_result$classification_df, 1)
    
    result[i, "test_material"] = df_name
    result[i, "prob_genotoxic"] = test_material_result[["prob_genotoxic"]]
    result[i, "prob_non_genotoxic"] = test_material_result[["prob_non_genotoxic"]]
    result[i, "classification"] = test_material_result[["classification"]]
    
    if (save_outputs) {
      # Default to the current working directory if out_dir is not provided
      if (is.null(out_dir)) {
        out_dir = getwd()
      }
      # Make the output directory if it doesn't exist
      if (!dir.exists(out_dir)) {
        dir.create(out_dir)
      }
      # Make sure the prefix doesn't contain any invalid file name characters
      outfile_prefix = fs::path_sanitize(df_name)

      # Create the outputs
      timestamp <- format(Sys.time(), "%d-%b-%Y")
      save_foldchange_file(out_dir, outfile_prefix, class_result$fc_matrix)
      save_gene_distance_file(out_dir, outfile_prefix, class_result$gene_dist)
      save_heatmap_with_cluster_training_pdf(
        fc_matrix = class_result$fc_matrix,
        classification_df = class_result$classification_df,
        gene_clust_dendo = class_result$gene_clust_dendo, 
        plot_title = plot_title, 
        timestamp = timestamp,
        out_dir = out_dir, 
        outfile_prefix = outfile_prefix
      )
    }
  }
  return(result)
}