% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch.R
\name{classify_batch_df}
\alias{classify_batch_df}
\title{Classify a batch of chemical results from a list of dataframes.}
\usage{
classify_batch_df(
  df_list,
  out_dir = NULL,
  data_start_column = 2,
  classifier_version = 2,
  plot_title = NULL,
  save_outputs = TRUE
)
}
\arguments{
\item{df_list}{List of dataframes to process.}

\item{out_dir}{The output directory to place the resulting figures and tables.
Defaults to the current working directory}

\item{data_start_column}{Index of the column where the log2 fold change data
starts. Defaults to 2.
The incoming dataframe may have metadata columns between the gene symbol column
and the fold change columns. Build a sequence of columns of just the FC columns.}

\item{classifier_version}{Version of the trained DDI classifier to use. Possible values are 1 or 2.}
}
\description{
Classify a batch of chemical results from a list of dataframes.
}
