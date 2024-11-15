% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\name{predict_classification}
\alias{predict_classification}
\title{Predict genotox probabilites and classifications for an input data frame.}
\usage{
predict_classification(valid_df)
}
\arguments{
\item{valid_df}{A list as returned by \link{validate_merged_input_df} describing
the validated input data}
}
\value{
A data frame with a row for each tested chemical and the following columns:
\itemize{
\item prob_genotoxic
\item prob_non_genotoxic
\item classification
}
}
\description{
Predict genotox probabilites and classifications for an input data frame.
}
