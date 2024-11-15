% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\name{calculate_score}
\alias{calculate_score}
\title{Calculate genotox score}
\usage{
calculate_score(probe_values, probe_scores, probe_score_sds)
}
\arguments{
\item{probe_values}{Vector of Log10 fold-change values for each probe for the
test chemical.}

\item{probe_scores}{Vector of training scores for each probe (either genotoxic
or non-genotoxic) provided by the classifier data frame.}

\item{probe_score_sds}{Vector of training standard deviation values for each
probe provided by the classifier data frame.}
}
\description{
Calculate genotox score
}
