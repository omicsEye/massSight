% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ml_match.R
\name{ml_match}
\alias{ml_match}
\title{ML Match}
\usage{
ml_match(ms1, ms2, mz_thresh = 15, rt_thresh = 1, prob_thresh = 0.5, seed = 72)
}
\arguments{
\item{ms1}{A \code{massSight} object representing the results of a preprocessed
LC-MS experiment}

\item{ms2}{A \code{massSight} object representing the results of a second
preprcoessed LC-MS experiment}

\item{mz_thresh}{\code{numeric} Mass to Charge threshold. Used to limit potential matches
between metabolites.}

\item{rt_thresh}{\code{numeric} Retention Time threshold. Used to limit potential matches
between metabolites}

\item{seed}{Seed value for reproducibility}
}
\value{
A dataframe consisting of predicted metabolite pairs
}
\description{
Trains, fits, and uses a machine learning model on known
metabolite data to predict unknown metabolite pair matches.
}
\note{
This function requires semi-annotated data (some metabolites must be
named)
}
\examples{
\dontrun{
ml_match(ms1, ms2, mz_thresh = 15, rt_thresh = .5, seed = 2)
}
}
