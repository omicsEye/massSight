% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auto_combine.R
\name{auto_combine}
\alias{auto_combine}
\title{Auto Combine}
\usage{
auto_combine(
  ms1,
  ms2,
  rt_lower = -0.5,
  rt_upper = 0.5,
  mz_lower = -15,
  mz_upper = 15,
  minimum_intensity = 10,
  iso_method = "manual",
  eps = 0.1,
  rt_iso_threshold = 0.1,
  mz_iso_threshold = 5,
  match_method = "unsupervised",
  smooth_method = "gam",
  weights = c(1, 1, 1),
  log = NULL,
  output = NULL
)
}
\arguments{
\item{ms1}{A \code{massSight} object representing the results of a preprocessed
LC-MS experiment.}

\item{ms2}{A \code{massSight} object representing the results of a second
preprocessed LC-MS experiment.}

\item{rt_lower}{A numeric indicating the lower bound of the RT
range to be considered for aligning two metabolites.}

\item{rt_upper}{A numeric indicating the upper bound of the RT
range to be considered for aligning two metabolites.}

\item{mz_lower}{A numeric indicating the lower bound of the m/z
range to be considered for aligning two metabolites.}

\item{mz_upper}{A numeric indicating the upper bound of the m/z
range to be considered for aligning two metabolites.}

\item{minimum_intensity}{A numeric indicating the minimum intensity
to be considered for alignment.}

\item{iso_method}{The isolation method used before modeling drift. Can
either be "manual" or "dbscan".}

\item{eps}{Epsilon value for dbscan algorithm. Only used if iso_method =
"dbscan"}

\item{rt_iso_threshold}{A numeric indicating the isolation
parameter for RT.}

\item{mz_iso_threshold}{A numeric indicating the isolation
parameter for m/z.}

\item{match_method}{A character indicating the initial matching method to
be used to detect inter-batch variability. Options are "unsupervised" and
"supervised".}

\item{smooth_method}{A character indicating the smoothing method to
be used. Options are "lowess", "spline", and "gaussian".}

\item{weights}{A numeric vector indicating the weights to be used for
the alignment.}

\item{log}{A character indicating the name of the log file.}

\item{output}{A character indicating the directory to save the output. If NULL,
the output will be saved in the current working directory.}
}
\value{
A \code{MergedMSObject} containing the combined data.
}
\description{
Combines two \code{massSight} objects, resulting in a single
\code{MergedMSObject}.
}
