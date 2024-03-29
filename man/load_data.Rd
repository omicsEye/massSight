% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{load_data}
\alias{load_data}
\title{Load Data}
\usage{
load_data(input, type = "known", sheet = 1, id = "Metabolite")
}
\arguments{
\item{input}{A string of the file path or a data frame of the data to be
loaded.}

\item{type}{A string indicating the type of data to be loaded. This can be
"known" or "unknown".}

\item{sheet}{A string or integer indicating the sheet number of the excel
file to be loaded.}

\item{id}{A string indicating the column name of the compound ID.}
}
\value{
A list of the sample metadata, feature metadata, and the data matrix.
}
\description{
This function will load LC-MS data from a file.
}
