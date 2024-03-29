% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{repetitions}
\alias{repetitions}
\title{Repetitions}
\usage{
repetitions(s, full = F)
}
\arguments{
\item{s}{Sequence vector.}

\item{full}{Optional. Boolean.}
}
\value{
If full returns a vector of length = length(s) with whether the ith item is a repetition of the last. Otherwise returns the mean repetition rate.
}
\description{
Repetitions
}
\examples{
s <- sample(1:10, 100, T)
repetitions(s)
repetitions(s, full=T)
}
