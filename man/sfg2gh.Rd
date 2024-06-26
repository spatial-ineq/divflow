% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph-and-ergm-fcns.R
\name{sfg2gh}
\alias{sfg2gh}
\title{sfg2gh}
\usage{
sfg2gh(sfg, directed = F, ...)
}
\arguments{
\item{sfg}{OD dataframe}

\item{directed}{whether to leave directed or make undirected (which aggregates
trips to/from the same tracts, rather than duplicating the O-D pair for each
direction.)}
}
\description{
From a origin-destination data.frame representing safegraph data, turn it into a
tidygraph object. Assumes \code{origin}, \code{dest}, and \code{n} columns. Also attaches
geometry and calculates "normalized connectedness," which will be different if
it's made undirected or left directed.
}
