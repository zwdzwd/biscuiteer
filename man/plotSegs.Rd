% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotSegs.R
\name{plotSegs}
\alias{plotSegs}
\title{plotting WGBScopy results (modified from HMMcopy)}
\usage{
plotSegs(cn, segs = NULL, chr = "chr22", name = "Segmented", ...)
}
\arguments{
\item{cn}{binned, corrected read depths from correctCoverage()}

\item{segs}{a segmentation result from WGBScopy()}

\item{chr}{what chromosome to display (chr22)}

\item{name}{a name for the plot ("Segmented")}

\item{...}{other arguments to pass on to plot}
}
\description{
plotting WGBScopy results (modified from HMMcopy)
}
