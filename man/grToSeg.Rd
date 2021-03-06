% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grToSeg.R
\name{grToSeg}
\alias{grToSeg}
\title{dump GRanges to .seg files}
\usage{
grToSeg(gr, filename = NULL, minAbs = NULL)
}
\arguments{
\item{gr}{A GRanges or GRangesList to dump as .seg}

\item{filename}{Where to save the result (unsaved if NULL, the default)}

\item{minAbs}{Optional minimum absolute gain/loss cutoff (default is NULL)}
}
\value{
\preformatted{     A data.frame(ID,chrom,loc.start,loc.end,num.mark,seg.mean)
}
}
\description{
dump GRanges to .seg files
}
