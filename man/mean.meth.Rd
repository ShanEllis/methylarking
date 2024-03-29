% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mean.meth.R
\name{mean.meth}
\alias{mean.meth}
\title{Mean Methylation Difference -- by case-control}
\usage{
\method{mean}{meth}(X, treatment)
}
\arguments{
\item{treatment}{vector with case control information \code{treatment}}

\item{methyl.mat}{methylation matrix \code{methyl.mat}}
}
\value{
mean methylation differences between cases and controls
}
\description{
Mean Methylation Difference -- by case-control
}
\examples{
#' mean.meth.diff=apply(methyl.obj,1,mean.meth)
}
\keyword{difference}
\keyword{mean}
\keyword{methyl.diff,}
\keyword{methylBase,}
\keyword{methylKit,}
\keyword{methylation}

