% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QQplot.R
\name{QQplot}
\alias{QQplot}
\title{QQ plot of all sites from DMS Analysis}
\usage{
QQplot(DMS, filename, cols = "dodgerblue")
}
\arguments{
\item{DMS}{output dataframe produced by runlimma \code{DMS}}

\item{filename}{filenames to attach \code{filename}}

\item{cols}{color of points \code{cols}}
}
\value{
QQplot single site analysis output
}
\description{
QQ plot of all sites from DMS Analysis
}
\examples{
QQplot(DMS,"test",cols="red")
}
\keyword{QQplot}
\keyword{limma,}

