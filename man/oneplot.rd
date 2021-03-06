\name{oneplot}
\alias{oneplot}
\alias{oneplot,OutlierDM-method}
\alias{oneplot,OutlierDM,numeric-method}
\title{
	Draw a dot-plot for a selected observation (peptide)
}
\description{
  This function draws a dot plot for a selected peptide based on the \code{OutlierDM} object.
}
\usage{
oneplot(object, i, ...)

\S4method{oneplot}{OutlierDM}(object, i = 1, pick = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted object}
  \item{i}{ a row number in order or drawing a dot-plot}
  \item{pick}{the number of locators to denote the names of the peptide}
  \item{...}{do not use at this term}
}
\seealso{
\code{\link{odm}}
}
\examples{
  data(lcms3)
  fit = odm(lcms3, method = "grubbs")
  oneplot(fit, i = 100)

  \dontrun{
  # Add row name
  oneplot(fit, i = 100, pick = 1)
  }
}
\keyword{subfunction}
