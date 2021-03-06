% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/callData.R
\docType{methods}
\name{callData,CNAclinicData-method}
\alias{callData,CNAclinicData-method}
\title{Call segments as 'gain', 'loss' or 'neutral'}
\usage{
\S4method{callData}{CNAclinicData}(object, callTypeLog2R = "summary",
  callThreshLog2R = c(-0.15, 0.15))
}
\arguments{
\item{object}{CNAclinicData.}
}
\value{
An object of class CNAclinicData
}
\description{
Call segments as 'gain', 'loss' or 'neutral'
}
\examples{
     \dontrun{
      vignette("CNAclinic")
     }

}
\author{
Dineika Chandrananda
}
