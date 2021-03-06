% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNAclinicData-accessors.R
\docType{methods}
\name{chromosomes,CNAclinicData-method}
\alias{chromosomes,CNAclinicData-method}
\title{An S4 class to contain the copy number, segmentation and calls from
multiple samples.}
\usage{
\S4method{chromosomes}{CNAclinicData}(object)
}
\description{
An object of this class contains several data.frame objects.
 \code{bins} contains information on binned genomic windows such as the
 chromosome, start and end of the bins as well as a condition to use or
 filter the bins. All other data.frames have rows corresponding to genomic
 bins and the columns to separate samples.
 Of these, \code{copyNumber} contains the bias corrected and normalized
 count data (log2 ratios), the others named \code{segCBS}, \code{segHMM},
 \code{segPLS} and \code{segLACBS} contain the value of the segment for
 each genomic bin calculated from multiple segmentation algorithms.
 The \code{segSummary} data.frame holds the summarised consensus segment
 values for each bin and input sample, while \code{calls} holds the
 respective copy number calls.
}
