% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processForSegmentation.R
\name{processForSegmentation}
\alias{processForSegmentation}
\title{Process reads counts from BAM files to prepare input
for segmentation algorithms}
\usage{
processForSegmentation(bamfiles = NULL, bamnames = NULL,
  refSamples = NULL, pathToBams = NULL, ext = "bam", binSize = NULL,
  genome = "hg19", outputType = "CNAclinicData",
  typeOfPreMadeBins = "SR50", userMadeBins = NULL,
  cache = getOption("QDNAseq::cache", FALSE), minMapq = 20,
  pairedEnds = NULL, isPaired = NA, isProperPair = NA,
  isUnmappedQuery = FALSE, hasUnmappedMate = NA, isMinusStrand = NA,
  isMateMinusStrand = NA, isFirstMateRead = NA, isSecondMateRead = NA,
  isSecondaryAlignment = NA, isDuplicate = FALSE, residualFilter = TRUE,
  blacklistFilter = TRUE, mappabilityFilter = 15,
  chromosomesFilter = c("X", "Y", "M", "MT"), spanForLoess = 0.65,
  familyForLoess = "symmetric", maxIterForCorrection = 1,
  cutoffForCorrection = 4, variablesForCorrection = c("gc", "mappability"),
  methodOfCorrection = "ratio", methodOfNormalization = "median",
  logTransformForSmoothing = TRUE, skipMedianNormalization = FALSE,
  skipOutlierSmoothing = FALSE, saveCountData = FALSE,
  filename = "corrected_QDNAseqCopyNumbers")
}
\arguments{
\item{bamfiles}{A \code{\link[base]{character}}  vector of BAM file names
with or without full path. If NULL (default), all files with extension
.bam, are read from directory path.}

\item{bamnames}{An optional \code{\link[base]{character}} vector of sample
names. Defaults to file names with extension \code{.bam} removed.
\code{bamnames} must be provided if \code{refSamples} is not NULL.}

\item{refSamples}{An optional \code{\link[base]{character}}  vector of the
reference sample names that are to be used in normalizing each sample in
\code{bamnames}.
If not NULL (default), \code{refSamples} must be the same length as
\code{bamnames} and should only include sample names contained in
\code{bamnames}. See vignette for further details.}

\item{pathToBams}{If \code{bamfiles} is NULL, all files ending with ".bam"
extension will be read from this path.
If NULL, defaults to the current working directory.}

\item{ext}{Input files extension. Defaults to "bam".}

\item{binSize}{A \code{\link[base]{numeric}} scalar specifying the width of
the bins in units of kbp (1000 base pairs), e.g. \code{binSize=50}
corresponds to 50 kbp bins.}

\item{genome}{Genome build used to align sequencing reads.
Currently, CNAclinic only allows \code{"hg19"} (default).
Also see: \code{userMadeBins}}

\item{outputType}{Return an object of class \code{"QDNAseqCopyNumbers"} or
\code{"CNAclinicData"} (default).}

\item{typeOfPreMadeBins}{A \code{\link[base]{character}} string to specify
the read type (single/paired) and length used to generate pre-made annotation.
e.g \code{"SR50"} (default) or \code{"PE100"}.}

\item{userMadeBins}{An optional data.frame or
an \code{\link[Biobase]{AnnotatedDataFrame}} object containing bin
annotations created using the \code{\link[QDNAseq]{createBins}} function.
Consult the \pkg{QDNAseq} vignette for further information.}

\item{cache}{Whether to read and write intermediate cache files, which
speeds up subsequent analyses of the same files. Requires packages
R.cache and digest (both available on CRAN) to be installed. Defaults
to getOption("QDNAseq::cache", FALSE)}

\item{minMapq}{If quality scores exists, the minimum quality score required
in order to keep a read (20, default).}

\item{pairedEnds}{A boolean value or vector specifying whether the BAM files
contain paired-end data or not.}

\item{isPaired}{A \code{\link[base]{logical}}(1) indicating whether unpaired
(FALSE), paired (TRUE), or any (NA, default) read should be returned.}

\item{isProperPair}{A \code{\link[base]{logical}}(1) indicating whether
improperly paired (FALSE), properly paired (TRUE), or any (NA, default) read
should be returned. A properly paired read is defined by the alignment algorithm
and might, e.g., represent reads aligning to identical reference sequences and
with a specified distance.}

\item{isUnmappedQuery}{A \code{\link[base]{logical}}(1) indicating whether
unmapped (TRUE), mapped (FALSE, default), or any (NA) read should be returned.}

\item{hasUnmappedMate}{A \code{\link[base]{logical}}(1) indicating whether
reads with mapped (FALSE), unmapped (TRUE), or any (NA, default) mate should
be returned.}

\item{isMinusStrand}{A \code{\link[base]{logical}}(1) indicating whether reads
aligned to the plus (FALSE), minus (TRUE), or any (NA, default) strand should
be returned.}

\item{isMateMinusStrand}{A \code{\link[base]{logical}}(1) indicating whether
mate reads aligned to the plus (FALSE), minus (TRUE), or any (NA, default)
strand should be returned.}

\item{isFirstMateRead}{A \code{\link[base]{logical}}(1) indicating whether the
first mate read should be returned (TRUE) or not (FALSE), or whether mate read
number should be ignored (NA, default).}

\item{isSecondMateRead}{A \code{\link[base]{logical}}(1) indicating whether the
 second mate read should be returned (TRUE) or not (FALSE), or whether mate
read number should be ignored (NA, default).}

\item{isSecondaryAlignment}{A \code{\link[base]{logical}}(1) indicating whether
 alignments that are primary (FALSE), are not primary (TRUE) or whose primary
status does not matter (NA, default) should be returned.}

\item{isDuplicate}{A \code{\link[base]{logical}}(1) indicating that
un-duplicated (FALSE, default), duplicated (TRUE), or any (NA) reads should be
 returned.}

\item{residualFilter}{Either a \code{\link[base]{logical}} specifying whether
to filter based on loess residuals of the calibration set or if a numeric, the
number of standard deviations to use as the cutoff. Default is TRUE, which
corresponds to 4.0 standard deviations.}

\item{blacklistFilter}{Either a \code{\link[base]{logical}} specifying whether
to filter based on overlap with ENCODE blacklisted regions, or if numeric,
the maximum percentage of overlap allowed. Default is @TRUE, which corresponds
 to no overlap allowed (i.e. value of 0).}

\item{mappabilityFilter}{A \code{\link[base]{numeric}} in \eqn{[0,100]} to
specify filtering out bins with mappabilities lower than the number specified
(15, default). FALSE will not filter based on mappability.}

\item{chromosomesFilter}{A \code{\link[base]{character}} vector specifying
which chromosomes to filter out. Defaults to the sex chromosomes and
mitochondrial reads, i.e. \code{c("X", "Y", "M", "MT")}. Use NA to use all
chromosomes.}

\item{spanForLoess}{For @see "stats::loess", the parameter alpha which
controls the degree of smoothing.}

\item{familyForLoess}{For @see "stats::loess", if "gaussian" fitting is by
least-squares, and if "symmetric" a re-descending M estimator is
used with Tukey's biweight function.}

\item{maxIterForCorrection}{An integer(1) specifying the maximum number of
iterations to perform, default is 1. If larger, after the first loess fit, bins
 with median residuals larger than \code{cutoffForCorrection} are removed, and the
 fitting repeated until the list of bins to use stabilizes or after
 \code{maxIter} iterations.}

\item{cutoffForCorrection}{A numeric(1) specifying the number of standard
deviations (as estimated with @see "matrixStats::madDiff") the cutoff for
removal of bins with median residuals larger than the cutoff. Not used if
 \code{maxIter=1} (default).}

\item{variablesForCorrection}{A character vector specifying which variables
to include in the correction. Can be c("gc", "mappability") (the default) or
"gc", or "mappability".}

\item{methodOfCorrection}{A \code{\link[base]{character}} string speficying
the correction method. \code{ratio} (default) divides \code{counts} with
\code{fit}. \code{median} calculates the median \code{fit}, and defines the
 correction for bins with GC content \code{gc} and mappability
 \code{map} as \code{median(fit) - fit(gc,map)}, which is added to
 \code{counts}. Method \code{none} leaves \code{counts} untouched.}

\item{methodOfNormalization}{A \code{\link[base]{character}} string specifying
the normalization method. Choices are "mean", "median" (default), or "mode".}

\item{logTransformForSmoothing}{If TRUE (default), data will be
log2-transformed for smoothing.}

\item{skipMedianNormalization}{Skip this step if TRUE.
Recommended when normalizing by refSamples}

\item{skipOutlierSmoothing}{Skip this specific step if TRUE.}

\item{saveCountData}{Save an object of class \linkS4class{QDNAseqCopyNumbers}
after the GC/mappability correction step. default is FALSE}

\item{filename}{Filename to save the before mentioned object.}
}
\value{
Returns an object of class \linkS4class{CNAclinicData} (default) or
\linkS4class{QDNAseqCopyNumbers}
}
\description{
\code{processForSegmentation} is a wrapper function that reads in BAM
files and carries out binning, filtering, bias correcting, smoothing and
normalizing of the read counts using functions of the \pkg{QDNAseq}
package.
}
\examples{
     \dontrun{
      vignette("CNAclinic")
     }

}
\seealso{
Internally, the following functions of the \pkg{QDNAseq}
package are used: \code{\link[QDNAseq]{getBinAnnotations}},
\code{\link[QDNAseq]{binReadCounts}}, \code{\link{applyFilters}},
\code{\link[QDNAseq]{estimateCorrection}}, \code{\link[QDNAseq]{correctBins}},
\code{\link[QDNAseq]{normalizeBins}}, \code{\link[QDNAseq]{smoothOutlierBins}} and
\code{\link[QDNAseq]{compareToReference}}
}
\author{
Dineika Chandrananda
}
