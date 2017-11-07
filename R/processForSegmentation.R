
#' Process reads counts from BAM files to prepare input
#' for segmentation algorithms
#'
#' \code{processForSegmentation} is a wrapper function that reads in BAM
#' files and carries out binning, filtering, bias correcting, smoothing and
#' normalizing of the read counts using functions of the \pkg{QDNAseq}
#' package.
#'
#' @param bamfiles A \code{\link[base]{character}}  vector of BAM file names
#'  with or without full path. If NULL (default), all files with extension
#'  .bam, are read from directory path.
#' @param bamnames An optional \code{\link[base]{character}} vector of sample
#' names. Defaults to file names with extension \code{.bam} removed.
#' \code{bamnames} must be provided if \code{refSamples} is not NULL.
#' @param refSamples An optional \code{\link[base]{character}}  vector of the
#' reference sample names that are to be used in normalizing each sample in
#' \code{bamnames}.
#' If not NULL (default), \code{refSamples} must be the same length as
#' \code{bamnames} and should only include sample names contained in
#' \code{bamnames}. See vignette for further details.
#' @param pathToBams If \code{bamfiles} is NULL, all files ending with ".bam"
#' extension will be read from this path.
#' If NULL, defaults to the current working directory.
#' @param ext Input files extension. Defaults to "bam".
#' @param binSize A \code{\link[base]{numeric}} scalar specifying the width of
#' the bins in units of kbp (1000 base pairs), e.g. \code{binSize=50}
#' corresponds to 50 kbp bins.
#' @param genome Genome build used to align sequencing reads.
#' Currently, CNAclinic only allows \code{"hg19"} (default).
#' Also see: \code{userMadeBins}
#' @param outputType Return an object of class \code{"QDNAseqCopyNumbers"} or
#' \code{"CNAclinicData"} (default).
#' @param typeOfPreMadeBins A \code{\link[base]{character}} string to specify
#' the read type (single/paired) and length used to generate pre-made annotation.
#' e.g \code{"SR50"} (default) or \code{"PE100"}.
#' @param userMadeBins An optional data.frame or
#' an \code{\link[Biobase]{AnnotatedDataFrame}} object containing bin
#' annotations created using the \code{\link[QDNAseq]{createBins}} function.
#' Consult the \pkg{QDNAseq} vignette for further information.
#' @param cache  Whether to read and write intermediate cache files, which
#' speeds up subsequent analyses of the same files. Requires packages
#' R.cache and digest (both available on CRAN) to be installed. Defaults
#' to getOption("QDNAseq::cache", FALSE)
#' @param minMapq If quality scores exists, the minimum quality score required
#' in order to keep a read (20, default).
#' @param pairedEnds A boolean value or vector specifying whether the BAM files
#' contain paired-end data or not.
#' @param isPaired A \code{\link[base]{logical}}(1) indicating whether unpaired
#' (FALSE), paired (TRUE), or any (NA, default) read should be returned.
#' @param isProperPair A \code{\link[base]{logical}}(1) indicating whether
#' improperly paired (FALSE), properly paired (TRUE), or any (NA, default) read
#' should be returned. A properly paired read is defined by the alignment algorithm
#' and might, e.g., represent reads aligning to identical reference sequences and
#' with a specified distance.
#' @param isUnmappedQuery A \code{\link[base]{logical}}(1) indicating whether
#' unmapped (TRUE), mapped (FALSE, default), or any (NA) read should be returned.
#' @param hasUnmappedMate A \code{\link[base]{logical}}(1) indicating whether
#' reads with mapped (FALSE), unmapped (TRUE), or any (NA, default) mate should
#' be returned.
#' @param isMinusStrand A \code{\link[base]{logical}}(1) indicating whether reads
#' aligned to the plus (FALSE), minus (TRUE), or any (NA, default) strand should
#' be returned.
#' @param isMateMinusStrand A \code{\link[base]{logical}}(1) indicating whether
#' mate reads aligned to the plus (FALSE), minus (TRUE), or any (NA, default)
#' strand should be returned.
#' @param isFirstMateRead A \code{\link[base]{logical}}(1) indicating whether the
#' first mate read should be returned (TRUE) or not (FALSE), or whether mate read
#' number should be ignored (NA, default).
#' @param isSecondMateRead A \code{\link[base]{logical}}(1) indicating whether the
#'  second mate read should be returned (TRUE) or not (FALSE), or whether mate
#' read number should be ignored (NA, default).
#' @param isSecondaryAlignment A \code{\link[base]{logical}}(1) indicating whether
#'  alignments that are primary (FALSE), are not primary (TRUE) or whose primary
#' status does not matter (NA, default) should be returned.
#' @param isDuplicate A \code{\link[base]{logical}}(1) indicating that
#' un-duplicated (FALSE, default), duplicated (TRUE), or any (NA) reads should be
#'  returned.
#' @param residualFilter Either a \code{\link[base]{logical}} specifying whether
#' to filter based on loess residuals of the calibration set or if a numeric, the
#' number of standard deviations to use as the cutoff. Default is TRUE, which
#' corresponds to 4.0 standard deviations.
#' @param blacklistFilter Either a \code{\link[base]{logical}} specifying whether
#' to filter based on overlap with ENCODE blacklisted regions, or if numeric,
#' the maximum percentage of overlap allowed. Default is @TRUE, which corresponds
#'  to no overlap allowed (i.e. value of 0).
#' @param mappabilityFilter A \code{\link[base]{numeric}} in \eqn{[0,100]} to
#' specify filtering out bins with mappabilities lower than the number specified
#' (15, default). FALSE will not filter based on mappability.
#' @param chromosomesFilter A \code{\link[base]{character}} vector specifying
#' which chromosomes to filter out. Defaults to the sex chromosomes and
#' mitochondrial reads, i.e. \code{c("X", "Y", "M", "MT")}. Use NA to use all
#' chromosomes.
#' @param spanForLoess For @see "stats::loess", the parameter alpha which
#' controls the degree of smoothing.
#' @param familyForLoess For @see "stats::loess", if "gaussian" fitting is by
#' least-squares, and if "symmetric" a re-descending M estimator is
#' used with Tukey's biweight function.
#' @param maxIterForCorrection An integer(1) specifying the maximum number of
#' iterations to perform, default is 1. If larger, after the first loess fit, bins
#'  with median residuals larger than \code{cutoffForCorrection} are removed, and the
#'  fitting repeated until the list of bins to use stabilizes or after
#'  \code{maxIter} iterations.
#' @param cutoffForCorrection A numeric(1) specifying the number of standard
#' deviations (as estimated with @see "matrixStats::madDiff") the cutoff for
#' removal of bins with median residuals larger than the cutoff. Not used if
#'  \code{maxIter=1} (default).
#' @param variablesForCorrection A character vector specifying which variables
#' to include in the correction. Can be c("gc", "mappability") (the default) or
#' "gc", or "mappability".
#' @param methodOfCorrection A \code{\link[base]{character}} string speficying
#' the correction method. \code{ratio} (default) divides \code{counts} with
#' \code{fit}. \code{median} calculates the median \code{fit}, and defines the
#'  correction for bins with GC content \code{gc} and mappability
#'  \code{map} as \code{median(fit) - fit(gc,map)}, which is added to
#'  \code{counts}. Method \code{none} leaves \code{counts} untouched.
#' @param methodOfNormalization A \code{\link[base]{character}} string specifying
#'  the normalization method. Choices are "mean", "median" (default), or "mode".
#' @param logTransformForSmoothing If TRUE (default), data will be
#' log2-transformed for smoothing.
#' @param skipMedianNormalization Skip this step if TRUE.
#' Recommended when normalizing by refSamples
#' @param skipOutlierSmoothing Skip this specific step if TRUE.
#' @param saveCountData Save an object of class \linkS4class{QDNAseqCopyNumbers}
#' after the GC/mappability correction step. default is FALSE
#' @param filename Filename to save the before mentioned object.
#'
#' @return Returns an object of class \linkS4class{CNAclinicData} (default) or
#' \linkS4class{QDNAseqCopyNumbers}
#'
#' @seealso Internally, the following functions of the \pkg{QDNAseq}
#' package are used: \code{\link[QDNAseq]{getBinAnnotations}},
#' \code{\link[QDNAseq]{binReadCounts}}, \code{\link{applyFilters}},
#' \code{\link[QDNAseq]{estimateCorrection}}, \code{\link[QDNAseq]{correctBins}},
#' \code{\link[QDNAseq]{normalizeBins}}, \code{\link[QDNAseq]{smoothOutlierBins}} and
#' \code{\link[QDNAseq]{compareToReference}}
#'
#' @importFrom QDNAseq getBinAnnotations binReadCounts applyFilters
#' estimateCorrection correctBins normalizeBins smoothOutlierBins
#' compareToReference
#'
#' @import QDNAseq.hg19
#' @importClassesFrom QDNAseq QDNAseqReadCounts QDNAseqCopyNumbers
#' @export processForSegmentation
#'
#' @author Dineika Chandrananda
#' @examples
#'      \dontrun{
#'       vignette("CNAclinic")
#'      }
#'


processForSegmentation = function(
    bamfiles=NULL, bamnames=NULL, refSamples=NULL, pathToBams=NULL,
    ext="bam",
    binSize=NULL, genome="hg19", outputType="CNAclinicData",
    typeOfPreMadeBins="SR50",
    userMadeBins=NULL,
    cache=getOption("QDNAseq::cache", FALSE),
    minMapq=20, pairedEnds=NULL, isPaired=NA, isProperPair=NA,
    isUnmappedQuery=FALSE, hasUnmappedMate=NA,
    isMinusStrand=NA, isMateMinusStrand=NA,
    isFirstMateRead=NA, isSecondMateRead=NA,
    isSecondaryAlignment=NA,
    isDuplicate=FALSE,
    residualFilter=TRUE, blacklistFilter=TRUE, mappabilityFilter=15,
    chromosomesFilter=c("X", "Y", "M", "MT"),
    spanForLoess=0.65, familyForLoess="symmetric",
    maxIterForCorrection=1, cutoffForCorrection=4.0,
    variablesForCorrection=c("gc", "mappability"), methodOfCorrection="ratio",
    methodOfNormalization="median",
    logTransformForSmoothing=TRUE,
    skipMedianNormalization=FALSE, skipOutlierSmoothing=FALSE,
    saveCountData=FALSE, filename="corrected_QDNAseqCopyNumbers"){

    ############################################################################
    # Check user specified arguments
    ############################################################################

    # if bamnames are present but not distinct complain!
    if(!is.null(bamnames)){
        if(!(length(unique(bamnames)) == length(bamfiles))){
            errorMsg <- "Please provide unique bamnames
            that match the order of the bamfiles"
            stop(errorMsg)
        }
    }

    if(!is.null(refSamples)){
        if(!all(refSamples %in% c(bamnames, NA, 'NA', 'drop'))
           || length(refSamples) != length(bamnames)){
            errorMsg <-
                "refSamples should contain reference sample names
            that are to be used in normalizing each sample in bamnames.
            If not NULL, refSamples must be the same length as bamnames
            and should only include sample names contained in bamnames, NA or
            'drop'.
            \n
            When 'drop', the corresponding sample from bamnames will be removed
            from the output. When NA, the corresponding sample will be
            normalized by its median.
            \n
            As an example, if bamnames = c('tumour1', 'tumour2', 'normal2') and
            refSamples = c(NA, 'normal2', 'drop'),
            tumour1 will be kept as is, since it does not have a matched normal
            tumour2 will be divided by its matched reference: normal2 and
            normal2 will be dropped from further analysis."

            stop(errorMsg)
        }else{
            refSampleIndex = vector("numeric", length = length(bamnames))
            for(i in 1:length(refSamples)){
                ref <- refSamples[i]
                if(ref %in% c(NA, 'NA')){
                    refSampleIndex[i] <- NA
                }else if(ref %in% c('drop')){
                    refSampleIndex[i] <- FALSE
                }
                else{
                    refSampleIndex[i] <- which(ref == bamnames)
                }
            }
        }
    }

    if(!is.null(userMadeBins)){
        if(!((genome %in% "hg19") &
             (binSize %in% c(1, 5, 10, 15, 30, 50, 100, 500, 1000)))){
            errorMsg <- "There are no pre-created annotations for the specified"
            if(genome != "hg19")
                errorMsg <- paste(errorMsg, "genome build;")
            if(!(binSize %in% c(1, 5, 10, 15, 30, 50, 100, 500, 1000)))
                errorMsg <- paste(errorMsg, "bin size. Available bin sizes are:
                                  hg19: 1, 5, 10, 15, 30, 50, 100, 500, 1000 Kbp.\n")

            stop(paste(errorMsg,
                       "CNAclinic does not generate the required annotation as this is
                       a time consuming step.
                       Please generate bin annotations using the QDNAseq package
                       and pass in the output via the userMadeBins argument."))
        }
        }else{

            if(is.null(binSize)){
                errorMsg <- "Default genome: hg19, please specify binSize argument"
                stop(errorMsg)
            }
        }

    outputType <- match.arg(outputType,
                            choices=c("CNAclinicData", "QDNAseqCopyNumbers"),
                            several.ok=FALSE)

    ############################################################################
    # Rename arguments to fit QDNAseq parameters
    ############################################################################
    binSize <- binSize[1]
    path <- pathToBams
    type <- typeOfPreMadeBins

    residual <- residualFilter
    blacklist <- blacklistFilter

    chromosomes <- chromosomesFilter
    chromosomes <- unique(c(chromosomes, "MT", "M"))

    span <- spanForLoess

    maxIter <- maxIterForCorrection
    cutoff <- cutoffForCorrection
    variables <- variablesForCorrection
    fit <- NULL

    logTransform <- logTransformForSmoothing
    ############################################################################



    if(is.null(userMadeBins)){

        # Read in pre-made QDNAseq bin annotation
        userMadeBins <- QDNAseq::getBinAnnotations(binSize=binSize,
                                                   genome=genome,
                                                   type=type, path = NULL)



    }

    # Bin read counts
    readCounts <- QDNAseq::binReadCounts(bins=userMadeBins,
                                         bamfiles=bamfiles,
                                         bamnames=bamnames,
                                         cache=cache, force=!cache,
                                         path=path, ext=ext,
                                         phenofile=NULL,
                                         isPaired=isPaired,
                                         isProperPair=isProperPair,
                                         isUnmappedQuery=isUnmappedQuery,
                                         hasUnmappedMate=hasUnmappedMate,
                                         isMinusStrand=isMinusStrand, isMateMinusStrand=isMateMinusStrand,
                                         isFirstMateRead=isFirstMateRead, isSecondMateRead=isSecondMateRead,
                                         isSecondaryAlignment=isSecondaryAlignment,
                                         isNotPassingQualityControls=FALSE,
                                         isDuplicate=isDuplicate,
                                         minMapq=minMapq,
                                         pairedEnds=pairedEnds)

    # Apply Blacklist filters
    readCountsFiltered <- QDNAseq::applyFilters(readCounts,
                                                residual=residual,
                                                blacklist=blacklist,
                                                mappability=mappabilityFilter,
                                                chromosomes=unique(c(chromosomes, "X", "Y")))

    # Estimate the correction for GC content and mappability
    readCountsFiltered <- QDNAseq::estimateCorrection(readCountsFiltered,
                                                      span=span,
                                                      family=familyForLoess,
                                                      adjustIncompletes=TRUE,
                                                      maxIter=maxIter,
                                                      cutoff=cutoff,
                                                      variables=variables)


    rm(list=c('readCounts'))
    gc(verbose=FALSE)

    if(!all(c("X", "Y") %in% chromosomes)){
        # Reverse the filtering for sex-chromosomes
        # to avoid confounding due to gender when estimating loess correction
        readCountsFiltered <- QDNAseq::applyFilters(readCountsFiltered,
                                                    residual=residual,
                                                    blacklist=blacklist,
                                                    mappability=mappabilityFilter,
                                                    chromosomes=NA)
    }



    # Apply the correction for GC content and mappability.
    # This returns a QDNAseqCopyNumbers object which is used to
    # normalise and smooth outliers
    copyNumbers <- QDNAseq::correctBins(readCountsFiltered,
                                        fit=fit,
                                        method=methodOfCorrection,
                                        adjustIncompletes=TRUE)



    if(saveCountData){

        saveRDS(copyNumbers, file=filename)

    }

    rm(list=c('readCountsFiltered'))
    gc(verbose=FALSE)

    #Apply median normalisation (Optional)

    if(!skipMedianNormalization)
        copyNumbers <- QDNAseq::normalizeBins(copyNumbers,
                                              method=methodOfNormalization)

    gc(verbose=FALSE)

    #Smooth outliers (Optional)
    if(!skipOutlierSmoothing)
        copyNumbers <- QDNAseq::smoothOutlierBins(copyNumbers,
                                                  logTransform=logTransform)

    gc(verbose=FALSE)

    # Normalize by reference samples
    if(!is.null(refSamples)){
        copyNumbers <- QDNAseq::compareToReference(
            copyNumbers,
            references=refSampleIndex)
    }

    message("Pre-processing steps are completed.")

    if(outputType == "CNAclinicData"){

        bins <- data.frame(
            chromosome=QDNAseq::chromosomes(copyNumbers),
            start=QDNAseq::bpstart(copyNumbers),
            end=QDNAseq::bpend(copyNumbers),
            usebin=QDNAseq:::binsToUse(copyNumbers))

        sampleNames <- Biobase::sampleNames(copyNumbers)

        copyNumber<-as.data.frame(
            QDNAseq:::log2adhoc(
                Biobase::assayDataElement(copyNumbers, "copynumber")))

        names(copyNumber) <- sampleNames


        return(CNAclinicData(bins=bins, copyNumber=copyNumber))


    }else{
        return(copyNumbers)
    }

    }
