
#' Run up to 4 different copy number segmentation algorithms.
#'
#' \code{runSegmentation} Segments binned and normalized
#' copy number values using multiple algorithms and allows the summarization
#' of the segmentation results. Automatically calculates the "gain", "loss" calls
#' for each segment.
#'
#' @param x Either a \code{\link{CNAclinicData}} object from the output of
#' \code{\link{processForSegmentation}} function or a \code{data.frame} with
#' specific columns (see details) or a \code{\link{QDNAseqCopyNumbers}} object.
#' @param isLogTransformed if TRUE (default) data is assumed to be log2
#' transformed. If FALSE, transformation is carried out prior to segmentation.
#' @param genome Genome build used to align sequencing reads.
#' @param segmentType One or more algorithms for segmentation.
#' i.e. \code{c("CBS", "HMM", "PLS", "LACBS")}
#' @param summaryMethod Summarization method for ensemble segmentation.
#' i.e. One of c("mean", "median", "min", "max", "Q1", "Q3")
#' @param segmentsToSummarise X
#' @param callTypeLog2R Segment type used to call CNAs (default "summary")
#' @param callThreshLog2R Thresholds used to call segments as a "loss" or "gain".
#' Defaults to c(-0.15, 0.15)
#' @param minimumBinsPerSegment Minimum number of bins in each segment,
#' default is 3. Argument is specific to PLS algorithm.
#' @param gapLength Minimum length in base-pairs between the two closest loci
#' to consider a region to be a "gap". (default 1000000)
#' Argument is specific to LACBS algorithm.
#' @param normalizeSegmentedBins Normalizes the segmented bins if TRUE (default).
#' @param inter The interval in which the function should search for the
#' normal level. Utilised in \code{normalizeSegmentedBins}.
#' @param alpha Significance levels for the test to accept change-points.
#' Default is 1e-10
#' @param undo.splits A character string specifying how change-points are to be
#' undone, if at all. Default is "sdundo", which undoes splits that are not at
#' least this many SDs apart. Other choices are "prune", which uses a sum of
#' squares criterion, and "none".
#' @param undo.SD The number of SDs between means to keep a split if
#' undo.splits="sdundo". Default is 1.0.
#' @param segmentStatistic Default is "seg.mean".
#'
#' @details If \code{x} is a data.frame it should contain the columns:
#' \code{chromosome}, \code{start}, \code{end} and the optional column
#' \code{usebin} followed by separate columns named after each different sample
#' contaning normalized and log-transformed copy number measurements.
#' \code{chromosome} should be a character vector containing only 1-22, X, Y
#' while all other columns should contain numerical values.
#' The genomic coordinates should belong to fixed-width bins.
#'
#' @return Returns an object of class \code{\link{CNAclinicData}} with
#' segmentation results from the chosen algorithms as well as the
#' summarised result.
#' @importFrom QDNAseq normalizeSegmentedBins
#' @importClassesFrom QDNAseq QDNAseqReadCounts QDNAseqCopyNumbers
#' @export runSegmentation
#'
#' @author Dineika Chandrananda
#' @examples
#'      \dontrun{
#'       vignette("CNAclinic")
#'      }
#'
runSegmentation=function(x,
    isLogTransformed=TRUE, genome=NULL,
    segmentType=c("CBS", "LACBS", "HMM", "PLS"),
    summaryMethod="mean",
    segmentsToSummarise=segmentType,
    callTypeLog2R="summary",
    callThreshLog2R=c(-0.15, 0.15),
    minimumBinsPerSegment=3,
    gapLength=1000000,
    normalizeSegmentedBins=TRUE, inter=c(-0.1, 0.1),
    alpha=1e-10, undo.splits="sdundo", undo.SD=1.0,
    segmentStatistic="seg.mean"){


#############################################################################

    # Check arguments

#############################################################################

    binsToUse <- NULL
    sampleNames <- NULL
    totalBins <- NULL

    if(missing(x))
        stop("x is missing")

    segmentType <- match.arg(segmentType,
        choices=c("CBS", "LACBS", "HMM", "PLS"),
        several.ok=TRUE)
    segmentsToSummarise <- match.arg(segmentsToSummarise,
        choices=c("CBS", "LACBS", "HMM", "PLS"),
        several.ok=TRUE)
    summaryMethod <- match.arg(summaryMethod,
        choices=c("mean", "median", "min", "max", "Q1", "Q3"),
        several.ok=FALSE)

    if(is.null(genome))
        stop("Please provide genome build used to map sequencing data.")
    if(!(genome %in% c("hg19", "hg18", "hg17", "hg16")) &
       ("PLS" %in% segmentType)){
        stop(paste("The PLS algorithm only runs on hg19, hg18, hg17 & hg16.",
        "Remove PLS from segmentType and re-run.", sep="\n"))

    }


    callTypeLog2R <- match.arg(callTypeLog2R,
        choices = c("CBS", "LACBS", "HMM", "PLS", "summary", "none"),
        several.ok = FALSE)

    if(!isLogTransformed %in% c(TRUE, FALSE))
        stop("isLogTransformed should be TRUE or FALSE")

#############################################################################

    # Check input data format

#############################################################################

    if(class(x) == "data.frame"){

        errorMsg <- ""

        # Check if the "chromosome", "start", "end" present in data.frame
        if(!all(c("chromosome", "start", "end") %in% names(x)))
            errorMsg <- paste(errorMsg,
                "\nx should contain the following columns: chromosome, start, end",
                sep="\n")

        if(!("usebin" %in% names(x))){
            warning(paste("\nx does not contain a column named usebin.",
            "Setting usebin=TRUE for all genomic bins", sep="\n"))
            usebin <- rep(TRUE, nrow(x))
            x <- data.frame(usebin, x, stringsAsFactors=FALSE)
        }

        sampleNames <- names(x)[!(names(x) %in%
            c("chromosome", "start", "end", "usebin"))]

        # Check if sample names are unique
        if(length(unique(sampleNames)) != length(sampleNames)){
            errorMsg <- paste(errorMsg,
                "\nThe following sample names in x are duplicated:",
                duplicated(sampleNames), sep="\n")
        }

        # Check if the start, end and sample columns have numeric values
        isNumeric <- apply(x[c("start", "end", sampleNames)], 2,
                           function(x) is.numeric(x))
        if(!all(isNumeric)){
            errorMsg <- paste(errorMsg,
                "\nThe following variables in x should be numeric:",
                paste(c("start", "end", sampleNames)[!isNumeric], collapse=" "),
                sep="\n")
        }

        if(!is.logical(x$usebin) | any(is.na(x$usebin))){
           errorMsg <- paste(errorMsg,
                "\nThe usebin variable in x should be logical with no NA values.",
                sep="\n")
        }

        # Are the genomic bins fixed width?
        # Check if (end - start) are fixed width for all bins

       # if(length(unique(as.numeric(x$start) - as.numeric(x$end))) != 1){
       #     errorMsg <- paste(errorMsg,
       #         "\nAll genomic bins should be of identical width.", sep="\n")
       # }

        if(paste0(errorMsg) != "")
            stop(paste0(errorMsg))

    }else if(class(x) == "QDNAseqCopyNumbers"){

        sampleNames <- Biobase::sampleNames(x)

        copyNumber<-as.data.frame(
            QDNAseq:::log2adhoc(Biobase::assayDataElement(x, "copynumber")))

        names(copyNumber) <- sampleNames

        tempDF <- data.frame(
            chromosome=QDNAseq::chromosomes(x),
            start=QDNAseq::bpstart(x),
            end=QDNAseq::bpend(x),
            usebin=QDNAseq:::binsToUse(x),
            copyNumber,
            stringsAsFactors=FALSE)

        x <- tempDF

        rm(list=c("tempDF", "copyNumber")); gc(FALSE)

    }else if(class(x) == "CNAclinicData"){

        if(any(c("bins", "copyNumber") %in% emptySlots(x)))
            stop("x should contain genomic bin and copy number data")

        tempDF <- data.frame(
            CNAclinic::bins(x),
            CNAclinic::copyNumber(x),
            stringsAsFactors=FALSE)

        x <- tempDF

        rm(list=c("tempDF")); gc(FALSE)

    }else{
            stop(paste("\nx should either be a CNAclinicData object,",
            "a data.frame with specific columns or a QDNAseqCopyNumbers object",
            sep="\n"))
    }

#############################################################################

    # Transform input data to correct format

#############################################################################

    # At this point x is a data.frame with the genomic regions & log2R values
    # with chromosome, start, end, usebin and sample columns

    # Check chromosome names, remove 'chr' suffix
    # convert to character 1:22, X, Y, MT and
    # set usebin=FALSE for other types of chromosomes

    binsToUse <- x$usebin
    tempDF <- .convertChrom(x$chromosome)
    x$chromosome <- tempDF$chromosome
    x$usebin <- binsToUse & tempDF$valid

    rm(list=c("binsToUse"))
    rm(list=c("tempDF")); gc(FALSE)

    # Remember to backtransform any logged counts as this is expected downstream
    # Transform  any unlogged counts

    logTransform = FALSE
    copyNumber <- NULL

    sampleNames <- names(x)[!(names(x) %in%
        c("chromosome", "start", "end", "usebin"))]


    x[!x$usebin, sampleNames] <- NA_real_

    if(!is.null(isLogTransformed)){
        if(!isLogTransformed){
            message(paste("\nisLogTransformed set to FALSE",
                "Converting ratios in x to log2 values", sep="\n"))

            x[sampleNames] <- .toLog2(x[sampleNames])
        }
    }

    # Get  totalBins
    totalBins <- nrow(x)

    x <- .reorderByChrom(x)
    binsToUse <- x$usebin

    #Set outlier counts to NA
    #tempDF <- x[sampleNames]
    #tempDF[abs(tempDF) >= outlierLog2R] <- NA
    #x[sampleNames] <- tempDF


#############################################################################

    # Initialize objects

#############################################################################

    # as.data.frame(setNames(replicate(length(sampleNames),numeric(0), simplify = F), sampleNames))
    segCBS <- segLACBS <- segHMM <- segPLS <- segSummary <- calls <- data.frame()

#############################################################################

    # Run the segmentation algorithms

#############################################################################

    if("HMM" %in% segmentType){
        message("*** ***")
        message("Running HMMcopy segmentation")
        flag <- TRUE
        tryCatch({
            segHMM <- .getHMMcopySegments(x)
        }, error=function(e) {
            message(paste(
            "\nThe following error occured in HMMcopy::HMMsegment()",
            e, "The HMM algorithm will be skipped", sep="\n"))
            flag <-FALSE
        })

        if(flag){
            if(!all(dim(segHMM) == c(totalBins, length(sampleNames)))){
                stop(paste("\nHMMcopy segmentation did not return the expected",
                "    number of bins or samples", sep="\n"))
            }
        }
    }


    if("CBS" %in% segmentType){
        message("*** ***")
        message("Running DNAcopy CBS segmentation")
        flag <- TRUE
        tryCatch({

            segDF <- .getCBSSegments(x,
                alpha=alpha,
                undo.splits=undo.splits,
                undo.SD=undo.SD,
                segmentStatistic=segmentStatistic,
                normalizeSegmentedBins=normalizeSegmentedBins,
                inter=inter)

        }, error=function(e) {
            message(paste(
            "\nThe following error occured in DNAcopy::segment()",
            e, "The CBS algorithm will be skipped", sep="\n"))

            flag <-FALSE
        })

        if(flag){
            segCBS <- as.data.frame(segDF, stringsAsFactors=FALSE)

            if(!all(dim(segCBS) == c(totalBins, length(sampleNames)))){
                stop(paste("\nCBS segmentation did not return the expected",
                "    number of bins or samples", sep="\n"))
            }
        }
    }

    if("LACBS" %in% segmentType){
        message("*** ***")
        message("Running Locus-Aware CBS segmentation")
        flag <- TRUE
        tryCatch({

            segLACBS <- x[sampleNames]
            segLACBS[sampleNames] <- NA

            for(i in seq_along(sampleNames)){
                message(paste0(sampleNames[i]))

                pscbsInput <- data.frame(
                    chromosome=as.character(x$chromosome),
                    x=as.numeric(x$start),
                    y=x[sampleNames[i]],
                    index=seq_len(nrow(x)),
                    stringsAsFactors=FALSE)

                names(pscbsInput) <- c("chromosome", "x", "y")

                # Convert chromosomes to integer
                pscbsInput$chromosome[pscbsInput$chromosome == "X"] = "23"
                pscbsInput$chromosome[pscbsInput$chromosome == "Y"] = "24"
                pscbsInput$chromosome[pscbsInput$chromosome == "MT"] = "25"

                pscbsInput$chromosome <- as.integer(pscbsInput$chromosome)

                pscbsInput <- PSCBS::dropSegmentationOutliers(pscbsInput)
                gaps <- PSCBS::findLargeGaps(pscbsInput, minLength=gapLength)

                # Segmenting multiple chromosomes at once
                knownSegments <- PSCBS::gapsToSegments(gaps)

                lacbsFit <- PSCBS::segmentByCBS(pscbsInput, undo=1.0,
                    knownSegments=knownSegments, verbose=0, seed=48879)

                rm(list=c("pscbsInput")); gc(FALSE)

                segResults <- na.omit(data.frame(
                    startRow=lacbsFit$segRows$startRow,
                    endRow=lacbsFit$segRows$endRow,
                    segValue=lacbsFit$output$mean,
                    stringsAsFactors=FALSE))

                rm(list=c("lacbsFit")); gc(FALSE)

                for(j in 1:nrow(segResults)){
                    startRow <- segResults$startRow[j]
                    endRow <- segResults$endRow[j]
                    segLACBS[startRow:endRow, sampleNames[i]] <- segResults$segValue[j]
                }
            }
        }, error=function(e) {
            message(paste(
            "\nThe following error occured in the PSCBS package",
            e, "The LACBS algorithm will be skipped", sep="\n"))

            flag <-FALSE
        })

        if(flag){
            if(!all(dim(segLACBS) == c(totalBins, length(sampleNames)))){
                stop(paste("\nLACBS segmentation did not return the expected",
                "number of bins or samples", sep="\n"))
            }
        }
    }

    if("PLS" %in% segmentType){
        message("*** ***")
        message("Running copynumber PLS segmentation")
        flag <- TRUE
        copynumberInput <- data.frame(
            chrs=x$chromosome,
            pos=x$start,
            x[sampleNames])
        names(copynumberInput) <- c("chrs", "pos", sampleNames)

        tryCatch({
            PLSvalue <- copynumber::pcf(copynumberInput, pos.unit = "bp",
                arms=NULL, Y=NULL, kmin = minimumBinsPerSegment,
                gamma=40,
                normalize=TRUE,
                fast=TRUE,
                assembly=genome,
                digits=4,
                return.est=TRUE,
                save.res=FALSE,
                file.names=NULL,
                verbose=FALSE)

        }, error=function(e) {
            message(paste(
            "\nThe following error occured in copynumber::pcf()",
            e, "The PLS algorithm will be skipped", sep="\n"))
            flag <-FALSE
            rm(list=c("copynumberInput")); gc(FALSE)
        })

        if(flag){
            segPLS <- PLSvalue$estimates[ ,
                        !(names(PLSvalue$estimates) %in% c("chrom", "pos")), drop=FALSE]
            rm(list=c('copynumberInput', 'PLSvalue'))
            gc(FALSE)

            if(!all(dim(segPLS) == c(totalBins, length(sampleNames)))){
                stop(paste("\nPLS segmentation did not return the expected",
                "number of bins or samples", sep="\n"))
            }
        }
    }

    message("Segmentation complete.")

    listSegDF = list()
    for(i in 1:length(segmentsToSummarise)){
        if(("CBS" %in% segmentsToSummarise) & (nrow(segCBS) != 0))
            listSegDF[["CBS"]] <- segCBS
        if(("LACBS" %in% segmentsToSummarise) & (nrow(segLACBS) != 0))
            listSegDF[["LACBS"]] <- segLACBS
        if(("HMM" %in% segmentsToSummarise) & (nrow(segHMM) != 0))
            listSegDF[["HMM"]] <- segHMM
        if(("PLS" %in% segmentsToSummarise) & (nrow(segPLS) != 0))
            listSegDF[["PLS"]] <- segPLS
    }

    listVec <- lapply(listSegDF, c, recursive=TRUE)
    listMatrix <- do.call(cbind, listVec); rm(list=c('listVec')); gc(FALSE)


    message(paste("Calculating", summaryMethod, "of",
        paste(names(listSegDF), collapse = ","), "segments per genomic bin"))

    segSummary <- switch(summaryMethod,
           mean = apply(listMatrix, 1, mean, na.rm=TRUE),
           median = apply(listMatrix, 1, median, na.rm=TRUE),
           min = apply(listMatrix, 1, function(x) {
               x[which.min(abs(x))]
               }),
           max = apply(listMatrix, 1, function(x) {
               x[which.max(abs(x))]
           }),
           Q1 = apply(listMatrix, 1, function(x) quantile(x, na.rm = TRUE)[2]),
           Q2 = apply(listMatrix, 1, function(x) quantile(x, na.rm = TRUE)[4]))

    segSummary <- matrix(segSummary, ncol=length(sampleNames))

    message("Creating CNAclinicData object as output...")

    bins <- data.frame(chromosome = x$chromosome,
                      start = x$start,
                      end = x$end,
                      usebin = binsToUse,
                      stringsAsFactors = FALSE)

    copyNumber <- x[sampleNames]

    segSummary <- as.data.frame(segSummary, col.names=sampleNames)
    segSummary[!binsToUse, ] <- NA
    names(segSummary) <- sampleNames


    if(("CBS" %in% segmentType) & (nrow(segCBS) != 0)){
        segCBS <- as.data.frame(segCBS, col.names=sampleNames)
        segCBS[!binsToUse, ] <- NA
    }
    if(("LACBS" %in% segmentType) & (nrow(segLACBS) != 0)){
        segLACBS <- as.data.frame(segLACBS, col.names = sampleNames)
        segLACBS[!binsToUse, ] <- NA
    }
    if(("HMM" %in% segmentType) & (nrow(segHMM) != 0)){
        segHMM <- as.data.frame(segHMM, col.names = sampleNames)
        segHMM[!binsToUse, ] <- NA
    }
    if(("PLS" %in% segmentType) & (nrow(segPLS) != 0)){
        segPLS <- as.data.frame(segPLS, col.names = sampleNames)
        segPLS[!binsToUse, ] <- NA
    }

    message("Done.")

    CNAdata <- CNAclinicData(
        bins=bins,
        copyNumber=copyNumber,
        segCBS=segCBS,
        segLACBS=segLACBS,
        segHMM=segHMM,
        segPLS=segPLS,
        segSummary=segSummary)

    if(callTypeLog2R != "none"){

        CNAdata <- callData(CNAdata,
            callTypeLog2R=callTypeLog2R,
            callThreshLog2R=callThreshLog2R)

    }

    return(CNAdata)

}

#############################################################################

# Helper function to smooth outlier bins post-segmentation

#############################################################################

.smoothOutliers <- function(x, bins){

    segValues <- as.matrix(x)

    # Extract annotation data
    #fData <- fData(object)

    # Sanity check
    #stopifnot(is.matrix(copynumber))

    # Log transform?
    #if (logTransform)
    #    copynumber <- log2adhoc(copynumber)

    # Filter
    condition <- bins$usebin

    #vmsg("Smoothing outliers ...", appendLF=FALSE)

    CNA.object <- DNAcopy::CNA(segValues[condition, , drop=FALSE],
                      chrom=bins$chromosome[condition],
                      maploc=bins$start[condition],
                      data.type="logratio", presorted=TRUE)
    CNA.object <- DNAcopy::smooth.CNA(CNA.object)
    CNA.object <- CNA.object[, -(1:2), drop=FALSE]
    segValues <- as.matrix(CNA.object)
    CNA.object <- NULL


    # Expand to full set of bins
    segValuesTemp <- matrix(NA_real_, nrow=nrow(x), ncol=ncol(x),
                          dimnames=list(row.names(x), names(x)))
    segValuesTemp[condition, ] <- segValues
    segValues <- NULL

   return(segValuesTemp)


}

#############################################################################

# Helper functions for HMM & PLS segmentation

#############################################################################


.splitCollatedSegs <- function(x, binSize=NULL, segmentType=NULL){


    if((segmentType == "HMM") & is.null(binSize))
        stop("Provide binSize argument")
    if(is.null(segmentType) || length(segmentType) > 1)
        stop("Provide 1 segmentType argument")

    if(segmentType == "PLS"){

        x <- data.frame(
            chr = x$chrom,
            start = x$start.pos,
            end = x$end.pos,
            n.probes = x$n.probes,
            median = x$mean)

    }

    chromosome <- as.character(x$chr)
    chromosome[which(chromosome == "X")] <- "23"
    chromosome[which(chromosome == "Y")] <- "24"
    chromosome[which(chromosome == "MT")] <- "25"
    x$chr <- as.numeric(chromosome)


    if(!all(is.numeric(unique(x$chr))))
       stop("Chromosome names are not numeric. e.g. 1,2,..,23,24,25")

    # Not assuming chromosomes have already been ordered

    x <- x[order(x[,"chr"], x[,"start"]), ]

    binState <- c()
    binSegLog2R <- c()

    for(i in 1:nrow(x)){
        if(segmentType == "HMM"){
            nBinsInSegment <- ceiling(((x$end[i] - x$start[i]) + 1) / binSize)
        }else if(segmentType == "PLS"){
            nBinsInSegment <- x$n.probes[i]
        }
        binSegLog2R <- c(binSegLog2R, rep(x$median[i], nBinsInSegment))

        if(segmentType == "HMM"){
            binState <- c(binState, rep(x$state[i], nBinsInSegment))
        }else{
            binState <- c(binState, rep(NA, nBinsInSegment))

        }
    }
    data.frame(binSegLog2R, binState, stringsAsFactors=FALSE)
}



.getHMMcopySegments = function(x){

    # x should be a data.frame()

    sampleNames <- names(x)[!(names(x) %in%
        c("chromosome", "start", "end", "usebin"))]

    copyNumber <- x[sampleNames]

    # Convert all log2R values greater than +/- 10 to NA

    copyNumber[abs(copyNumber) >= 10] <- NA_real_
    copyNumber[!x$usebin, ] <- NA_real_

    binSize <- unique(diff(x$end[1:10]))

    if(length(binSize) != 1)
        stop("bin size in x is not fixed.")

    hmmSegsCalls <- as.data.frame(
        matrix(NA, nrow=nrow(x), ncol = length(sampleNames) * 2))

    for(i in 1:length(sampleNames)){
        message(sampleNames[i])
        sampleRangedData <- IRanges::RangedData(
            ranges=IRanges::IRanges(start=x$start, width=(x$end - x$start) + 1),
            space=x$chromosome,
            copy=copyNumber[ , sampleNames[i]])

        sampleHMMData = suppressMessages(HMMcopy::HMMsegment(sampleRangedData, verbose = TRUE))

        hmmSegsCalls[ ,((i*2)-1):(i*2)] <- .splitCollatedSegs(
            x=sampleHMMData$segs, binSize=binSize, segmentType="HMM")
    }

    hmmSegsCalls <- data.frame(hmmSegsCalls, stringsAsFactors=FALSE)
    names(hmmSegsCalls) <- paste(rep(sampleNames, each = 2),
        c("", "_hmmCall"), sep = "")

    hmmSegsCalls[ , !grepl("_hmmCall", names(hmmSegsCalls)), drop=FALSE]
}



#############################################################################

# Wrapper function for CBS segmentation

#############################################################################

.getCBSSegments <- function(x, alpha=1e-10, undo.splits="sdundo", undo.SD=1.0,
    segmentStatistic="seg.mean", normalizeSegmentedBins=FALSE,
    inter=c(-0.1, 0.1)){

    condition <- x$usebin

    sampleNames <- names(x)[!(names(x) %in%
        c("chromosome", "start", "end", "usebin"))]

    copyNumber <- x[sampleNames]
    copyNumber[!condition, ] <- NA_real_

    # To Stabilize variance for CBS
    copyNumber <- .toLog2(copyNumber, inv=TRUE)
    copyNumber <- .stabilisedTrans(copyNumber, inv=FALSE)
    copyNumber <- as.matrix(copyNumber)

    chromosome <- x$chromosome
    start <- x$start

    # Create a list of CNA objects that can be analyzed with *lapply()
    cna <- lapply(sampleNames, function(x)
        DNAcopy::CNA(
            genomdat=copyNumber[condition, x, drop=FALSE],
            chrom=factor(chromosome[condition], levels=unique(chromosome),
                         ordered=TRUE),
            maploc=start[condition],
            data.type="logratio",
            sampleid=x,
            presorted=TRUE))

    rm(list=c("chromosome", "start")); gc(FALSE)

    # Create a vector of messages to be printed
    msgs <- paste0(sampleNames)

    # Use sample names for indexing, they are available in the CNA objects
    # as the name of the third column
    names(msgs) <- sampleNames
    segments <- .flapply(cna, FUN=function(x, ...) {
        message(msgs[colnames(x)[3]])
        DNAcopy::segment(x, alpha=alpha, undo.splits=undo.splits,
                         undo.SD=undo.SD, verbose=0, ...)
    })

    segmentStatisticCol <- grep(segmentStatistic,
        colnames(DNAcopy::segments.summary(segments[[1]])))

    segDF <- matrix(NA_real_, nrow=nrow(copyNumber), ncol=ncol(copyNumber))

    segDF[condition, ] <- do.call(cbind, lapply(segments, function(x)
        rep(DNAcopy::segments.summary(x)[,segmentStatisticCol],
            x$output$num.mark)))

    colnames(segDF) <- sampleNames
    row.names(segDF) <- NULL
    segDF[is.na(copyNumber)] <- NA_real_
    segDF <- .toLog2(.stabilisedTrans(segDF, inv=TRUE), inv=FALSE)

    rm(list=c("cna")); gc(FALSE)

    if(normalizeSegmentedBins){

        x <- data.frame(
            chromosome=x$chromosome,
            start=x$start,
            end=x$end,
            usebin=x$usebin,
            segDF,
            stringsAsFactors=FALSE)
        names(x) <- c("chromosome", "start", "end", "usebin", sampleNames)


        seg <- .cghSegConversion(x, copyNumber)
        rm(list=c("copyNumber")); gc(FALSE)

        # Normalize
        postseg <- CGHcall::postsegnormalize(seg, inter=inter)

        # Return the normalized segments and copynumber values in a list
        segDF <- matrix(NA_real_, nrow=nrow(segDF), ncol=ncol(segDF))
        segDF[condition, ] <- CGHbase::segmented(postseg)

        #copyNumber[condition, ] <- CGHbase::copynumber(postseg)

        row.names(segDF) <- NULL
        #row.names(copyNumber) <- NULL

        colnames(segDF) <- sampleNames
        #colnames(copyNumber) <- sampleNames

    }

    return(segDF)
    #return(list(segDF=segDF, copyNumber=copyNumber))

}

#############################################################################

# Helper functions needed for chromosome formatting

#############################################################################

.convertChrom <- function(chromosome){

    # If factor or numeric; need to convert to character first
    chromosome <- as.character(chromosome)

    if(all(grepl("chr", chromosome))){
        chromosome <- unlist(lapply(strsplit(chromosome, "chr"), function(x) x[2]))
    }else if(!sum(grepl("chr", chromosome)) %in% c(0, length(chromosome))){
        stop("Chromosome names in dataset should not contain `chr` suffix")
    }

    # Replace 23 by X:
    chromosome[which(chromosome == "23")] <- "X"

    # Replace 24 by Y
    chromosome[which(chromosome == "24")] <- "Y"

    # Replace 25 by MT
    chromosome[which(chromosome == "25")] <- "MT"

    # Replace M by MT
    chromosome[which(chromosome == "M")] <- "MT"


    # TRUE if chromosome is acceptable
    valid <- (chromosome %in% c(1:22, "X", "Y"))

    return(data.frame(chromosome, valid, stringsAsFactors=FALSE))
}

.reorderByChrom <- function(x){

    # x is a data.frame(chromosome, start, end, usebin, samp1, samp2 ...)
    # The chromosomes are assumed to be character and only inluclude X, Y, MT
    # Order chromosomes & the data from 1:22, 23, 24, 25
    # Re-convert the chromosomes back to 1:22, X, Y, MT


    chromosome <- as.character(x$chromosome)
    chromosome[which(chromosome == "X")] <- "23"
    chromosome[which(chromosome == "Y")] <- "24"
    chromosome[which(chromosome == "MT")] <- "25"

    x$chromosome <- as.numeric(chromosome)

    x <- x[order(x["chromosome"], x["start"]), ]

    x$chromosome <- as.character(x$chromosome)

    # Replace 23 by X:
    x$chromosome[which(x$chromosome == "23")] <- "X"

    # Replace 24 by Y
    x$chromosome[which(x$chromosome == "24")] <- "Y"

    # Replace 25 by MT
    x$chromosome[which(x$chromosome == "25")] <- "MT"

    return(x)

}


#############################################################################

# Helper functions needed for CBS segmentation adapted from QDNAseq

#############################################################################


# Create a cghSeg object for post CBS normalization

.cghSegConversion = function(x, copyNumber){

    sampleNames <- names(x)[!(names(x) %in%
                                  c("chromosome", "start", "end", "usebin"))]


    x_original <- x

    className <- 'cghSeg'

    condition <- x$usebin
    x <- x[condition, ]
    copynumber <- copyNumber[condition, , drop=FALSE]
    condition <- NULL; gc(FALSE)

    chromosome <- as.character(x$chromosome)
    chromosome[chromosome == "X"] <- "23"
    chromosome[chromosome == "Y"] <- "24"
    chromosome[chromosome == "MT"] <- "25"
    x$chromosome <- as.integer(chromosome)

    if (any(is.na(x$chromosome)))
        stop(paste0("Unknown chromosome names:\n",
                    paste(unique(x$chromosome[is.na(x$chromosome)]), collapse=", ")))

    # Update column names
    names <- colnames(x)
    names[names == 'chromosome'] <- 'Chromosome'
    names[names == 'start'] <- 'Start'
    names[names == 'end'] <- 'End'
    colnames(x) <- names

    phenodata <- data.frame(name=sampleNames,
                            row.names=sampleNames,
                            stringsAsFactors=FALSE)

    featuredata <- data.frame(Chromosome=x$Chromosome,
                              Start=x$Start,
                              End=x$End,
                              stringsAsFactors=FALSE)

    segmented=x[sampleNames]

    scipen <- options("scipen")
    options(scipen=10)
    featureNames <- paste(x$Chromosome, ":", x$Start, "-", x$End, sep="")

    row.names(copynumber) <- featureNames
    row.names(featuredata) <- featureNames
    row.names(segmented) <- featureNames
    row.names(x) <- featureNames


    # Instantiate a cghSeg object
    seg <- new('cghSeg',
               copynumber=copynumber,
               segmented=segmented,
               featureData=Biobase::AnnotatedDataFrame(featuredata),
               phenoData=Biobase::AnnotatedDataFrame(phenodata))

    # Remove unwanted objects
    rm(list=c("copynumber", "featuredata", "segmented", "phenodata"))
    #gc(FALSE)
    options(scipen=scipen)

    return(seg)

}

.flapply <- function(X, FUN, ..., seeds=NULL) {

    if (!is.null(seeds)) {
        if (length(seeds) < length(X))
            seeds <- rep_len(seeds, length.out=length(X))
    } else {
        seeds <- rep.int(NA_integer_, times=length(X))
    }

   # if(length(X) == 1 || !"future" %in% loadedNamespaces()) {
        res <- list()
        for (ii in seq_along(X)) {
            if (!is.na(seeds[ii]))
                set.seed(seeds[ii])
            res[[ii]] <- FUN(X[[ii]], ...)
        }
        names(res) <- names(X)
        return(res)
   # }
}

.fapply <- function(X, MARGIN, FUN, ...) {
    return(apply(X, MARGIN, FUN, ...))
}

.toLog2 = function(x, inv=FALSE) {

    offset=.Machine$double.xmin
    offset <- as.double(offset)
    stopifnot(is.finite(offset))

    if (!inv) {
        x[x < 0] <- 0
        x <- x + offset
        log2(x)
    } else {
        x <- 2^x
        x - offset
    }

}

.stabilisedTrans = function(x, inv=FALSE){

    factor=as.double(3/8)
    offset=as.double(0)

    if (!inv) {
        x <- x + offset
        sqrt(x * factor)
    } else {
        x <- x^2 * (factor^-1)
        x - offset
    }

}

# EOF

