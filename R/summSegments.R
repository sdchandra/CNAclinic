
#' Summarise segments from multiple algorithms at every genomic bin
#'
#' @param object CNAclinicData.
#'
#' @return An updated CNAclinicData object.
#' @export summSegments
#'
#' @examples
#'      \dontrun{
#'       vignette("CNAclinic")
#'      }
#'
#'
setMethod("summSegments", signature(object="CNAclinicData"),
    function(object, segmentsToSummarise=c("CBS", "LACBS", "HMM", "PLS"),
    summaryMethod="mean"){


    segmentsToSummarise <- match.arg(segmentsToSummarise,
        choices = c("CBS", "LACBS", "HMM", "PLS"), several.ok = TRUE)
    summaryMethod <- match.arg(summaryMethod,
        choices=c("mean", "median", "min", "max", "Q1", "Q3"),
        several.ok=FALSE)

    emptySlots <- emptySlots(object)
    if(any(segmentsToSummarise %in% emptySlots)){
        stop(paste("Some values provided for segmentsToSummarise are not available in your data.",
        "Re-do runSegmentation() or provide different segmentsToSummarise.", sep="\n"))
    }

    sampleNames <- CNAclinic::sampleNames(object)

    listSegDF = list()
    for(i in 1:length(segmentsToSummarise)){
       if(("CBS" %in% segmentsToSummarise) & (nrow(segCBS) != 0))
           listSegDF[["CBS"]] <- segCBS(object)
       if(("LACBS" %in% segmentsToSummarise) & (nrow(segLACBS) != 0))
           listSegDF[["LACBS"]] <- segLACBS(object)
       if(("HMM" %in% segmentsToSummarise) & (nrow(segHMM) != 0))
           listSegDF[["HMM"]] <- segHMM(object)
       if(("PLS" %in% segmentsToSummarise) & (nrow(segPLS) != 0))
           listSegDF[["PLS"]] <- segPLS(object)
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

    message("Updating segSummary in the CNAclinicData object output ...")

    segSummary <- matrix(segSummary, ncol = length(sampleNames))
    segSummary <- as.data.frame(segSummary, col.names = sampleNames)
    segSummary[!CNAclinic::usebin(object), ] <- NA
    names(segSummary) <- sampleNames

    CNAclinic::segSummary(object) <- segSummary

    message("Done.")

    return(object)

    })
