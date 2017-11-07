
#' An S4 class to contain the copy number, segmentation and calls from
#' multiple samples.
#'
#'  An object of this class contains several data.frame objects.
#'  \code{bins} contains information on binned genomic windows such as the
#'  chromosome, start and end of the bins as well as a condition to use or
#'  filter the bins. All other data.frames have rows corresponding to genomic
#'  bins and the columns to separate samples.
#'  Of these, \code{copyNumber} contains the bias corrected and normalized
#'  count data (log2 ratios), the others named \code{segCBS}, \code{segHMM},
#'  \code{segPLS} and \code{segLACBS} contain the value of the segment for
#'  each genomic bin calculated from multiple segmentation algorithms.
#'  The \code{segSummary} data.frame holds the summarised consensus segment
#'  values for each bin and input sample, while \code{calls} holds the
#'  respective copy number calls.
#'
#'  @name CNAclinicData
#'
#'  @slot bins A \code{\link[base]{data.frame} object containing bin annotations
#'   with columns: chromosome, start, end & usebin.
#'  @slot copyNumber A \code{\link[base]{data.frame} object containing bias
#'  corrected and normalized binned log2 ratios.
#'  @slot segCBS A \code{\link[base]{data.frame} object containing results of
#'  running Circular Binary Segmentation algorithm.
#'  @slot segLACBS A \code{\link[base]{data.frame} object containing results of
#'  running Locus-Aware Circular Binary Segmentation algorithm.
#'  @slot segHMM A \code{\link[base]{data.frame} object containing results of
#'  running a Hidden Markov Model segmentation algorithm.
#'  @slot segPLS A \code{\link[base]{data.frame} object containing results of
#'  running Penalized Least Squares segmentation algorithm.
#'  @slot segSummary A \code{\link[base]{data.frame} object containing
#'  summarised segment values.
#'  @calls A data.frame object containing the copy number calls as -1 = loss,
#'  0 = neutral & 1 = gain.
#'
#'  @author Dineika Chandrananda
#'
#'  @examples
#'      showClass("CNAclinicData")
#'      \dontrun{
#'      data(CNAData)
#'      show(CNAData)
#'      }
#'  }
#'
#' @exportClass CNAclinicData
#' @export CNAclinicData
#' @export chromosomes bpstart bpend usebin sampleNames bins copyNumber
#' @export segCBS segLACBS segHMM segPLS segSummary calls
#' @export emptySlots filledSlots
#' @export usebin<- sampleNames<- copyNumber<-
#' @export segHMM<- segCBS<- segLACBS<- segPLS<- segSummary<- calls<-
#' @export subsetData combineData
#' @export plotSampleData plotMultiSampleData getGeneInfo callData summSegments
#'
#' @rdname CNAclinicData
#'

# Access methods

setMethod("chromosomes", signature=c(object="CNAclinicData"),
    definition=function(object) {
        (object@bins)$chromosome
    })

setMethod("bpstart", signature=c(object="CNAclinicData"),
    definition=function(object) {
        (object@bins)$start
    })

setMethod("bpend", signature=c(object="CNAclinicData"),
    definition=function(object) {
        (object@bins)$end
    })

setMethod("usebin", signature=c(object="CNAclinicData"),
    definition=function(object) {
        (object@bins)$usebin
    })

setMethod("sampleNames", signature=c(object="CNAclinicData"),
    definition=function(object) {
        nonEmpty <- (slotNames(object))[slotNames(object) %in% filledSlots(object)]
        names(slot(object, (nonEmpty[!(nonEmpty %in% "bins")])[1]))

    })

setMethod('bins', signature=c(object="CNAclinicData"),
          definition=function(object) {
              object@bins
          })

setMethod('copyNumber', signature=c(object="CNAclinicData"),
    definition=function(object) {
      object@copyNumber
    })

setMethod('segHMM', signature=c(object="CNAclinicData"),
    definition=function(object) {
        object@segHMM
    })

setMethod('segCBS', signature=c(object="CNAclinicData"),
    definition=function(object) {
      object@segCBS
    })

setMethod('segLACBS', signature=c(object="CNAclinicData"),
          definition=function(object) {
              object@segLACBS
          })

setMethod('segPLS', signature=c(object="CNAclinicData"),
    definition=function(object) {
      object@segPLS
    })

setMethod('segSummary', signature=c(object="CNAclinicData"),
    definition=function(object) {
      object@segSummary
    })

setMethod('calls', signature=c(object="CNAclinicData"),
    definition=function(object) {
      object@calls
    })


setMethod('emptySlots', signature=c(object="CNAclinicData"),
    definition=function(object) {
        slotNames <- slotNames(object)
        rowsOfSlots <- "c("

        for(s in 1:length(slotNames)){
            rowsOfSlots <- paste(rowsOfSlots, "nrow(object@", slotNames[s], ")",
                                 sep = "")
            if(s != length(slotNames)){
                rowsOfSlots <- paste(rowsOfSlots, ",", sep = "")
            }else{
                rowsOfSlots <- paste(rowsOfSlots, ")", sep = "")
            }
        }

        rowN <- eval(parse(text = rowsOfSlots))
        slotNames(object)[rowN == 0]
    })


setMethod('filledSlots', signature=c(object="CNAclinicData"),
    definition=function(object) {
        slotNames <- slotNames(object)
        rowsOfSlots <- "c("

        for(s in 1:length(slotNames)){

            rowsOfSlots <- paste(rowsOfSlots, "nrow(object@", slotNames[s], ")",
                                 sep = "")
            if(s != length(slotNames)){
                rowsOfSlots <- paste(rowsOfSlots, ",", sep = "")
            }else{
                rowsOfSlots <- paste(rowsOfSlots, ")", sep = "")
            }
        }

        rowN <- eval(parse(text = rowsOfSlots))
        slotNames(object)[rowN != 0]
    })


# Replacement methods

setMethod('usebin<-', 'CNAclinicData',
    function(object, value) {

        stopifnot(all(is.logical(value)))

        object@bins <- data.frame(
            bins(object)[ , c("chromosome", "start", "end")],
            usebin=value, stringsAsFactors=FALSE)

        if(validObject(object))
            return(object)

    })



# setReplaceMethod('usebin',
#                  signature=c(object='CNAclinicData', value='logical'),
#                  definition=function(object, value) {
#                      usebin <- usebin(object)
#
#                      bins(object) <-
#                          data.frame(bins(object)[ , c("chromosome", "start", "end")],
#                                     usebin, stringsAsFactors=FALSE)
#
#                      if(validObject(object))
#                          return(object)
#                  })


#setReplaceMethod('bins',
#                 signature=c(object='CNAclinicData', value='data.frame'),
#                 definition=function(object, value){

#                    object@bins <- value

#                    if(validObject(object))
#                         return(object)
#                 })

setReplaceMethod("sampleNames",
    signature(object='CNAclinicData', value="ANY"),
    definition=function(object, value){

        if(!is.null(value) && (length(value) != dim(object@copyNumber)[[2]])){
            stop("number of new names (", length(value), ") ",
                "should equal number of columns in copyNumber (",
                dim(object@copyNumber)[[2]], ")")
        }else{
            names(object@copyNumber) <- value
            names(object@segCBS) <- value
            names(object@segHMM) <- value
            names(object@segPLS) <- value
            names(object@segSummary) <- value
            names(object@calls) <- value

            if(validObject(object))
                return(object)
        }
    })


setReplaceMethod('copyNumber',
    signature=c(object='CNAclinicData', value='data.frame'),
    definition=function(object, value) {
        object@copyNumber <- value
        if(validObject(object))
            return(object)
    })


setReplaceMethod('segHMM',
                 signature=c(object='CNAclinicData', value='data.frame'),
                 definition=function(object, value) {
                     object@segHMM <- value
                     if(validObject(object))
                         return(object)
                 })


setReplaceMethod('segCBS',
                 signature=c(object='CNAclinicData', value='data.frame'),
                 definition=function(object, value) {
                     object@segCBS <- value
                     if(validObject(object))
                         return(object)
                 })

setReplaceMethod('segLACBS',
                 signature=c(object='CNAclinicData', value='data.frame'),
                 definition=function(object, value) {
                     object@segLACBS <- value
                     if(validObject(object))
                         return(object)
                 })


setReplaceMethod('segPLS',
                 signature=c(object='CNAclinicData', value='data.frame'),
                 definition=function(object, value) {
                     object@segPLS <- value
                     if(validObject(object))
                         return(object)
                 })


setReplaceMethod('segSummary',
                 signature=c(object='CNAclinicData', value='data.frame'),
                 definition=function(object, value) {
                     object@segSummary <- value
                     if(validObject(object))
                         return(object)
                 })


setReplaceMethod('calls',
                 signature=c(object='CNAclinicData', value='data.frame'),
                 definition=function(object, value) {
                     object@calls <- value
                     if(validObject(object))
                         return(object)
                 })

