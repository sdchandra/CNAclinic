
#' Extract data belonging to a subset of samples in a CNAclinic object
#'
#' This function allows splitting up a CNAclinic object into different sample
#' subsets.
#'
#' @param object a \code{\link{CNAclinicData}} object
#' @param sampleNames a subset of sample names in \code{x}
#'
#' @return Returns an object of class \code{\link{CNAclinicData}} with data
#' belonging to the \code{sampleNames}
#'
#' @export subsetData
#'
#' @examples
#'      \dontrun{
#'       vignette("CNAclinic")
#'      }
#'

setMethod("subsetData",
    signature(object="CNAclinicData"),
    function(object, sampleNames=CNAclinic::sampleNames(object)){

    # Check if sampleNames isn't duplicated
        if(any(duplicated(sampleNames)))
            stop("sampleNames must be unique.")

        allsampleNames <- sampleNames(object)
        badNames <- !(sampleNames %in% allsampleNames)
        if(any(badNames))
            stop("CNAclinic object does not contain the sampleNames requested.")

        dropNames<- allsampleNames[!allsampleNames %in% sampleNames]
        if (0 == length(dropNames)){
            message("Specify sampleNames. No subsetting done.")

            return(object)
        }

        x <- CNAclinicData()

        for(s in filledSlots(object)){
            # Use list subsetting in order to keep slots as data.frame even
            # with 1 sample

            if(s %in% "bins")
                slot(x, s) <- slot(object, s)
            else
                slot(x, s) <- (slot(object, s))[as.character(sampleNames)]

        }
        return(x)
    })
