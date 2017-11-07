#' Combine the slots of two CNAclinicData objects
#'
#' This method uses cbind to create a new CNAcliniData object that contains the
#' sample data of both arguments x and y. Both CNAcliniData objects to be
#' combined (x, y)  must have the same collection of slots filled and have
#' identical bin annotations.
#'
#' @param x a \code{\link{CNAclinicData}} object.
#' @param y a \code{\link{CNAclinicData}} object.
#'
#' @return Returns an object of class \code{\link{CNAclinicData}} with
#' sample data of both \code{x} and \code{y}.
#'
#' @export combineData
#'
#' @examples
#'      \dontrun{
#'       vignette("CNAclinic")
#'      }
#'

setMethod("combineData", signature(x="CNAclinicData", y="CNAclinicData"),
    function(x, y) {

        if(class(x) != class(y)) {
            errorMsg <- "'x & y objects have different classes '%s', '%s'"
            stop(sprintf(errorMsg, class(x), class(y)))
        }

        sampleNames = c(sampleNames(x), sampleNames(y))
        duplicated <- duplicated(sampleNames)

        if(sum(duplicated) != 0){
            stop("The following sample names are duplicated between x & y :",
                sampleNames[duplicated])
        }

        equalBins <- all.equal(bins(x), bins(y))
        if(!isTRUE(equalBins))
            stop("The genomic bins of the CNAclinicObjects x & y are not equal\n
                Check these values using the bins() function.")

        filledSlotsX <- filledSlots(x)
        filledSlotsY <- filledSlots(y)

        if((length(filledSlotsX) != length(filledSlotsY)) ||
            !all(filledSlotsX %in% filledSlotsY)){

            stop("The two CNAclinicObjects should have the same type of
            filled slots to be combined.
            Use the filledSlots() function to investigate.")

        }

        z <- CNAclinicData(bins=bins(x))

        filledSlotsX <- filledSlotsX[!(filledSlotsX %in% "bins")]

        for(s in filledSlotsX){

            df <- cbind(slot(x, s), slot(y, s))
            names(df) <- sampleNames

            slot(z, s) <- df

        }

        return(z)
    })

# EOF
