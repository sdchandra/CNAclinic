
#' Call segments as 'gain', 'loss' or 'neutral'
#'
#' @param object CNAclinicData.
#'
#' @return An object of class CNAclinicData
#'
#' @export callData
#'
#' @author Dineika Chandrananda
#' @examples
#'      \dontrun{
#'       vignette("CNAclinic")
#'      }
#'
setMethod('callData', signature(object="CNAclinicData"),
    definition=function(object, callTypeLog2R="summary",
    callThreshLog2R=c(-0.15, 0.15)) {

    callTypeLog2R <- match.arg(callTypeLog2R,
        choices = c("CBS", "LACBS", "HMM", "PLS", "summary"),
        several.ok = FALSE)

    slots <- c(CBS="segCBS", LACBS="segLACBS", HMM="segHMM", PLS="segPLS",
        summary="segSummary")

    for(s in names(slots)){
        if((s == callTypeLog2R) & (slots[s] %in% emptySlots(object))){
            errorMsg<- paste(
                paste("No", s, "values in object."),
                paste("Choose a different callTypeLog2R or re-do
                      runSegmentation()"), sep = "\n")
            stop(errorMsg)
        }else if(s == callTypeLog2R){
            cnValue <- slot(object, slots[s])

            calls <- apply(cnValue, 2,
                function(x, callThresh){
                    .callThreshValidity(x, callThresh)
                }, callThresh=callThreshLog2R)
        }
    }

    message("Updating calls in the CNAclinicData object output ...")

    CNAclinic::calls(object) <- as.data.frame(calls)

    message("Done.")

    return(object)

})
