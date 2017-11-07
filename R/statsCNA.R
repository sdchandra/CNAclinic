#' Quantify and quality check copy-number profiles.
#'
#' Calculate the fraction of genome altered (FGA) and/or the
#' median of the absolute values of all pairwise differences (MAPD) between
#' adjacent genomic bin values.
#'
#' @param object a CNAclinicData object containing copy number calls
#' @param measure Either "FGA", "MAPD" or both
#'
#' @return a matrix containing the values for the chosen measure(s) for each
#' sample
#'
#' @export statsCNA
#'
#' @examples
#'      \dontrun{
#'       vignette("CNAclinic")
#'      }
#'


setMethod('statsCNA', signature(object="CNAclinicData"),
    definition=function(object, measure=c("FGA", "MAPD"),
    digits=3) {

    measure <- match.arg(measure, choices=c("FGA", "MAPD"),
                                   several.ok=TRUE)

    #measure <- na.omit(c(measure, ifelse("MAPD" %in% measure, "SD", NA)))

    stats <- matrix(NA, nrow=length(sampleNames(object)), ncol=length(measure))
    colnames(stats) <- measure
    row.names(stats) <- NULL

    for(i in 1:length(measure)){

        if(measure[i] == "FGA"){
            dataValues <- CNAclinic::calls(object)
            if("calls" %in% emptySlots(object)){
                warning("CNAclinicData object does not contain CNA calls")
                next()
            }
        }else if(measure[i] == "MAPD"){
            dataValues <- CNAclinic::copyNumber(object)
            if("copyNumber" %in% emptySlots(object)){
                warning("CNAclinicData object does not contain log2 ratios")
                next()

            }
        }

        stats[ , i] <- apply(dataValues, 2, function(x, stat=measure[i]){

            x <- na.omit(x)

            if(stat == "FGA"){
                round(sum(x != 0)/length(x), digits=digits)
            }else if(stat == "MAPD"){
                round(median(abs(diff(x))), digits=digits)
            }else if(stat == "SD"){
                round(sd(x), digits=digits)
            }
            })
    }

    stats <- cbind(sampleNames=sampleNames(object),
        as.data.frame(stats, stringsAsFactors=FALSE) , stringsAsFactors=FALSE)

    return(stats)
})

# EOF
