CNAclinicData <- setClass("CNAclinicData",
    slots = c(
        bins="data.frame",
        copyNumber="data.frame",
        segCBS="data.frame",
        segLACBS="data.frame",
        segHMM="data.frame",
        segPLS="data.frame",
        segSummary="data.frame",
        calls="data.frame"))

setValidity("CNAclinicData", function(object) {
    errorMsg <- NULL
    valid <- TRUE

    slotNames <- slotNames(object)

    rowsOfSlots <- "c("
    colsOfSlots <- "c("
    sampleNames <- "c("

    for(s in 1:length(slotNames)){

        rowsOfSlots <-
            paste(rowsOfSlots, "nrow(object@", slotNames[s], ")", sep = "")
        if(slotNames[s] != "bins"){
            colsOfSlots <-
            paste(colsOfSlots, "ncol(object@", slotNames[s], ")", sep = "")
            sampleNames <- paste(sampleNames, "names(object@",
                slotNames[s], ")", sep = "")
        }

        if(s != length(slotNames)){
            rowsOfSlots <- paste(rowsOfSlots, ",", sep = "")
            if(slotNames[s] != "bins"){
                colsOfSlots <- paste(colsOfSlots, ",", sep = "")
                sampleNames <- paste(sampleNames, ",", sep = "")
            }

        }else{
            rowsOfSlots <- paste(rowsOfSlots, ")", sep = "")
            colsOfSlots <- paste(colsOfSlots, ")", sep = "")
            sampleNames <- paste(sampleNames, ")", sep = "")
        }
    }

    rowN <- eval(parse(text=rowsOfSlots))
    colN <- eval(parse(text=colsOfSlots))
    sampleNames <- eval(parse(text=sampleNames))

    if(nrow(object@bins) != 0){
        if(!all(names(object@bins) %in% c("chromosome", "start", "end", "usebin"))){
            valid <- FALSE
            errorMsg <- paste(c(errorMsg,
                        "Column names in data.frame bins should be named:
                        chromosome, start, end, usebin", sep="\n"))
        }
    }

    if(!(length(unique(rowN[rowN != 0])) %in% c(0, 1))){
        valid <- FALSE
        errorMsg <- paste(c(errorMsg,
            paste("Number of rows in", paste(slotNames[rowN[rowN != 0]],
            collapse = ", "),
            "need to be identical")), sep = "\n")
    }
    if(!(length(unique(colN[colN != 0])) %in% c(0,1))){
        valid <- FALSE
        slotNames <- slotNames[slotNames != "bins"]
        errorMsg <- paste(c(errorMsg,
                    paste("Number of columns in",
                    paste(slotNames[colN[colN != 0]], collapse=", "),
                    "need to be identical")), sep="\n")
    }

    if(valid & (length(sampleNames) != 0)){

        if(unique(colN[colN != 0]) != length(unique(sampleNames))){
            valid <- FALSE
            slotNames <- slotNames[slotNames != "bins"]
            errorMsg <- paste(c(errorMsg,
                        paste("Column names (i.e. the sampleNames) in",
                        paste(slotNames[colN[colN != 0] != length(unique(sampleNames))], collapse=", "),
                        "need to be unique for each sample
                        and identical between the slots")), sep="\n")
        }

    }

    if (valid) TRUE else errorMsg
})

setMethod("show",
    signature=c(object="CNAclinicData"),
    definition=function(object){
        cat("An object of class ", class(object), " containing multiple data.frames\n", sep = "")

        cat(" ", "(i) bins is a data.frame containing ",
          nrow(object@bins),
          " genomic bins with the following columns:\n\t ",
          paste(names(object@bins), collapse = ", "), "\n", sep = "")

        cat(" ", "(ii) copyNumber is a data.frame containing ",
          nrow(object@copyNumber), " copy number measurements for ",
          ncol(object@copyNumber), " samples.\n\t", sep = "")

        if(!is.null(nrow(object@copyNumber))){
          cat(" ", "sample names: ",
              paste(sampleNames(object), collapse = ", "), "\n", sep = "")
        }

        cat(" ", "(iii) segCBS is a data.frame containing ",
          nrow(object@segCBS), " CBS segments for ",
          ncol(object@segCBS), " samples.\n", sep = "")

        cat(" ", "(iv) segLACBS is a data.frame containing ",
            nrow(object@segLACBS), " LACBS segments for ",
            ncol(object@segLACBS), " samples.\n", sep = "")

        cat(" ", "(v) segHMM is a data.frame containing ",
          nrow(object@segHMM), " HMM segments for ",
          ncol(object@segHMM), " samples.\n", sep = "")

        cat(" ", "(vi) segPLS is a data.frame containing ",
          nrow(object@segPLS), " PLS segments for ",
          ncol(object@segPLS), " samples.\n", sep = "")

        cat(" ", "(viii) segSummary is a data.frame containing ",
          nrow(object@segSummary), " summarised segment values for ",
          ncol(object@segSummary), " samples.\n", sep = "")

        cat(" ", "(ix) calls is a data.frame containing ",
          nrow(object@calls), " copy number calls for ",
          ncol(object@calls), " samples.\n", sep = "")

        invisible(NULL)})


setMethod("summary",
    signature=c(object="CNAclinicData"),
    definition=function(object, ...){
        show(object)
        invisible(NULL)
})

# EOF
