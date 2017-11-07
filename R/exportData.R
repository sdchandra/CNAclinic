#' Export data
#'
#' @param x X
#' @param dataType X
#' @param geneInfo X
#' @param fileName X
#' @param fileType X
#' @param filterNA X
#' @param digits X
#' @param ... X
#'
#' @return NULL
#' @export exportData
#' @import IRanges
#'
#' @author Dineika Chandrananda
#' @examples
#'      \dontrun{
#'       vignette("CNAclinic")
#'      }
#'

exportData <- function(x, dataType, geneInfo=NULL,
                       fileName=NULL, fileType=c("csv", "igv"),
                       filterNA=TRUE, digits=3, ...) {


    if(class(x) != "CNAclinicData")
        stop("exportData() needs a CNAclinicData object as input")

    if(is.null(fileName))
        stop("Provide a fileName")

    fileType <- match.arg(fileType)
    dataType <- match.arg(dataType,
        c("copynumber", "summary", "CBS", "LACBS", "HMM", "PLS", "calls"))

    slots <- c(copynumber="copyNumber", CBS="segCBS", LACBS="segLACBS",
               HMM="segHMM", PLS="segPLS", summary="segSummary", calls="calls")
    dataType <- slots[dataType]

    if(!is.null(geneInfo)){

        if(!all(c("geneSymbol", "chromosome", "start", "end") %in% names(geneInfo))){
            stop(paste("geneInfo should be a data.frame that contains the following columns:",
                       paste("geneSymbol", "chromosome", "start", "end", sep=","), sep="\n"))
        }

        geneInfo <- geneInfo[ , c("geneID", "geneSymbol", "chromosome", "start", "end")]
        geneInfo <- na.omit(geneInfo)

        if(nrow(geneInfo) < 1)
            stop("geneInfo data.frame is empty")

    }


    if(dataType %in% CNAclinic::emptySlots(x))
        stop("Provided dataType not present in CNAclinicData object")

    data <- slot(x, dataType)
    # Convert to data.frame if only one sample is present
    if(!is.data.frame(data)){
        data <- data.frame(data)
        colnames(data) <- sampleNames(x)
    }

    if(filterNA) {
        data <- data[CNAclinic::usebin(x), ,drop=FALSE]
        bins <-CNAclinic::bins(x)[CNAclinic::usebin(x), ]
    }else{
        bins <-CNAclinic::bins(x)
    }

    if(is.numeric(digits)) {
        data <- round(data, digits=digits)
    }

    tmp <- options("scipen")
    options(scipen=15)

    if(is.null(geneInfo)){

        if(fileType == "csv") {
            out <- data.frame(bins, data, check.names=FALSE, stringsAsFactors=FALSE)
            write.csv(out, file=fileName,
                quote=FALSE, na="", row.names=FALSE, ...)
        }else if (fileType == "igv") {
            out <- data.frame(bins, data, check.names=FALSE, stringsAsFactors=FALSE)
            cat('#type=COPY_NUMBER\n#track coords=1\n', file=fileName)
            suppressWarnings(write.table(out, file=fileName, append=TRUE,
                quote=FALSE, sep="\t", na="", row.names=FALSE, ...))
        }
    }else{

        #tempDF <- cbind(chromosome, bpstart, bpend, featureName, heatmapDF)

        # Create 1-based GRanges objects
        gene_tx_GRange <- with(geneInfo,
            GenomicRanges::GRanges(chromosome, IRanges::IRanges(start, end, names=geneSymbol)))

        # Get the overlap of genes with copy number segment values for each sample
        geneSegmentCalls_df <- matrix(NA,
            nrow=length(gene_tx_GRange),
            ncol=ncol(data))

        colnames(geneSegmentCalls_df) <- CNAclinic::sampleNames(x)

        for(j in 1:ncol(data)){

            # Create 1-based GRanges objects for the copy number segments and/or calls
            segmentCalls_GRange <-
                GenomicRanges::GRanges(paste("chr", bins[,"chromosome"], sep = ""),
                                       IRanges::IRanges(bins[,"start"], bins[,"end"], names=1:nrow(bins)),
                    strand=rep("+", nrow(bins)),
                    value=data[ , j])

            # Find the overlaps
            # a) ignoring strand
            # b) report all segments that overlap gene

            overlaps <- IRanges::findOverlaps(gene_tx_GRange,
                segmentCalls_GRange, select = "all", ignore.strand=TRUE)

            if(class(overlaps) == "SortedByQueryHits"){

                avValue <- aggregate(segmentCalls_GRange$value[slot(overlaps, "to")],
                                     by = list(slot(overlaps, "from")), mean)

                value <- rep(NA, length(gene_tx_GRange))
                value[avValue$Group.1] <- avValue$x

                geneSegmentCalls_df[ ,j] <- value

            }else if(class(overlaps) == "integer"){

                geneSegmentCalls_df[ ,j] <- segmentCalls_GRange$value[overlaps]

            }
        }

        colnames(geneSegmentCalls_df) <- CNAclinic::sampleNames(x)

        bins <- cbind(chromosome=as.character(GenomicRanges::seqnames(gene_tx_GRange)),
            start=start(gene_tx_GRange),
            end=end(gene_tx_GRange),
            featureName=names(gene_tx_GRange))

        if(fileType == "csv") {
            out <- data.frame(bins, geneSegmentCalls_df,
                              check.names=FALSE,
                              stringsAsFactors=FALSE)
            colnames(out) <- c(colnames(bins), colnames(geneSegmentCalls_df))
            write.csv(out, file=fileName, quote=FALSE, na="",
                      row.names=FALSE, ...)
        }else if (fileType == "igv") {
            out <- data.frame(bins, geneSegmentCalls_df, check.names=FALSE, stringsAsFactors=FALSE)
            colnames(out) <- c(colnames(bins), colnames(geneSegmentCalls_df))
            cat('#type=COPY_NUMBER\n#track coords=1\n', file=fileName)
            suppressWarnings(write.table(out, file=fileName, append=TRUE,
                                         quote=FALSE, sep="\t", na="", row.names=FALSE, ...))
        }

    }

    options(scipen=tmp)
}
