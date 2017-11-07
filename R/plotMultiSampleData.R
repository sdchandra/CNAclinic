
#' Plots heat map of a cohort of samples
#'
#' @param object CNAclinicData.
#'
#' @return a ggplot2 object
#' @export plotMultiSampleData
#'
#' @examples
#'      \dontrun{
#'       vignette("CNAclinic")
#'      }
#'
setMethod("plotMultiSampleData",
          signature(object="CNAclinicData"),
          function(object,
                   geneInfo=NULL,
                   clusterGenes=FALSE,
                   clusterSamples=FALSE,
                   dataType="summary",
                   chromosomesFilter=c("X", "Y", "M", "MT"),
                   outlierLog2R=6,
                   showCalls=TRUE,
                   callThreshLog2R=c(-0.15, 0.15),
                   callLegend=FALSE,
                   legendSize=12, legendKeySize=2,
                   legendPosition="bottom",
                   legendHoriz=TRUE,
                   lossColor="orange2",
                   gainColor="blue",
                   pointColor="black",
                   neutralColor="white", missingColor="grey20",
                   xlab=ifelse(is.null(geneInfo), "Chromosome", "Gene"),
                   xaxt="s", xaxSize=10,
                   ylab=NULL, yaxSize=10,
                   main=NULL, mainSize=10){


              # Check if the gene information is in the correct format
              .checkGeneInfo(geneInfo)

              chromosomesFilter <- as.character(chromosomesFilter)
              chromosomesFilter <- match.arg(chromosomesFilter,
                                             choices = c(1:22, "X", "Y", "M", "MT", NA),
                                             several.ok = TRUE)

              dataType <- match.arg(dataType,
                                    choices = c("CBS", "LACBS", "HMM", "PLS", "HAAR", "summary", "calls"),
                                    several.ok = FALSE)

              emptySlots <- emptySlots(object)

              slots <- c(CBS="segCBS",
                         LACBS="segLACBS",
                         HMM="segHMM",
                         PLS="segPLS",
                         HAAR="segHaar",
                         summary="segSummary",
                         calls="calls")

              if(slots[dataType] %in% emptySlots){

                  errorMsg<- paste(
                      paste("No", dataType, "values in object."),
                      paste("Choose a different dataType or use
                            runSegmentation()/callData() functions"), sep = "\n")
                  stop(errorMsg)
              }

              if(dataType == "calls")
                  showCalls=TRUE


              # Get the data values specified by the user
              heatmapDF <- slot(object, slots[dataType])

              condition <- usebin(object)
              condition[chromosomes(object) %in% chromosomesFilter] <- FALSE

              if(showCalls){
                  heatmapDF <- apply(heatmapDF, 2,
                                 function(x, callThresh){
                                     .callThreshValidity(x, callThresh)
                                 }, callThresh=callThreshLog2R)
              }


              if (is.null(main))
                  main <- "Heatmap of all samples"



              heatmapDF <- heatmapDF[condition, ]
              heatmapDF[abs(heatmapDF) >= outlierLog2R] <- NA

              bpstart <- bpstart(object)[condition]
              bpend <- bpend(object)[condition]

              all.chrom <- chromosomes(object)
              all.chrom[all.chrom == "X"] <- 23
              all.chrom[all.chrom == "Y"] <- 24
              all.chrom[all.chrom == "MT"] <- 25

              chrom <- all.chrom[condition]

              if (is.null(ylab))
                  ylab <- expression(log[2]~Ratio)

              if (is.null(main))
                  main <- CNAclinic::sampleNames(object)

              if (length(ylab) != 1)
                ylab <- ylab[1]

              if(!is.null(geneInfo)){

                  featureName <- rep("", length(bpstart))
                  tempDF <- cbind(chrom, bpstart, bpend, featureName, heatmapDF)


                  # Create 1-based GRanges objects
                  gene_tx_GRange <- with(geneInfo,
                                         GenomicRanges::GRanges(chromosome, IRanges::IRanges(start, end, names = geneSymbol)))


                  # Get the overlap of genes with copy number segment values for each algorithm

                  geneSegmentCalls_df <- matrix(NA,
                                                nrow=length(gene_tx_GRange),
                                                ncol=ncol(heatmapDF))

                  for(j in 1:ncol(heatmapDF)){

                      # Create 1-based GRanges objects for the copy number segments and/or calls
                      segmentCalls_GRange <-
                          GenomicRanges::GRanges(paste("chr", tempDF[,"chrom"], sep = ""),
                                                 IRanges::IRanges(bpstart, bpend, names=tempDF[,"featureName"]),
                                                 strand=rep("+", nrow(tempDF)),
                                                 value=heatmapDF[ , j])


                      # Find the overlaps
                      # a) ignoring strand
                      # b) report all segments that overlap gene

                      overlaps <- IRanges::findOverlaps(gene_tx_GRange,
                                                        segmentCalls_GRange, select = "all")

                      if(class(overlaps) == "SortedByQueryHits"){  #

                          avValue <- aggregate(segmentCalls_GRange$value[slot(overlaps, "to")],
                                               by = list(slot(overlaps, "from")), mean)

                          value <- rep(NA, length(gene_tx_GRange))
                          value[avValue$Group.1] <- avValue$x

                          geneSegmentCalls_df[ ,j] <- value


                      }else if(class(overlaps) == "integer"){

                          geneSegmentCalls_df[ ,j] <- segmentCalls_GRange$value[overlaps]

                      }
                  }

                  names(geneSegmentCalls_df) <- sampleNames(object)
                  heatmapDF <- cbind(
                      chrom=as.character(GenomicRanges::seqnames(gene_tx_GRange)),
                      bpstart=GenomicRanges::start(gene_tx_GRange),
                      bpend=GenomicRanges::end(gene_tx_GRange),
                      featureName=names(gene_tx_GRange),
                      geneSegmentCalls_df)

                  colnames(heatmapDF)[5:length(colnames(heatmapDF))] <- sampleNames(object)

                  rm(list=c('tempDF', 'geneSegmentCalls_df')); gc(FALSE)

              }else{

                  featureName <- rep("", length(bpstart))
                  heatmapDF <- cbind(chrom, bpstart, bpend, featureName, heatmapDF)


              }


              p2 <- .plotHeatmap(
                  x=heatmapDF, type=ifelse(is.null(geneInfo), "chrom", "gene"),
                  clusterGenes=clusterGenes,
                  clusterSamples=clusterSamples,
                  showCalls=showCalls,
                  main=main[i],
                  mainSize=mainSize,
                  xlab=xlab,
                  xaxt=xaxt,
                  xaxSize=xaxSize,
                  yaxSize=yaxSize,
                  legendSize=legendSize,
                  legendKeySize=legendKeySize,
                  lossColor=lossColor,
                  gainColor=gainColor,
                  neutralColor=neutralColor,
                  missingColor=missingColor)


                 p2 <- p2 + ggtitle(main) +
                 theme(plot.title = element_text(face="bold",
                                                  size=mainSize, hjust = 0.5))

                 return(p2)

})
# EOF
