
#' plotSampleData
#'
#' plotSampleData XXXXX
#'
#' @param object CNAclinicData object
#'
#' @return a list containing ggplot objects (1 per sample, list named by sampleNames)
#' @export plotSampleData
#'
#' @importFrom data.table as.data.table melt
#' @importFrom dplyr group_by mutate %>%
#' @import ggplot2 GenomicRanges chemometrics
#'
#' @examples
#'      \dontrun{
#'       vignette("CNAclinic")
#'      }
#'

setMethod("plotSampleData",
    signature(object="CNAclinicData"),
    function(object,
        geneInfo=NULL, clusterGenes=FALSE,
        chromosomesFilter=c("X", "Y", "M", "MT"),
        outlierLog2R=5,
        showLog2RPlot=TRUE,
        showHeatmap=FALSE,
        showDataPoints=TRUE,
        showSegments=TRUE,
        segmentType=c("summary"),
        segHMMColor="green4", segCBSColor="green",
        segPLSColor="purple", segLACBSColor="hotpink",
        segSummaryColor="red2",
        segmentLineWidth=rep(1, length(segmentType)),
        segmentLineType=rep(1, length(segmentType)),
        segmentLegend=showSegments,
        legendSize=12, legendKeySize=2, segmentLegendPosition="bottom",
        segmentLegendHoriz=TRUE,
        showCalls=TRUE, callTypeLog2R="summary",
        callThreshLog2R=c(-0.15, 0.15),
        callLegend=FALSE,
        callLegendPosition="bottom", callLegendHoriz=TRUE,
        pointShape=".", pointSize=1,
        lossColor="orange2",
        gainColor="blue",
        pointColor="black",
        neutralColor="white", missingColor="grey20",
        xlab=ifelse(is.null(geneInfo), "Chromosome", "Gene"), xaxt="s",
        xaxSize=7,
        ylab=NULL, yaxSize=10,ylim=NULL, yaxp=NULL,
        main=NULL, mainSize=10, ...){


        # Check if the gene information is in the correct format
        .checkGeneInfo(geneInfo)

        chromosomesFilter <- as.character(chromosomesFilter)
        chromosomesFilter <- match.arg(chromosomesFilter,
            choices = c(1:22, "X", "Y", "M", "MT", NA),
            several.ok = TRUE)
        segmentType <- match.arg(segmentType,
            choices = c("CBS", "LACBS", "HMM", "PLS", "summary"),
            several.ok = TRUE)
        callTypeLog2R <- match.arg(callTypeLog2R,
            choices = c("CBS", "LACBS", "HMM", "PLS", "summary", "calls"),
            several.ok = FALSE)

        if(length(segmentLineWidth) != length(segmentType))
            segmentLineWidth <- rep(
                segmentLineWidth, length(segmentType))[1:length(segmentType)]

        if(length(segmentLineType) != length(segmentType))
            segmentLineType <- rep(
                segmentLineType, length(segmentType))[1:length(segmentType)]

        names(segmentLineWidth) <- segmentType
        names(segmentLineType) <- segmentType

        #outlierSD <- abs(outlierSD)

        emptySlots <- emptySlots(object)
        if(any(c("bins", "copyNumber") %in% emptySlots))
            stop("object does not contain genomic bin and copy number data")

        if(showSegments){
            slots <- c(
                CBS="segCBS",
                LACBS="segLACBS",
                HMM="segHMM",
                PLS="segPLS",
                summary="segSummary")

            userText <- slots[segmentType]
            userText <- userText[userText %in% filledSlots(object)]

            if(length(names(slots)[slots %in% userText])==0){
                stop(paste0("Provided segmentType: ",
                paste(names(slots)[!(slots %in% slots[segmentType])], collapse=", "),
                " not present in object"))
            }else{
                segmentType <- names(slots)[slots %in% userText]
            }
        }



        condition <- usebin(object)
        copyNumber <- copyNumber(object)
        condition[chromosomes(object) %in% chromosomesFilter] <- FALSE
        calls=NULL

        if(showCalls){

            slots <- c(
                CBS="segCBS",
                LACBS="segLACBS",
                HMM="segHMM",
                PLS="segPLS",
                summary="segSummary",
                calls="calls")

            userText <- slots[callTypeLog2R]
            userText <- userText[userText %in% filledSlots(object)]

            if(length(callTypeLog2R)==0){
                stop(paste0("Provided callTypeLog2R: ",
                            paste(names(slots)[!(slots %in% slots[callTypeLog2R])], collapse=", "),
                            " not present in object. Run callData()"))
            }else{
                callTypeLog2R <- names(slots)[slots %in% userText]
            }

            for(s in names(slots)){
                if((s == callTypeLog2R) & (slots[s] %in% emptySlots)){
                    errorMsg<- paste(
                        paste("No", s, "values in object."),
                        paste("Choose a different callTypeLog2R or use
                        runSegmentation()/callData() functions"), sep = "\n")
                    stop(errorMsg)
                }else if(s == callTypeLog2R){
                    cnValue <- slot(object, slots[s])
                    #cnValue[!condition, ] = NA
                    calls <- apply(cnValue, 2,
                        function(x, callThresh){
                            .callThreshValidity(x, callThresh)
                        }, callThresh=callThreshLog2R)
                }
            }
        }


        if(!is.null(geneInfo) & showLog2RPlot){
            message("Gene information only shown as heatmap. Setting showLog2RPlot=FALSE")
            showLog2RPlot <- FALSE
            showHeatmap <- TRUE
        }
        if(showHeatmap & showLog2RPlot){
            message("showHeatmap and showLog2RPlot both set as TRUE. Setting showLog2RPlot=FALSE")
            showLog2RPlot <- FALSE

        }else if(!showLog2RPlot & !showHeatmap){
            message('showLog2RPlot & showHeatmap are both FALSE. Setting showLog2RPlot=TRUE')
            showLog2RPlot <- TRUE
        }

        #---------------------------------------------------------
        # Plot read counts and copy number segments and calls
        #---------------------------------------------------------


        if (is.null(ylim))
            ylim <- c(-3, 3)
        if (is.null(ylab))
            ylab <- expression(log[2]~Ratio)
        if (is.null(yaxp))
            yaxp <- c(ylim[1], ylim[2], ylim[2]-ylim[1])
        if (is.null(main))
            main <- CNAclinic::sampleNames(object)
        if (length(ylab) == 1)
            ylab <- rep(ylab, times=ncol(copyNumber(object)))
        baseLine <- 0

        plotList <- list()

        for(i in seq_len(ncol(copyNumber))){

            condition <- usebin(object)
            condition[chromosomes(object) %in% chromosomesFilter] <- FALSE

            cn <- copyNumber[, i]

            # Remove outlier values with user supplied cutoff
            condition[abs(cn) >= outlierLog2R] <- FALSE
            cn[abs(cn) >= outlierLog2R] <- NA
            cn <- cn[condition]

            # Remove outliers through trimmed SD
            #trimmedSD <- chemometrics::sd_trim(na.omit(cn),trim=outlierTrim,
            #const=TRUE)
            #condition[abs(cn) >= (outlierSD * trimmedSD)] <- FALSE

            message("Plotting sample ", main[i],
                " (", i, " of ", ncol(copyNumber), ") ...",
                appendLF=FALSE)

            pointcol <- rep(pointColor, length(na.omit(cn)))

            if(!is.null(calls)){
                call <- calls[condition, i]
                pointcol[call == 0] <- pointColor
                pointcol[call == -1] <- lossColor
                pointcol[call == 1] <- gainColor
            }

            bpstart <- bpstart(object)[condition]
            bpend <- bpend(object)[condition]

            all.chrom <- chromosomes(object)
            all.chrom[all.chrom == "X"] <- 23
            all.chrom[all.chrom == "Y"] <- 24
            all.chrom[all.chrom == "MT"] <- 25

            chrom <- all.chrom[condition]
            uni.chrom <- unique(chrom)
            chrom.num <- as.integer(factor(chrom, levels=uni.chrom, ordered=TRUE))
            uni.chrom.num <- unique(chrom.num)

            pos <- pos2 <- 1:sum(condition)
            chrom.ends <- aggregate(pos, by=list(chromosome=chrom), FUN=max)$x


            if(showLog2RPlot){



                # Get actual positions (Q1, median, Q3) per chromosome

                if(length(uni.chrom.num) <= 2){

                    which_quantiles <- c(2, 3, 4)
                    xlabels2 <- format(bpstart, digits=1, scientific=TRUE)
                    xlabels2pos <-
                        round(as.vector(
                            aggregate(pos, by=list(chromosome=as.integer(chrom)),
                                      FUN=function(x){
                                          quantile(x, type=1)[which_quantiles]
                                      })$x))
                }else{

                    xlabels2 = rep(" ", length(bpstart))

                }
                names(xlabels2) <- pos

                plotDF <- data.frame(chrom, pos, bpstart, cn, pointcol,
                    stringsAsFactors=FALSE)

                p1 <- ggplot(plotDF, aes(x=pos, y=cn)) + theme_bw()
                p1 <- p1 +
                    scale_y_continuous(ylab, limits=ylim) +
                    ggtitle(main[i]) +
                    theme(plot.title = element_text(face="bold", size=mainSize, hjust = 0.5),
                        axis.text.x=element_text(size=xaxSize),
                        panel.grid.major.y = element_line(colour="grey90"),
                        panel.grid.minor.y = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank())


                if(showDataPoints)
                    p1 <- p1 +
                        geom_point(
                            size=pointSize,
                            color=pointcol,
                            shape=pointShape)

                rm(list=c('plotDF')); gc(FALSE)

                if(!is.na(xaxt) && xaxt != "n") {
                    chrom.ends <- chrom.ends[order(chrom.ends)]
                    ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)]))/2

                    if(length(uni.chrom.num) <= 2){
                        ax <- c(ax, xlabels2pos)
                        uni.chrom <- c(paste("\n\n", uni.chrom, sep = ""),
                        xlabels2[xlabels2pos])
                    }

                    uni.chrom[uni.chrom == 23] <- "X"
                    uni.chrom[uni.chrom == 24] <- "Y"
                    uni.chrom[uni.chrom == 25] <- "MT"

                    p1 <- p1 +
                        scale_x_continuous(
                            name=xlab,
                            breaks=ax,
                            labels=uni.chrom)

                }

                chrom.ends <- sort(chrom.ends)

                xintercept <- chrom.ends[!(seq_along(chrom.ends) %% 2)]

               # geom_vline(xintercept=chrom.ends+0.5, colour="black") +

                p1 <- p1 +
                    geom_vline(
                        xintercept = c(0, chrom.ends),
                        colour = "grey90",
                        linetype = 1) +
                    geom_hline(yintercept = baseLine, colour = "black")
            }

            if(showSegments & showLog2RPlot){
                legendCol = c()
                legendLab = c()
                segmented <- NULL
                gc(FALSE)

                for(s in 1:length(segmentType)){
                    segColor <- "white"
                    segmented <- NULL
                    if("CBS" %in% segmentType[s]){
                        segmented <- segCBS(object)[condition, i]
                        segColor <- segCBSColor
                        legendLab = c(legendLab, segmentType[s])
                    }else if("LACBS" %in% segmentType[s]){
                        segmented <- segLACBS(object)[condition, i]
                        segColor <- segLACBSColor
                        legendLab = c(legendLab, segmentType[s])
                    }else if("HMM" %in% segmentType[s]){
                        segmented <- segHMM(object)[condition, i]
                        segColor <- segHMMColor
                        legendLab = c(legendLab, segmentType[s])
                    }else if("PLS" %in% segmentType[s]){
                        segmented <- segPLS(object)[condition, i]
                        segColor <- segPLSColor
                        legendLab = c(legendLab, segmentType[s])
                    }else if("summary" %in% segmentType[s]){
                        segmented <- segSummary(object)[condition, i]
                        segColor <- segSummaryColor
                        legendLab = c(legendLab, segmentType[s])
                    }

                    if(showLog2RPlot){
                        tempDF <- data.frame(chrom, segmented, stringsAsFactors=FALSE)
                        # Note: Should NAs be removed from tempDF?
                        # Should we threshold using outlierLog2R here as well?

                        if(nrow(tempDF) == 0) next
                        flag <- TRUE

                        tryCatch({
                            segment <- .makeSegments(tempDF$segmented, tempDF$chrom)

                        }, error=function(e) {
                            warning(paste("An error has occured in .makeSegments()",
                            paste(segmentType[s],
                            "segments will not be plotted."), sep="\n"))
                            flag <- FALSE
                        })
                        if (!flag){
                            legendLab <- legendLab[-c(length(legendLab))]
                            next
                        }

                        p1 <- p1 +
                            geom_segment(x=pos[segment$start], y=segment$values,
                                xend=segment$end, yend=segment$values,
                                color=segColor,
                                linetype=segmentLineType[s],
                                size=segmentLineWidth[s],
                                data=as.data.frame(segment, stringsAsFactors=FALSE))

                        legendCol <- c(legendCol, segColor)
                    }
                }

                rm(list=c('tempDF')); gc(FALSE)

                if(segmentLegend & (length(legendCol) > 0) & showLog2RPlot){

                    fakeDF <-data.frame(x=rep(rep(0,20), length(legendLab)),
                                     y=rep(rep(0,20), length(legendLab)),
                                     legendLab2=rep(legendLab, 20))

                   # suppressMessages(suppressWarnings(

                    p1 <- p1 +
                        geom_line(data=fakeDF, aes(x, y, color=legendLab2)) +
                        scale_colour_manual(
                            name="Segments: ",
                            values=setNames(legendCol, legendLab),
                            guide='legend') +
                        guides(colour=guide_legend(override.aes=
                            list(linetype=segmentLineType[legendLab]))) +
                        theme(legend.position=segmentLegendPosition,
                            legend.box=ifelse(segmentLegendHoriz,
                            "horizontal", "vertical"),
                            legend.key = element_blank(),
                            legend.key.size=unit(legendKeySize, "lines"),
                            legend.text=element_text(size=legendSize),
                            legend.title =element_text(size=legendSize),
                            axis.title=element_text(size=12))
                   # ))
                }
            }

            if(showHeatmap){

                heatmapDF <- NULL
                segmentTypeHeat <- filledSlots(object)
                segmentTypeHeat <-
                    segmentTypeHeat[!(segmentTypeHeat %in%
                        c("bins", "copyNumber", "calls"))]

                segmentTypeColumns <- c()
                for(s in 1:length(segmentTypeHeat)){
                    segColor <- "white"
                    segmented <- NULL

                    if(grepl("LACBS", segmentTypeHeat[s], ignore.case = TRUE)){
                        segmented <- segLACBS(object)[condition, i]
                        segmentTypeColumns <- c(segmentTypeColumns, "LACBS")
                    }else if(grepl("CBS", segmentTypeHeat[s], ignore.case = TRUE)){
                        segmented <- segCBS(object)[condition, i]
                        segmentTypeColumns <- c(segmentTypeColumns, "CBS")
                    }else if(grepl("HMM", segmentTypeHeat[s], ignore.case = TRUE)){
                        segmented <- segHMM(object)[condition, i]
                        segmentTypeColumns <- c(segmentTypeColumns, "HMM")
                    }else if(grepl("PLS", segmentTypeHeat[s], ignore.case = TRUE)){
                        segmented <- segPLS(object)[condition, i]
                        segmentTypeColumns <- c(segmentTypeColumns, "PLS")
                    }else if(grepl("summary", segmentTypeHeat[s], ignore.case = TRUE)){
                        segmented <- segSummary(object)[condition, i]
                        segmentTypeColumns <- c(segmentTypeColumns, "summary")
                    }

                    if(showHeatmap){
                        if(s == 1)
                            heatmapDF <- data.frame(segmented,
                                stringsAsFactors = FALSE)
                        else
                            heatmapDF <- cbind(heatmapDF, segmented)
                    }
                }

                names(heatmapDF) = segmentTypeColumns

                # Should be no NAs due to filtered bins in matrix.
                #heatmapDF[abs(heatmapDF) >= outlierLog2R]  <- NA


                if(showCalls){
                    heatmapDF <- apply(heatmapDF, 2,
                                function(x, callThresh){
                                    .callThreshValidity(x, callThresh)
                                }, callThresh=callThreshLog2R)
                }

                if(!is.null(geneInfo)){

                    featureName <- rep("", length(bpstart))
                    tempDF <- cbind(chrom, bpstart, bpend, featureName, heatmapDF)


                    # Create 1-based GRanges objects
                    gene_tx_GRange <- with(geneInfo,
                        GenomicRanges::GRanges(chromosome, IRanges::IRanges(start, end, names = geneSymbol)))


                    # Get the overlap of genes with copy number segment values for each algorithm/sample

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

                    names(geneSegmentCalls_df) = segmentTypeColumns
                    heatmapDF <- cbind(
                        chrom=as.character(GenomicRanges::seqnames(gene_tx_GRange)),
                        bpstart=GenomicRanges::start(gene_tx_GRange),
                        bpend=GenomicRanges::end(gene_tx_GRange),
                        featureName=names(gene_tx_GRange),
                        geneSegmentCalls_df)

                    colnames(heatmapDF)[5:length(colnames(heatmapDF))] <- segmentTypeColumns

                    rm(list=c('tempDF', 'geneSegmentCalls_df')); gc(FALSE)

                }else{

                    featureName <- rep("", length(bpstart))
                    heatmapDF <- cbind(chrom, bpstart, bpend, featureName, heatmapDF)


                }



                #return(heatmapDF)

                p2 <- .plotHeatmap(
                    x=heatmapDF,
                    type=ifelse(is.null(geneInfo), "chrom", "gene"),
                    clusterGenes=clusterGenes,
                    clusterSamples=FALSE,
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

                if(!showLog2RPlot)
                    p2 <- p2 + ggtitle(main[i]) +
                        theme(plot.title = element_text(face="bold",
                            size=mainSize, hjust = 0.5))

            }

            if(callLegend)
                print("The call type legend is not currently plotted.")


            if(showLog2RPlot & !showHeatmap)
                plotList[[sampleNames(object)[i]]] <- p1
            else if(!showLog2RPlot & showHeatmap)
                plotList[[sampleNames(object)[i]]] <- p2

        }

        return(plotList)
    })




.plotHeatmap = function(x, type, clusterGenes, clusterSamples, showCalls, main, mainSize, xlab, xaxt,
                        xaxSize,
                        yaxSize,
                        legendSize,
                        legendKeySize,
                        lossColor,
                        gainColor,
                        neutralColor,
                        missingColor){

    if(type != "gene"){
        x[x[ ,"chrom"]=="X", "chrom"] <- 23
        x[x[ ,"chrom"]=="Y", "chrom"] <- 24
        x[x[ ,"chrom"]=="MT", "chrom"] <- 25
        x[ ,"chrom"]<- as.numeric(x[ ,"chrom"])

        # x is now a numeric matrix but featureName has gone from empty strings to NAs
        x <- apply(x, 2, as.numeric)
        # Order by chromosome & bpstart
        x <- x[order(x[,"chrom"],x[,"bpstart"],decreasing=FALSE),]
    }

    xDF_T <- t(x[,-c(1:4)])
    xDF_T <- cbind(rownames(xDF_T), xDF_T)
    xDF_T <- data.frame(xDF_T, stringsAsFactors=FALSE)

    if(type == "gene"){

        colnames(xDF_T) <- c("groupVar", x[ , "featureName"])

    }else{

        colnames(xDF_T) <- c("groupVar", 1:(ncol(xDF_T)-1))
    }

    # Get the order for clustering rows and columns
    # Convert NA values to numeric
    xDF_T[is.na(xDF_T)] <- -999

    column_order <- hclust(dist(t(xDF_T[ ,-1]), method = "euclidean"), method = "complete")$order

    row_order <- hclust(dist(xDF_T[ ,-1], method = "euclidean"),
                        method = "complete")$order

    # Re-insert NA values
    xDF_T[xDF_T == -999] <- NA

    rownames(xDF_T) <- NULL

    xDF_TM <- data.table::melt(data.table::as.data.table(xDF_T),
        id = "groupVar", variable.factor=FALSE)
    xDF_TM <- as.data.frame(xDF_TM)
    xDF_TM$value <- as.numeric(xDF_TM$value)

    if(type=="gene" & clusterGenes){

        xDF_TM$variable <- factor(xDF_TM$variable, levels=names(xDF_T)[-1][column_order])

    }else if(type=="gene" & !clusterGenes){

        xDF_TM$variable <- factor(xDF_TM$variable, levels=names(xDF_T))

    }else{

        xDF_TM$variable <- factor(xDF_TM$variable, levels=1:nrow(x))

        if(clusterSamples){
            xDF_TM$groupVar <- factor(xDF_TM$groupVar, levels=xDF_T$groupVar[row_order])
        }
    }

    if(showCalls){

        names(xDF_TM) <- c("groupVar", "variable", "calls")

        xDF_TM$calls[xDF_TM$calls == 1] <- c("GAIN")
        xDF_TM$calls[xDF_TM$calls == 0] <- c("NEUTRAL")
        xDF_TM$calls[xDF_TM$calls == -1] <- c("LOSS")

        p2 <- ggplot2::ggplot(xDF_TM, aes(variable, groupVar)) +
            geom_raster(aes(fill=calls)) +
            scale_fill_manual(values=c(LOSS=lossColor,
                NEUTRAL=neutralColor, GAIN=gainColor),
                na.value=missingColor)
    }else{

        minValue <- min(xDF_TM$value, na.rm=T)
        maxValue <- max(xDF_TM$value, na.rm= T)

        xDF_TM <- xDF_TM %>% group_by(groupVar) %>%
            mutate(log2R=scales::rescale(value,
                to=c(minValue, maxValue)))

        p2 <- ggplot2::ggplot(xDF_TM, aes(variable, groupVar)) +
            geom_raster(aes(fill=log2R)) +
            scale_fill_gradient2(
                low=lossColor, mid=neutralColor, high=gainColor,
                na.value=missingColor,
                midpoint=0)
                #,
                #breaks=c(round(minValue,2), 0, round(maxValue, 2)),
                #labels=as.character(c(round(minValue,2), 0, round(maxValue, 2))))

    }

    base_size <- 9
    p2 = p2 + theme_bw(base_size=base_size) +
        geom_hline(yintercept=c(1:length(unique(xDF_TM$groupVar)))+0.5, alpha=0.5, size=1.5) +
        labs(x = "", y = "") +
        scale_y_discrete(expand = c(0, 0)) +
        theme(legend.position = "bottom", legend.key.width=unit(3,"lines"),
                text = element_text(size=8))

    if(type == "gene") {
        p2 = p2 + xlab(xlab) +
            geom_vline(xintercept=c(0:length(unique(xDF_TM$variable)))+0.5, alpha=0.5, size=1) +
       #     labs(fill="") +
            theme(legend.position = "bottom",
                legend.text=element_text(size=legendSize),
                legend.key.size=unit(legendKeySize, "lines"),
                legend.title=element_text(size=legendSize),
                legend.key=element_rect(colour = "black"),
                text=element_text(size=10),
                axis.text.x=element_text(size=xaxSize, angle=90),
                axis.text.y=element_text(size=yaxSize),
                axis.title=element_text(size=12),
                panel.grid.major = element_blank())
    }else{

        chrom <- x[ ,"chrom"]
        #chrom[chrom == "X"] <- 23
        #chrom[chrom == "Y"] <- 24
        #chrom[chrom == "MT"] <- 25

        uni.chrom <- unique(as.character(chrom))
        chrom.num <- as.integer(factor(chrom, levels=uni.chrom, ordered=TRUE))
        uni.chrom.num <- unique(chrom.num)


        pos <- pos2 <- 1:nrow(x)
        chrom.ends <- aggregate(pos, by=list(chromosome=as.integer(chrom)), FUN=max)$x

        ### CHECK
        chrom.ends <- chrom.ends[order(chrom.ends)]

        ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)]))/2


        if(length(uni.chrom.num) <= 2){

            which_quantiles <- c(2, 3, 4)
            xlabels2 <- format(as.numeric(x[,"bpstart"]), digits=1, scientific=TRUE)
            xlabels2pos <-
                round(as.vector(
                    aggregate(pos, by=list(chromosome=as.integer(chrom)),
                              FUN=function(x){
                                  quantile(x, type=1)[which_quantiles]
                              })$x))
            ### CHECK
            xlabels2pos <- xlabels2pos[order(xlabels2pos)]

        }else{

            xlabels2 = rep(" ", nrow(x))

        }

        names(xlabels2) <- pos

        uni.chrom[uni.chrom == "23"] <- "X"
        uni.chrom[uni.chrom == "24"] <- "Y"
        uni.chrom[uni.chrom == "25"] <- "MT"

        if(length(uni.chrom.num) <= 2){
            ax <- c(ax, xlabels2pos)
            uni.chrom <- c(paste("\n\n", uni.chrom, sep = ""),
                           xlabels2[xlabels2pos])
        }


        p2 <- p2 +
            scale_x_discrete(expand = c(0, 0),
                name=xlab,
                breaks=as.integer(ax),
                labels=uni.chrom) +
                    #factor(uni.chrom, ordered=TRUE))
            geom_vline(xintercept=chrom.ends+0.5, colour="black") +

            labs(fill="") +

            theme(legend.position = "bottom",
                  legend.text=element_text(size=legendSize),
                  legend.key.size=unit(legendKeySize, "lines"),
                  legend.title =element_text(size=legendSize),
                  legend.key=element_rect(colour="black"),
                  text=element_text(size=10),
                  axis.text.x=element_text(size=xaxSize, angle=90),
                  axis.text.y=element_text(size=yaxSize),
                  axis.title=element_text(size=12),
                  panel.grid.major = element_blank())

    }

}

.callThreshValidity = function(cnValue, callThresh){

    msg <- c("callThresh needs to either be one value which is a
             probability between (0, 0.5] or
             two values (1 positive, 1 negative) between (0, +/-0.5]
             to use as log2 value cutoffs for gain and deletion")

    if(all(unique(as.numeric(as.character(cnValue))) %in% c(-1L,0L,1L, NA))){
        calls <- cnValue
    }else if(!any(is.na(callThresh))){
        if(all(callThresh > 0.5) | all(callThresh <= 0)){
            stop(msg)
        }else{

            if(length(callThresh) == 1){
                callThresh <- abs(callThresh)
                calls <- rep(0L, length(cnValue))
                calls[cnValue <= -callThresh] <- -1L
                calls[cnValue >= callThresh] <- 1L

            }else if(!((callThresh[1] > 0) & (callThresh[2] < 0) |
                       ((callThresh[2] > 0) & (callThresh[1] < 0)))){

                stop(msg)

            }else{
                callThresh <- callThresh[1:2]

                if(!all(abs(callThresh) <= 0.5))
                    stop(msg)

                calls <- rep(0L, length(cnValue))
                calls[cnValue <= callThresh[callThresh < 0]] = -1L
                calls[cnValue >= callThresh[callThresh > 0]] = 1L

            }
        }
    }else{
        stop("callThresh needs to either be one value which is a
             probability between (0, 0.5] or
             two values (1 positive, 1 negative) between
             +/- (0, 0.5] to use as log2 value cutoffs
             for gain and deletion")
    }

    calls[is.na(cnValue)] <- NA
    return(calls)
}


.makeSegments <- function(data,chrdata) {

    counter <- 0
    allends <- c()
    allstarts <- c()
    values <- c()

    # If factor or numeric; need to convert to character first
    chromosome <- as.character(chrdata)
    # Replace 23 by X:
    chromosome[which(chromosome == "23")] <- "X"
    # Replace 24 by Y
    chromosome[which(chromosome == "24")] <- "Y"
    # Replace 25 by MT
    chromosome[which(chromosome == "25")] <- "MT"
    # Replace M by MT
    chromosome[which(chromosome == "M")] <- "MT"

    chromosome[which(chromosome == "X")] <- "23"
    # Replace 24 by Y
    chromosome[which(chromosome == "Y")] <- "24"
    # Replace 25 by MT
    chromosome[which(chromosome == "MT")] <- "25"

    chrdata <- as.numeric(chromosome)
    rm(list=c("chromosome"))

    for(chrom in unique(chrdata)){
        counter <- counter + 1

        datasubset <- data[chrdata == chrom]
        boundaries <- diff(datasubset)

        if(counter == 1){
            allstarts <- starts <- c(1, which(boundaries != 0)+1)
            allends <- ends <- c(which(boundaries != 0), length(datasubset))

        }else{

            starts <- c(1, which(boundaries != 0)+1)
            ends <- c(which(boundaries != 0), length(datasubset))

            #allends <- c(allends, ends + allstarts[length(allstarts)])
            #allstarts <- c(allstarts, starts + allstarts[length(allstarts)])

            allstarts <- c(allstarts, starts + allends[length(allends)])
            allends <- c(allends, ends + allends[length(allends)])


        }

        values <- c(values, datasubset[starts])

    }

    result  <- data.frame(values, start=allstarts, end=allends, stringsAsFactors=FALSE)

    return(result)
}

.checkGeneInfo <-function(geneInfo=NULL){

    if(!is.null(geneInfo)){

        if(!all(c("geneID", "geneSymbol", "chromosome", "start", "end") %in% names(geneInfo))){
            stop(paste("geneInfo should be a data.frame that contains the following columns:",
                paste("geneID", "geneSymbol", "chromosome", "start", "end", sep=","), sep="\n"))
        }

        geneInfo <- geneInfo[ , c("geneID", "geneSymbol", "chromosome", "start", "end")]
        geneInfo <- na.omit(geneInfo)

        if(nrow(geneInfo) < 1)
            stop("geneInfo data.frame is empty")
    }
}
