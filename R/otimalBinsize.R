#' Assess optimal genomic bin size to partition read counts.
#'
#' Calculate Akaike's information criterion (AIC) and cross-validation (CV)
#' log-likelihood to infer the optimal bin size to partition read counts across
#' genome.
#'
#' @details As a guidance, choose bin sizes which have low AIC and/or high CV values
#' but also contain 30-180 read counts on average. This strikes a reasonable
#' balance between error variability and bias of CNA. Using a much smaller bin
#' size may result in many genomic regions with zero read count and make the
#' overall analysis non-informative. At the other extreme, using a much bigger
#' bin size will 'smooth out' some pattern of alteration (i.e. increasing bias).
#' The process of estimating the optimal bin size is in the context of
#' low-coverage sequence data, so use sensible values for the binSizes argument
#' when the input data is not of shallow whole-genome depth (<10 million reads).
#'
#' @param bamfiles A \code{\link[base]{character}} vector of BAM file names
#' with or without full path. If NULL (default), all files with extension
#' .bam, are read from directory path.
#' @param bamnames An optional \code{\link[base]{character}} vector of sample
#' names. Defaults to file names with extension \code{.bam} removed.
#' @param pathToBams If \code{bamfiles} is NULL, all files ending with ".bam"
#' extension will be read from this path.
#' @param binSizes A \code{\link[base]{numeric}} vector of genomic bin sizes,
#' in units of kilo base pairs (1000 base pairs),
#' e.g. \code{binSizes = c(10, 30, 50)} corresponds to bins of 10, 30 and 50 kbp
#' bins.
#' @param measure The goodness of fit criteria (AIC or CV). Defaults to "CV".
#' @param lineColor Line color to use in plot.
#' @param chromosomesFilter A \code{\link[base]{character}} vector specifying
#' which chromosomes to filter out.
#' Defaults to the sex chromosomes and mitochondrial reads,
#' i.e. \code{c("X", "Y", "M", "MT")}. Use NA to use all chromosomes.
#' @param savePlot if TRUE (default) saves plots of each sample to
#' working directory.
#' @param plotPrefix Prefix for plot title and pdf file name. Defaults to "optimalBinsize".
#' @param minMapq If quality scores exists, the minimum quality score required
#' in order to keep a read (20, default).
#' @param isPaired A \code{\link[base]{logical}}(1) indicating whether unpaired
#' (FALSE), paired (TRUE), or any (NA, default) read should be returned.
#' @param isProperPair A \code{\link[base]{logical}}(1) indicating whether
#' improperly paired (FALSE), properly paired (TRUE), or any (NA, default) read
#' should be returned.
#' @param isUnmappedQuery A \code{\link[base]{logical}}(1) indicating whether
#' unmapped (TRUE), mapped (FALSE, default), or any (NA) read should be returned.
#' @param hasUnmappedMate A \code{\link[base]{logical}}(1) indicating whether
#' reads with mapped (FALSE), unmapped (TRUE), or any (NA, default) mate should
#' be returned.
#' @param isMinusStrand A \code{\link[base]{logical}}(1) indicating whether
#' reads aligned to the plus (FALSE), minus (TRUE), or any (NA, default) strand
#' should be returned.
#' @param isMateMinusStrand A \code{\link[base]{logical}}(1) indicating whether
#' mate reads aligned to the plus (FALSE), minus (TRUE), or any (NA, default)
#' strand should be returned.
#' @param isFirstMateRead A \code{\link[base]{logical}}(1) indicating whether
#' the first mate read should be returned (TRUE) or not (FALSE), or whether mate
#' read number should be ignored (NA, default).
#' @param isSecondMateRead A \code{\link[base]{logical}}(1) indicating whether
#' the second mate read should be returned (TRUE) or not (FALSE), or whether
#' mate read number should be ignored (NA, default).
#' @param isSecondaryAlignment A \code{\link[base]{logical}}(1) indicating
#' whether alignments that are primary (FALSE), are not primary (TRUE) or whose
#' primary status does not matter (NA, default) should be returned.
#' @param isDuplicate A \code{\link[base]{logical}}(1) indicating that
#' un-duplicated (FALSE, default), duplicated (TRUE), or any (NA) reads should
#' be returned.
#'
#' @return Returns a list. The first element is a data.frame holding information of the
#' average read counts per bin size, the other elements are sample-specific
#' \code{\link[ggplot2]{ggplot}} objects.
#'
#' @author Dineika Chandrananda
#'
#' @seealso Internally, the function \code{\link[NGSoptwin]{opt.win.onesample}}
#' of the \pkg{NGSoptwin} package is used.
#'
#' @import Rsamtools
#' @export optimalBinsize
#'
#' @examples
#'      \dontrun{
#'       vignette("CNAclinic")
#'      }
#'

optimalBinsize = function(
    bamfiles=NULL, bamnames=NULL, pathToBams=NULL,
    binSizes=c(10, 30, 50, 100, 250, 500, 750, 1000),
    measure="CV", lineColor="red4",
    chromosomesFilter=c("X", "Y", "M", "MT"),
    savePlot=FALSE, plotPrefix="optimalBinsize",
    minMapq=20, isPaired=NA, isProperPair=NA,
    isUnmappedQuery=FALSE, hasUnmappedMate=NA,
    isMinusStrand=NA, isMateMinusStrand=NA,
    isFirstMateRead=NA, isSecondMateRead=NA,
    isSecondaryAlignment=NA,
    isDuplicate=FALSE){

    scipen <- options("scipen")
    options(scipen=100000)

    if (is.null(bamfiles))
        bamfiles <- list.files(ifelse(is.null(pathToBams), '.', pathToBams),
                               pattern=sprintf('%s$', 'bam'), full.names=TRUE)
    if (length(bamfiles) == 0L)
        stop('No files to process.')
    if (is.null(bamnames)) {
        bamnames <- basename(bamfiles)
        bamnames <- sub(sprintf('[\\.]?%s$', 'bam'), '', bamnames)
    } else if (length(bamfiles) != length(bamnames)) {
        stop('bamfiles and bamnames have to be of same length.')
    }

    chromosomesFilter <- as.character(chromosomesFilter)
    chromosomesFilter <- match.arg(chromosomesFilter,
        choices=c(1:22, "X", "Y", "M", "MT", NA), several.ok=TRUE)

    measure <- match.arg(measure, choices=c("AIC", "CV"),
        several.ok=FALSE)

    if(length(binSizes) < 3)
        stop("The length of binSizes vector is less than 3")

    binSizes <- binSizes * 1000

    bamfiles <- normalizePath(bamfiles)
    fullname <- sub('\\.[^.]*$', '', basename(bamfiles))

    optwinList <- list()

    QDNAseq:::vmsg('counting reads ...', appendLF=FALSE)

    for (i in seq_along(bamfiles)) {
        QDNAseq:::vmsg("    ", bamnames[i], " (", i, " of ", length(bamfiles), "): \n",
             appendLF=FALSE)

        flag <- Rsamtools::scanBamFlag(isPaired=isPaired,
            isProperPair=isProperPair, isUnmappedQuery=isUnmappedQuery,
            hasUnmappedMate=hasUnmappedMate, isMinusStrand=isMinusStrand,
            isMateMinusStrand=isMateMinusStrand,
            isFirstMateRead=isFirstMateRead,
            isSecondMateRead=isSecondMateRead,
            isSecondaryAlignment=isSecondaryAlignment,
            isNotPassingQualityControls=FALSE,
            isDuplicate=isDuplicate)

        params <- Rsamtools::ScanBamParam(flag=flag, what=c('rname', 'pos', 'mapq'))
        reads <- Rsamtools::scanBam(bamfiles[i], param=params)
        reads <- reads[[1L]]

        hasMapq <- any(is.finite(reads[['mapq']]))
        if (hasMapq) {
            keep <- which(reads[['mapq']] >= minMapq)
            if (length(keep) == 0L)
                next
            reads <- lapply(reads, FUN=function(x) x[keep])
        }
        reads[['mapq']] <- NULL

        chrs <- unique(reads[['rname']])
        chrs <- sub('^chr', '', chrs)
        chrs <- chrs[(chrs %in% c(1:22, "X", "Y", "M", "MT"))]
        chrs <- chrs[!(chrs %in% chromosomesFilter)]
        reads[['rname']] <- sub('^chr', '', reads[['rname']])

        readDF <- data.frame()
        c <- 0
        for(chr in chrs) {
            message(paste("chr", chr, sep=""))
            c <- c + 1
            keep <- which(reads[['rname']] == chr)

            if(length(keep) == 0L)
                next

            if(c == 1){
                readDF <- data.frame(
                    paste('chr', reads[['rname']][keep], sep=""),
                    reads[['pos']][keep], stringsAsFactors=FALSE)
                names(readDF) <- c("chrom", "pos")
            }else{
                readDF <- rbind(readDF,
                    data.frame(chrom=paste('chr', reads[['rname']][keep], sep=""),
                        pos=reads[['pos']][keep], stringsAsFactors=FALSE))
            }
        }

        rm(list=c('reads')); gc(FALSE)
        names(readDF) <- c("chrom", "pos")

        optwin <- opt.win.onesample(readDF, win.size=binSizes)

        optwinList[[i+1]] <- .plotPerSampOptwin(optwin,
            measure=measure,
            main=bamnames[i],
            colorAIC=lineColor,
            colorCV=lineColor)

        if(i == 1){

            avReadsPerBin <- (attr(optwin, "winstat"))$rpw
            medianReadsPerBin <- (attr(optwin, "winstat"))$rpw

            allSampDF <- data.frame(
                binSizes=(attr(optwin, "winsize"))/1000,
                Value=round((attr(optwin, "winstat"))$rpw),
                stringsAsFactors=FALSE)
            names(allSampDF) <- c("binSizes", bamnames[i])


        }else{

            avReadsPerBin <- avReadsPerBin + (attr(optwin, "winstat"))$rpw
            medianReadsPerBin <- rbind(medianReadsPerBin, (attr(optwin, "winstat"))$rpw)

            tempNames <- names(allSampDF)
            allSampDF <- cbind(allSampDF,
                data.frame(Value=round((attr(optwin, "winstat"))$rpw),
                stringsAsFactors=FALSE))
            names(allSampDF) <- c(tempNames, bamnames[i])
        }

        if(i == length(bamfiles)){

            avReadsPerBin <- round(avReadsPerBin/length(bamfiles), 2)
            if(class(medianReadsPerBin) == "numeric"){
                medianReadsPerBin <- round(median(medianReadsPerBin), 2)
            }else{
                medianReadsPerBin <- round(apply(medianReadsPerBin, 2, median), 2)
            }

            allSampDF <- cbind(allSampDF, medianReadsPerBin, avReadsPerBin)

        }
    }

    if(savePlot){
        pdf(paste(plotPrefix, "_", measure, ".pdf", sep=""),
            onefile = TRUE)
            invisible(lapply(optwinList, print))
        dev.off()
        message("Plots saved to working directory.")
    }

    names(allSampDF) <- c("binsize", bamnames, "median", "average")
    averageCounts <- .plotAverageCount(allSampDF)

    optwinList[[1]] <- averageCounts

    options(scipen=scipen)
    return(optwinList)
}


.plotPerSampOptwin = function(optwin, measure="AIC", colorAIC="orange2",
                              colorCV="deeppink2", main=""){

    if(measure == "AIC")
        x <- data.frame(binSizes=(attr(optwin, "winsize"))/1000,
            measure=optwin$aic,
            readsPerBin=(attr(optwin, "winstat"))$rpw)
    else
        x <- data.frame(binSizes=(attr(optwin, "winsize"))/1000,
            measure=optwin$cv,
            readsPerBin=(attr(optwin, "winstat"))$rpw)

    x_labels <- as.numeric(unique(x$binSizes))
    x.melted <- reshape2::melt(x, id = c("binSizes", "readsPerBin"))

    names(x.melted) <- c("binSizes", "readsPerBin", "Measure", "Value")
    g <- ggplot(x.melted, aes(binSizes, Value)) + theme_bw() +
            ggtitle(main) +
            geom_line(color=ifelse(measure=="AIC", colorAIC, colorCV), size=1.5) +
            geom_point(shape=19, size=3, color=ifelse(measure=="AIC", colorAIC, colorCV)) +
            scale_x_continuous(name="Bin size in Kbp (logged)
                \n(mean number of reads per bin)",
                trans='log10', breaks=x_labels,
                labels=paste(paste("(", round(x$readsPerBin), ")", sep=""),
                x_labels, sep="  ")) +
            ylab(ifelse(measure=="AIC", "Akaike's information criterion",
                "Cross-validation log-likelihood")) +
        theme(plot.title=element_text(face="bold", hjust = 0.5),
              axis.title=element_text(size = 15),
              axis.text.x=element_text(angle = 90, hjust = 1, size=10),
              axis.text.y=element_text(size=12),
              panel.grid.major= element_line(colour="grey90"),
              panel.grid.minor=element_blank())

    g
}

.plotAverageCount = function(df){

    ncol <- ncol(df)
    meltDF <- reshape::melt(df[,-c(ncol, ncol-1)], "binsize")
    meltDF$binsize <- factor(meltDF$binsize)

    g <- ggplot(meltDF,
                aes(x=binsize, y=value)) +
        geom_point(color="navy") +
        geom_hline(yintercept=c(30, 180), linetype=2, col="red") +
        labs(x="Bin size (Kbp in logged scale)",
             y="Average read count in genomic bins") +
        theme_bw() +
        theme(legend.position="none",
              text = element_text(size=10),
              axis.text.x=element_text(size=12, angle=90),
              axis.text.y=element_text(size=12),
              axis.title = element_text(size = 15))

}

# Relevant code from NGSoptwin for ease of distribution
# Author: Arief Gusnanto with some contribution from Ibrahim Nafisah

# A wrapper for hist() to count reads
# per window, by chromosome
# The attribute "density" contains the associated density

hist.ch = function(x,d){

        list.ch <- unique(x[,1])
        mat.count <- NULL # final matrix for counts
        mat.dens <- NULL # final matrix for density
        for(j in list.ch){
            pos <- as.numeric(x[x[,1]==j,2])# take pos by chr
            max.pos.for.bins <- (max(pos, na.rm=T)%/%d)*d
            pos <- pos[pos<=max.pos.for.bins]
            res.hist <- hist(pos, breaks=seq(0,max.pos.for.bins,by=d), plot=F)
            temp <- data.frame(ch=j,count=res.hist$count)
            mat.count <- rbind(mat.count,temp)
            temp <- data.frame(ch=j,dens=res.hist$density)
            mat.dens <- rbind(mat.dens,temp)
        }
        attr(mat.count,"density") <- mat.dens
        return(mat.count)
}

# function to join two read-count matrices
# x is the 'test' data and y is the normal/control data
join.mat = function(x,y){

        list.ch <- unique(x[,1])
        mat <- NULL # final matrix for counts
        for(j in list.ch){
            count.x <- as.numeric(x[x[,1]==j,2])# take count by chr
            count.y <- as.numeric(y[y[,1]==j,2])# take count by chr
            nbin <- min(length(count.x),length(count.y))
            mat <- rbind(mat,data.frame(Chr=j,Test=count.x[1:nbin], Norm=count.y[1:nbin]))
        }
        return(mat)
}

# function to calculate cross-validated log-likelihood
# and AIC across different values of win.size
# INPUT :
# x : reads positions for the test sample
# win.size : a vector of window sizes in UNIT base pair wide
#            e.g. for 10kbp wide, size=10000

opt.win.onesample = function(x, win.size=NULL){

        if(length(win.size)<3) stop("The length of win.size vector is less than 3 (including NULL)")

        test.cv <- c() # final vector of cross-validated log likelihood
        # across different window sizes
        test.aic <- c() # final vector of aic
        # across different window sizes
        test.nwin <- c()
        test.rpw <- c()

        for(k in win.size){
            cat("Window size=",k,"\n")
            m.count.x <- hist.ch(x,k)
            test.nwin <- c(test.nwin, nrow(m.count.x))
            test.rpw <- c(test.rpw, mean(m.count.x[,2], na.rm=T))
            opt.results <- optimal.onesample(m.count.x)
            test.cv <- c(test.cv, opt.results$ML.Test)
            test.aic <- c(test.aic, opt.results$AIC.Test)
        } # end for k win.size

        result <- list(cv=test.cv, aic=test.aic)
        win.stat <- list(nwin=test.nwin, rpw=test.rpw)
        attr(result,"winsize") <- win.size
        attr(result,"winstat") <- win.stat
        attr(result,"class") <- c("optwin.onesample")
        result
}

# function to calculate cross-validated log-likelihood
# and AIC across different values of win.size
# INPUT :
# x : reads positions for the test sample
# y : reads positions for the normal/control sample
# win.size : a vector of window sizes in UNIT base pair wide
#            e.g. for 10kbp wide, size=10000

opt.win = function(x,y, win.size=NULL){


        if(length(win.size)<3) stop("The length of win.size vector is less than 3 (including NULL)")

        test.cv <- c() # final vector of cross-validated log likelihood
        # across different window sizes
        test.aic <- c() # final vector of aic
        # across different window sizes
        control.cv <- c() # final vector of cross-validated log likelihood
        # across different window sizes
        control.aic <- c() # final vector of aic
        # across different window sizes
        test.nwin <- c()
        control.nwin <- c()
        test.rpw <- c()
        control.rpw <- c()

        for(k in win.size){
            cat("Window size=",k,"\n")
            m.count.x <- hist.ch(x,k)
            test.nwin <- c(test.nwin, nrow(m.count.x))
            test.rpw <- c(test.rpw, mean(m.count.x[,2], na.rm=T))
            m.count.y <- hist.ch(y,k)
            control.nwin <- c(control.nwin, nrow(m.count.y))
            control.rpw <- c(control.rpw, mean(m.count.y[,2], na.rm=T))
            m <- join.mat(m.count.x, m.count.y)
            opt.results <- optimal(m)
            test.cv <- c(test.cv, opt.results$ML.Test)
            test.aic <- c(test.aic, opt.results$AIC.Test)
            control.cv <- c(control.cv, opt.results$ML.Normal)
            control.aic <- c(control.aic, opt.results$AIC.Normal)
        } # end for k win.size

        result <- list(test.cv=test.cv, test.aic=test.aic, control.cv=control.cv, control.aic=control.aic)
        win.stat <- list(test.win=test.nwin, test.rpw=test.rpw, control.nwin=control.nwin, control.rpw=control.rpw)
        attr(result,"winsize") <- win.size
        attr(result,"winstat") <- win.stat
        attr(result,"class") <- c("optwin")
        result
}

# Calculating the AIC and cross validation log-likelihood
# from input x (read counts) in a single sample

optimal.onesample = function(x){


        read.count = x[,2]


        cross = function(x){
            # this function applies leave-one-out cross-validation and calculate the likelihood value
            # x is reads count
            # L is window size

            np = nrow(x)     # Windows number
            h = 1/np         # Fraction
            T = x[x>1]       # Removing empty and single read windows
            N = sum(T)       # Total reads number without empty and single read windows
            n = nrow(T)    # non-single read windows number
            f = T * log(T - 1) - T * log((N-1)*h)
            f2 = T * log(T) - T * log((N)*h)
            ml = sum(f)       # The log likelihood value
            ml2 = sum(f2)
            return(list(ml = ml, ml2=ml2, np = np))
        } # end cross functions.


        ml = c()
        ml2 = c()
        np = c()
        res = cross(data.frame(read.count))
        ml = res$ml
        ml2 = res$ml2
        np = res$np


        aic. = - 2 * ml2 + 2 * np


        return(list(ML.Test = ml, AIC.Test = aic.))

}

# function to read the starting positions
# Note, if argument select is given (non NULL),
# the argumet excludeXY is ignored.
# Also, if select is null and excludeXY=F,
# we effectively take all reads starting positions

read.pos = function(x, excludeXY=T, header=F, select=NULL, ...){

        file.pos <- read.table(x,header=header,...)

        if(is.null(select)){
            if(excludeXY){
                select = c("chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
            } # end if excludeXY
        } # end if is.null select

        file.pos <- file.pos[file.pos[,1]%in%select,]
        return(file.pos)
    }
