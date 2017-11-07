
#' Get transcript and position information for NCBI Entrez gene IDs
#'
#' A helper function to generate gene information as input for plotSampleData
#' function.
#'
#' @param geneID a \code{\link[base]{numeric}} vector containing NCBI Entrez
#' gene IDs
#' @param genome Genome build used to align sequencing reads.
#' @param select When select is \code{"first"} (default), the first
#' transcript of the input genes are chosen. If \code{"all"}, information of
#' all possible transcripts will be contained in output.
#'
#' @return a data.frame containing the information for the \code{geneID} input.
#' The data.frame has the columns named \code{chromosome}, \code{start}, \code{end}
#' @export
#' @import annotate org.Hs.eg.db
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#'
#' @examples
#'      \dontrun{
#'       vignette("CNAclinic")
#'      }
#'

getGeneInfo = function(geneID, genome,
    transcript="first"){

    genome <- match.arg(genome, c("hg19", "hg38"))
    transcript <- match.arg(transcript, c("first", "all"))

    # Map GENE ID to GENE SYMBOL and Position (one to one)

    a <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
        keys = as.character(geneID),
        columns=c("ENTREZID", "SYMBOL"),
        keytype="ENTREZID")

    # Map GENE ID to transcript information (one to many)

    if(genome == "hg19"){
        b <- AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg19.knownGene,
            keys = as.character(geneID),
            columns=c("GENEID", "TXCHROM", "TXSTRAND", "TXSTART", "TXEND", "TXID"),
            keytype="GENEID")
    }else if(genome == "hg38"){
        b <- AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg38.knownGene,
            keys = as.character(geneID),
            columns=c("GENEID", "TXCHROM", "TXSTRAND", "TXSTART", "TXEND", "TXID"),
            keytype="GENEID")
    }


    names(b) <- c("ENTREZID", "TXID", "TXCHROM", "TXSTRAND", "TXSTART", "TXEND")
    df <- merge(a, b, "ENTREZID")

    if(transcript == "first"){
        df <- df[!duplicated(df$ENTREZID), ]
    }else if(transcript == "all"){
        df <- df[!duplicated(df), ]
    }

    df <- df[ , c("ENTREZID", "SYMBOL", "TXID", "TXCHROM", "TXSTRAND", "TXSTART", "TXEND")]
    names(df) <- c("geneID", "geneSymbol", "txID", "chromosome", "strand", "start", "end")
    return(df)
}


