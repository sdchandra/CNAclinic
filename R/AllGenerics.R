# Access the genomic bin information
setGeneric("chromosomes", function(object) standardGeneric("chromosomes"))

setGeneric("bpstart", function(object) standardGeneric("bpstart"))
setGeneric("bpend", function(object) standardGeneric("bpend"))
setGeneric("usebin", function(object) standardGeneric("usebin"))
setGeneric("bins", function(object) standardGeneric("bins"))

# Access sample names
setGeneric("sampleNames", function(object) standardGeneric("sampleNames"))

# Access log2R and segment values as well as CNA calls
setGeneric("copyNumber", function(object) standardGeneric("copyNumber"))
setGeneric("segHMM", function(object) standardGeneric("segHMM"))
setGeneric("segCBS", function(object) standardGeneric("segCBS"))
setGeneric("segLACBS", function(object) standardGeneric("segLACBS"))
setGeneric("segPLS", function(object) standardGeneric("segPLS"))
setGeneric("segSummary", function(object) standardGeneric("segSummary"))
setGeneric("calls", function(object) standardGeneric("calls"))

# Get the names of empty & non-empty slots in CNAclinicData objects
setGeneric("emptySlots", function(object) standardGeneric("emptySlots"))
setGeneric("filledSlots", function(object) standardGeneric("filledSlots"))

# Get summary and overview of CNAclinic object
setGeneric("summary", function(object, ...)
    standardGeneric("summary"))
# Get summary and overview of CNAclinic object
setGeneric("show", function(object)
    standardGeneric("show"))


# Plotting data
setGeneric("plotSampleData", function(object, ...)
    standardGeneric("plotSampleData"))
setGeneric("plotMultiSampleData", function(object, ...)
    standardGeneric("plotMultiSampleData"))

# Calculate consensus segments
setGeneric("summSegments", function(object, ...)
    standardGeneric("summSegments"))

# Call the copy number segments
setGeneric("callData", function(object, ...)
    standardGeneric("callData"))

# Manipulate CNAclinicData objects
setGeneric("subsetData", function(object, ...)
    standardGeneric("subsetData"))
setGeneric("combineData", function(x, y)
    standardGeneric("combineData"))

# Get gene information
setGeneric("getGeneInfo", function(geneID, ...)
    standardGeneric("getGeneInfo"))

# Get QC and quantification statistics
setGeneric("statsCNA", function(object, ...)
    standardGeneric("statsCNA"))

# Export data
setGeneric("exportData", function(x, ...)
    standardGeneric("exportData"))




setGeneric("usebin<-", function(object, value)
    standardGeneric("usebin<-"))
setGeneric("sampleNames<-", function(object, value)
    standardGeneric("sampleNames<-"))
setGeneric("copyNumber<-", function(object, value)
    standardGeneric("copyNumber<-"))
setGeneric("segHMM<-", function(object, value)
    standardGeneric("segHMM<-"))
setGeneric("segCBS<-", function(object, value)
    standardGeneric("segCBS<-"))
setGeneric("segLACBS<-", function(object, value)
    standardGeneric("segLACBS<-"))
setGeneric("segPLS<-", function(object, value)
    standardGeneric("segPLS<-"))
setGeneric("segSummary<-", function(object, value)
    standardGeneric("segSummary<-"))
setGeneric("calls<-", function(object, value)
    standardGeneric("calls<-"))

# EOF
