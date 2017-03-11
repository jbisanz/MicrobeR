#' \code{Merge.Replicates} Sum counts together for replicates (ex. pcr/extraction/sequencing run)
#'
#' @description Collapses a table via aggregate FUN=sum based on a provided column in the datatable.
#'
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param METADATA A fraction to use to filter
#' @param CATEGORY Minimum number of reads to pass filter as a string, ex. "RepID"
#' @return A list of [[1]] new Feature/OTU/SV table and [[2]] collapsed metadata file
#' @export

Merge.Replicates<-function(OTUTABLE, METADATA, CATEGORY){
print(paste("Found",length(METADATA[,CATEGORY])-length(unique(METADATA[,CATEGORY])), "samples with replicates. Suming reads."))
OTUTABLE<-as.data.frame(t(OTUTABLE))
OTUTABLE$TEMPID<-METADATA[rownames(OTUTABLE),CATEGORY]
rownames(OTUTABLE)<-NULL
OTUTABLE<-aggregate(. ~ TEMPID, data=OTUTABLE, FUN=sum)
rownames(OTUTABLE)<-OTUTABLE$TEMPID
OTUTABLE$TEMPID<-NULL
rownames(METADATA)<-NULL
METADATA<-METADATA[!duplicated(METADATA[,CATEGORY]),]
rownames(METADATA)<-METADATA[,CATEGORY]
OTUTABLE<-t(OTUTABLE)
return(list(OTUTABLE,METADATA))
}
