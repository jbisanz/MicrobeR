#' \code{Read.Filter} Remove samples with less than a user specified number of reads.
#'
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names.
#' @param METADATA Metadata file where samples are row names.
#' @param READCUTOFF Minimum number of reads to pass filter
#' @param VERBOSE Should lists of which samples were removed or kept be printed to screen? (defaults to TRUE)
#' @param PLOT Should a plot of read depth be made? (defaults to TRUE)
#' @usage filter<-Read.Filter(svtable, metadata, 5000)
#' newtable<-filter$Features
#' newmetadata<-filter$Metadata
#' @return A named list containing new feature table and metadata.
#' @export

Read.Filter<-function(OTUTABLE,METADATA,READCUTOFF, VERBOSE, PLOT){
  if(missing(VERBOSE)){VERBOSE=T}
  if(missing(PLOT)){PLOT=T}

  if(sum(!colnames(OTUTABLE) %in% rownames(METADATA))>0){stop("There are samples in feature table with no corresponding metadata!")}

  if(VERBOSE==T){
    message("Removed the following ", sum(colSums(OTUTABLE)<READCUTOFF), " samples with less than ", READCUTOFF, " reads:")
    print(colnames(OTUTABLE[,colSums(OTUTABLE)<READCUTOFF]))
  }
  OTUTABLE<-OTUTABLE[,colSums(OTUTABLE)>=READCUTOFF]
  if(VERBOSE==T){
    message("This left the following samples, check to ensure that no control samples remain:")
    print(colnames(OTUTABLE)[order(colnames(OTUTABLE))])
  }
  METADATA<- METADATA[colnames(OTUTABLE),]

  if(PLOT==T){
  readcounts<-as.data.frame(colSums(OTUTABLE))
  colnames(readcounts)[1]<-"reads"
  print(ggplot(data=readcounts, aes(x=reads))
        + geom_freqpoly(binwidth=2000)
        + ggtitle("Read Depth Post-filter")
        + theme(legend.title=element_blank(),
                axis.text.x=element_text(size=6, angle=45, hjust=1),
                axis.text.y=element_text(size=6),
                legend.text=element_text(size=6),
                axis.title.x=element_text(size=7),
                axis.title.y=element_text(size=7),
                plot.title=element_text(size=8))
        + xlab("Number of Reads")
        + ylab("Number of Samples")
        + theme_bw()
      )}

  return(list(Features=OTUTABLE, Metadata=METADATA))
}
