#' \code{Read.Filter} Remove samples with less than a user specified number of reads.
#'
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param METADATA A fraction to use to filter
#' @param READCUTOFF Minimum number of reads to pass filter
#' @return A list containing new OTUtable and Metadata. Also creates a histogram of read counts
#' @export

Read.Filter<-function(OTUTABLE,METADATA,READCUTOFF){
  message("Removed the following ", sum(colSums(OTUTABLE)<READCUTOFF), " samples with less than ", READCUTOFF, " reads:")
  print(colnames(OTUTABLE[,colSums(OTUTABLE)<READCUTOFF]))
  OTUTABLE<-OTUTABLE[,colSums(OTUTABLE)>=READCUTOFF]
  message("This left the following samples, check to ensure that no control samples remain:")
  print(colnames(OTUTABLE)[order(colnames(OTUTABLE))])
  METADATA<- METADATA[colnames(OTUTABLE),]
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
        + theme_bw()
      )

  return(list(OTUTABLE, METADATA))
}
