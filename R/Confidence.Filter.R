#' Filter features such that they are present in at least X samples with at least a total of Y reads
#'
#' Takes a table of features and filters for such that they are present in at least X samples with at least a total of Y reads
#'
#' @param FEATURES Table of feature/OTU/SV counts where Samples are columns, and feature IDs are row names
#' @param MINSAMPS A minimum number of samples for the feature/OTU/SV to be observed in
#' @param MINREADS A minimum number of reads for a feature/OTUs/SV to be kept across all samples
#' @param VERBOSE Should summary be printed? (T/F)
#' @return filtered feature table
#' @export

Confidence.Filter<-function(FEATURES,MINSAMPS,MINREADS,VERBOSE){
  FEATURES<-TidyConvert.ToMatrix(FEATURES, colnames(FEATURES)[1])
  total_observations<-sum(FEATURES)  #get total number of reads

  if(missing(VERBOSE)){VERBOSE=T}
  if(VERBOSE==T){
    message(paste("Filtering features such that they are present in at least", MINSAMPS, "samples with a total of at least", MINREADS, "reads."))
    message(paste("...There are", total_observations, "reads and", nrow(FEATURES), " features"))
  }

  filtered<-FEATURES[apply(FEATURES,1,function(x){length(grep("TRUE",x!=0))}>=MINSAMPS),]
  filtered<-filtered[((rowSums(filtered))>=MINREADS),]

  if(VERBOSE==T){
    message(paste("...After filtering there are", sum(filtered), "reads and", nrow(filtered), "features"))
  }

  return(filtered)
}
