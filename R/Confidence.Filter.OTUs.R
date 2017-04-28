#' Filter Features (OTUS/SVs) such that they are present in at least X samples with at least a total of Y reads
#'
#' \code{Confidence.Filter.OTUs} Takes a table of FeaturesOTUs/SVs and filters for such that they are present in at least X samples with at least a total of Y reads
#'
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and feature IDs are row names
#' @param MINSAMPS A minimum number of samples for the feature/OTU/SV to be observed in
#' @param MINREADS A minimum number of reads for a feature/OTUs/SV to be kept
#' @return filtered.OTUTABLE
#' @export

Confidence.Filter.OTUs<-function(OTUTABLE,MINSAMPS,MINREADS){
  total_observations<-sum(OTUTABLE)  #get total number of reads

  print(paste("Filtering OTUs such that they are present in at least", MINSAMPS, " samples with a total of at least ", MINREADS, " reads"))
  print(paste("...There are", total_observations, "reads and", nrow(OTUTABLE), " OTUs"))

  filtered.OTUTABLE<-OTUTABLE[apply(OTUTABLE,1,function(x){length(grep("TRUE",x!=0))}>=MINSAMPS),]
  filtered.OTUTABLE<-filtered.OTUTABLE[((rowSums(filtered.OTUTABLE))>=MINREADS),]

  print(paste("...After filtering there are", sum(filtered.OTUTABLE), "reads and", nrow(filtered.OTUTABLE), "OTUs"))
  return(filtered.OTUTABLE)
}
