#' Filter.OTUs.Fraction
#'
#' Remove Features/OTUs/SVs representing less than a certain specified fraction of the total dataset. Similar to filter_otus_from_otutable.py in QIIME.
#'
#' @param FEATURES Table of feature/OTU/SV counts where Samples are columns, and IDs are row names.
#' @param MINFRACTION A fraction to use to filter.
#' @param VERBOSE Should summary be printed? (T/F)
#' @return Filtered Table
#' @export

Fraction.Filter<-function(FEATURES,MINFRACTION, VERBOSE){ #filters OTUs in an OTU table such that they are removed if they are less than MINFRACTION of reads, similar to filter_otus_from_otutable.py in QIIME
  FEATURES<-TidyConvert.ToMatrix(FEATURES, colnames(FEATURES)[1])
  total_observations<-sum(FEATURES)  #get total number of reads

  if(missing(VERBOSE)){VERBOSE=T}
  if(VERBOSE==T){
    print(paste("Filtering table at a min fraction of", MINFRACTION, "of feature table..."))
    print(paste("...There are", total_observations, "reads and", nrow(FEATURES), " features"))
  }

  feature_sums<-rowSums(FEATURES)
  filtered<-FEATURES[((feature_sums/total_observations)>=MINFRACTION),]

  if(VERBOSE==T){
    print(paste("...After filtering there are", sum(filtered), "reads and", nrow(filtered), "OTUs"))
  }
    return(filtered)
}
