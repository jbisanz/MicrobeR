#' \code{Filter.OTUs.Fraction} Remove Features/OTUs/SVs representing less than a certain specified fraction.
#'
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names.
#' @param MINFRACTION A fraction to use to filter.
#' @return Filtered Table
#' @export


Filter.OTUs.Fraction<-function(OTUTABLE,MINFRACTION){ #filters OTUs in an OTU table such that they are removed if they are less than MINFRACTION of reads, similar to filter_otus_from_otutable.py in QIIME
  total_observations<-sum(OTUTABLE)  #get total number of reads

  print(paste("Filtering OTU table at a min fraction of", MINFRACTION, "of OTU table..."))
  print(paste("...There are", total_observations, "reads and", nrow(OTUTABLE), " OTUs"))
  otu_sums<-rowSums(OTUTABLE)

  filtered.OTUTABLE<-OTUTABLE[((otu_sums/total_observations)>=MINFRACTION),]

  print(paste("...After filtering there are", sum(filtered.OTUTABLE), "reads and", nrow(filtered.OTUTABLE), "OTUs"))
  return(filtered.OTUTABLE)
} #take an otutable and minimum fraction and return a filtered otutable
