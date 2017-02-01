#' \code{Make.CLR} Make a CLR-normalized Feature/OTU/SV table
#' @description  Takes a table of Features/OTUs/SVs and converts centered log ratio abundances. See:
#'  Compositional analysis: a valid approach to analyze microbiome high-throughput sequencing data
#'   Gregory B. Gloor, Gregor Reid; Canadian Journal of Microbiology, 2016, 62:692-703, 10.1139/cjm-2015-0821
#'
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @return Table of CLR-normalized Abundances
#' @export

Make.CLR<-function(TABLE){
  TABLE<-TABLE+0.5
  TABLE<-apply(log2(TABLE),2, function(column){(column-mean(column))})
  return(TABLE)
}
