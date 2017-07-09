#' \code{Make.Percent} Takes a table of Features/OTUs/SVs and converts to Percentages (Sum to 100).
#'
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @usage PercentTab<-Make.Percent(FeatureTable)
#' @return Table of Percents
#' @export

Make.Percent<-function(TABLE){
  Percent<-apply(TABLE,2, function(x){100*(x/sum(x))})
  return(Percent)
}
