#' Make.Percent
#'
#' Takes a table of Features/OTUs/SVs and converts to Percentages (Sum to 100).
#'
#' @param FEATURES Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @usage PercentTab<-Make.Percent(FeatureTable)
#' @return Table of Percents
#' @export

Make.Percent<-function(FEATURES){
  FEATURES<-TidyConvert.ToMatrix(FEATURES, colnames(FEATURES)[1])
  Percent<-apply(FEATURES,2, function(x){100*(x/sum(x))})
  return(Percent)
}
