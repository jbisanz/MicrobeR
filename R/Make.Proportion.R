#' Make.Proportion
#'
#' Takes a table of Features/OTUs/SVs and converts to Proportion (Sum to 1).
#'
#' @param FEATURES Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @usage ProportionTab<-Make.Proportion(FeatureTable)
#' @return Table of Proportions
#' @export

Make.Proportion<-function(FEATURES){
  FEATURES<-TidyConvert.ToMatrix(FEATURES, colnames(FEATURES)[1])
  Proportion<-apply(FEATURES,2, function(x){(x/sum(x))})
  return(Proportion)
}
