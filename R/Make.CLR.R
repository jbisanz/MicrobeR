#' Make.CLR
#'
#' Takes a table of Features/OTUs/SVs and converts to centered log2 ratio abundances. See:
#'  Compositional analysis: a valid approach to analyze microbiome high-throughput sequencing data
#'   Gregory B. Gloor, Gregor Reid; Canadian Journal of Microbiology, 2016, 62:692-703, 10.1139/cjm-2015-0821
#' For information on CZM replacement see : Martin-Fernandez, J.A., Hron, K., Templ, M., Filzmoser, P., Palarea-Albaladejo, J. Bayesian-multiplicative treatment of count zeros in compositional data sets. Statistical Modelling 2015; 15 (2): 134-158.
#' @param FEATURES Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param PRIOR A numeric value to add before log transformation (default 0.5). *Ignored if CZM=TRUE
#' @param CZM Should count zero multiplicative method be used instead of adding simple prior (TRUE/FALSE, defaults to FALSE)
#' @return Table of CLR-normalized Abundances (log2)
#' @export

Make.CLR<-function(FEATURES, PRIOR, CZM){
  FEATURES<-TidyConvert.ToMatrix(FEATURES, colnames(FEATURES)[1])

  if(nrow(FEATURES)<100){message("WARNING: CLR being applied with relatively few features.")}

  if(missing(CZM)){CZM=FALSE}
  if(missing(PRIOR)){PRIOR=0.5}

  if(CZM==TRUE){
    FEATURES <- t(cmultRepl(t(FEATURES),  label=0, method="CZM", output="counts", suppress.print=T))
  } else {
    FEATURES<-FEATURES + PRIOR
  }

  FEATURES<-apply(FEATURES,2, function(column){log2(column)-mean(log2(column))})
  return(FEATURES)
}
