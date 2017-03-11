#' \code{Subsample.Table} Rarify/Subsample a Feature
#' @description Takes a Feature/OTU/SV table where samples are columns and feature names are row names, subsamples using VEGAN implementation using randomization seed if provided. Note: samples without replacement
#'
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param DEPTH Count depth, defaults to min(colSums(OTUTABLE)) if not passed
#' @param SEED A randomization SEED, defaults to not specificed for VEGAN.
#' @return Subsampled Table
#' @export

Subsample.Table<-function(OTUTABLE,DEPTH, SEED){

  #test=data.frame(Samp1=sample(1:100,20), Samp2=sample(1:100,20))

  if(missing(DEPTH)){DEPTH=min(colSums(OTUTABLE))}
  if(!missing(SEED)){set.seed(SEED)}

  print(paste("Subsampling OTU table to", DEPTH))
  subsampled.OTUTABLE<-t(rrarefy(t(OTUTABLE),DEPTH)) #expecting the transpose for otu table layout so transpose then transpose back
  subsampled.OTUTABLE<-subsampled.OTUTABLE[rowSums(subsampled.OTUTABLE)>0,] #remove now zero count OTUs
  print(paste("...sampled to",DEPTH, "reads with", nrow(subsampled.OTUTABLE), "taxa"))
  return(subsampled.OTUTABLE)
}
