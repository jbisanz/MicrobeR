#' Subsample.Table
#'
#' @description Takes a Feature/OTU/SV table where samples are columns and feature names are row names. This was formerly using the VEGAN rrarefy function but was found to be slow on large numbers of samples (~3000). Suspected this was due to serial/vectorized implimentations. Now this function is really just an alias for rarefy_even_depth from phyloseq. NOTE: Sampling with replacement!
#'
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param DEPTH Count depth, defaults to min(colSums(OTUTABLE)) if not passed
#' @param SEED A randomization SEED, defaults to 182.
#' @param VERBOSE Should progress and metrics be printed to screen via message()? Default=TRUE
#' @return Subsampled Table
#' @export

Subsample.Table<-function(FEATURES,DEPTH, SEED, VERBOSE){
  if(missing(VERBOSE)){VERBOSE=T}
  if(missing(DEPTH)){DEPTH=min(colSums(FEATURES))}
  if(missing(SEED)){SEED=182}

  if(VERBOSE==T){message(paste("Subsampling feature table to", DEPTH, ", currently has ", nrow(FEATURES), " taxa."))}
  subsampled.FEATURES<-as.data.frame(phyloseq::rarefy_even_depth(otu_table(FEATURES, taxa_are_rows = T), sample.size=DEPTH, rngseed=SEED, verbose = FALSE)) #expecting the transpose for otu table layout so transpose then transpose back
  if(VERBOSE==T){message(paste("...sampled to",DEPTH, "reads with", nrow(subsampled.FEATURES), "taxa"))}


  return(subsampled.FEATURES)
}
