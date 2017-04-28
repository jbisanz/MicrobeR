#' \code{Subsample.Table} Rarify/Subsample a Feature
#' @description Takes a Feature/OTU/SV table where samples are columns and feature names are row names. This was formerly using the VEGAN rrarefy function but was found to be slow on large numbers of samples (~3000). It is coded using a loop while the implementation in phyloseq is vectorized it would appear. Now this function is really just an alias for rarefy_even_depth from phyloseq. NOTE: Sampling with replacement!
#'
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param DEPTH Count depth, defaults to min(colSums(OTUTABLE)) if not passed
#' @param SEED A randomization SEED, defaults to 182 as tribute to scarface and blink.
#' @return Subsampled Table
#' @export

Subsample.Table<-function(OTUTABLE,DEPTH, SEED){

  if(missing(DEPTH)){DEPTH=min(colSums(OTUTABLE))}
  if(missing(SEED)){SEED=182}

  print(paste("Subsampling OTU table to", DEPTH, ", currently has ", nrow(OTUTABLE), " taxa."))
  subsampled.OTUTABLE<-phyloseq::rarefy_even_depth(otu_table(OTUTABLE, taxa_are_rows = T), sample.size=DEPTH, rngseed=SEED, verbose = FALSE) #expecting the transpose for otu table layout so transpose then transpose back
  print(paste("...sampled to",DEPTH, "reads with", nrow(subsampled.OTUTABLE), "taxa"))


  return(subsampled.OTUTABLE)
}
