#' Subsample.Table
#'
#' @description Takes a Feature/OTU/SV table of counts where samples are columns and feature names are row names and returns a table with even coverage on a per-sample basis. Now this function is really just an alias for rarefy_even_depth from phyloseq or rtk from rtk. NOTE: Sampling with replacement for sinlge rarefraction method!
#'
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param DEPTH Count depth, defaults to min(colSums(OTUTABLE)) if not passed
#' @param SEED A randomization SEED, defaults to 182. This is ignored for multiple subsamples and instead the seeds 1:NSAMPS is used.
#' @param NSAMPS The number of samples that should be taken of the table for multiple subsamples. Use an odd number to avoid decimals. Default=101
#' @param THREADS Number of CPUs to use for multiple rarefraction, defaults to 2/3 of available.
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

Multiple.Subsample.Table<-function(FEATURES,DEPTH, VERBOSE, NSAMPS, THREADS, AGGFUNCTION){
  if(missing(VERBOSE)){VERBOSE=T}
  if(missing(DEPTH)){DEPTH=min(colSums(FEATURES))}
  if(missing(NSAMPS)){NSAMPS=101}
  
  if(NSAMPS%%2==0){stop("NSAMPS must be an odd number to prevent fractions in read counts.")}
  
  if(missing(THREADS)){THREADS=round(parallel::detectCores()*2/3, 0)}
  if(missing(AGGFUNCTION)){AGGFUNCTION="median"}
  
  if(!AGGFUNCTION %in% c("median","mean")){stop("Aggregation function (AGGFUNCTION) must be either mean or median")}
  
  if(VERBOSE){message(paste("Subsampling feature table to", DEPTH, ", currently has ", nrow(FEATURES), " taxa."))}
  
  tables<-rtk::rtk(FEATURES, repeats=NSAMPS, depth=DEPTH, threads=THREADS, ReturnMatrix = NSAMPS, verbose=FALSE, margin=2)$raremat
  
  if(VERBOSE){message("Merging individual features tables:")}
  subsampled.FEATURES<-matrix(nrow=nrow(FEATURES), ncol=ncol(FEATURES))
  rownames(subsampled.FEATURES)<-rownames(FEATURES)
  colnames(subsampled.FEATURES)<-colnames(FEATURES)
  for(i in 1:ncol(FEATURES)){
    if(VERBOSE){message(i/ncol(FEATURES)*100, "%") }
    subsampled.FEATURES[,i]<-  
      lapply(tables, function(x) x[,i]) %>%
      do.call(cbind, .) %>%
      apply(., 1, AGGFUNCTION)
  }
  
  subsampled.FEATURES<-subsampled.FEATURES[rowSums(subsampled.FEATURES)>0,]
  
  if(VERBOSE){message(paste("...multiply sampled to",DEPTH, "with the median feature count reported. A total of",nrow(subsampled.FEATURES), "taxa have been returned."))}
  if(VERBOSE){print(summary(colSums(subsampled.FEATURES)))}
  
  return(subsampled.FEATURES)
}
