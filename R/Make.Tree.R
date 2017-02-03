#' \code{Make.Tree.R} Generate a phylogenetic tree to and plot
#'
#' @description Takes a series of sequences and uses muscle and make_phylogeny (QIIME). Can be used with MacQIIME or QIIME. Will default to MacQIIME locations. A file called Feature_Seqs.tree will be created in your working directory for future use.  Uses midpoint rooting, and fasttree. NOTE: This is using system commands... might only work on mac and linux? If you are using windows reconsider your life choices.
#' @param FEATURENAMES The names of your features/OTUs/SVs. Will be the row names of your input table in most cases.
#' @param FEATURESEQS The corresponding sequences for each feature. If this is a dada2 table, then feature names are also feature seqs.
#' @param MUSCLE (optional) The location of muscle, defaults to /macqiime/bin/muscle
#' @param MAKEPHY (optional) Location of qiime make_phylogeny.py, defatuls to /macqiime/anaconda/bin/make_phylogeny.py
#' @param PRINTTREE (optional) Would you like a tree to be printed TRUE/FALSE? defaults to TRUE
#'
#' @return Returns tree object and prints a phylogenetic tree
#' @export

Make.Tree<-function(FEATURENAMES,FEATURESEQS, MUSCLE, MAKEPHY){

  if(missing(MUSCLE)){MUSCLE="/macqiime/bin/muscle"}
  if(missing(MAKEPHY)){MAKEPHY="/macqiime/anaconda/bin/make_phylogeny.py"}
  if(missing(PRINTTREE)){PRINTTREE="yes"}

  write.fasta(as.list(FEATURENAMES), FEATURENAMES, "Feature_Seqs.fasta")

  system(paste(MUSCLE,"-in Feature_Seqs.fasta -out Feature_Seqs.mfa"))
  system(paste("source /macqiime/configs/bash_profile.txt;",MAKEPHY,"--tree_method fasttree --root_method=midpoint -i Feature_Seqs.mfa -o Feature_Seqs.tree"))


  TREE<-ape::read.tree("Feature_Seqs.tree")

if (PRINTTREE==TRUE){
  plot.phylo(TREE, cex=0.3, direction="downwards")
}

  return(TREE)

}
