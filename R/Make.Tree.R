#' Make.Tree
#'
#' Generate a phylogenetic tree and optional plot.
#'
#' @description Takes a series of sequences and uses muscle and make_phylogeny (Make.Tree.QIIME) OR DECIPHER and phangorn (Make.Tree.R) as per https://f1000research.com/articles/5-1492/v2. The R-method will take considerably longer and becomes unusable for >a few hundred features. For QIIME method: can be used with MacQIIME or QIIME. Will default to MacQIIME locations. A file called Feature_Seqs.tree will be created in your working directory for future use.  Uses midpoint rooting, and fasttree. NOTE: This is using system commands... probabily will only work on mac and linux. If you are using windows reconsider your life choices and consider dual booting Ubuntu.
#' @param FEATURENAMES The names of your features/OTUs/SVs. Will be the row names of your input table in most cases.
#' @param FEATURESEQS The corresponding sequences for each feature. If this is a dada2 table, then feature names are also feature seqs.
#' @param MUSCLE (optional) The location of muscle, defaults to /macqiime/bin/muscle
#' @param MAKEPHY (optional) Location of qiime make_phylogeny.py, defatuls to /macqiime/anaconda/bin/make_phylogeny.py
#' @param PRINTTREE (optional) Would you like a tree to be printed TRUE/FALSE? defaults to TRUE
#'
#' @return Returns tree object and prints a phylogenetic tree
#' @export

Make.Tree.QIIME<-function(FEATURENAMES,FEATURESEQS, MUSCLE, MAKEPHY, PRINTTREE){

  if(missing(MUSCLE)){MUSCLE="/macqiime/bin/muscle"}
  if(missing(MAKEPHY)){MAKEPHY="/macqiime/anaconda/bin/make_phylogeny.py"}
  if(missing(PRINTTREE)){PRINTTREE=F}

  fasta<-DNAStringSet(FEATURESEQS)
  names(fasta)<-FEATURENAMES
  writeXStringSet(fasta, "Feature_Seqs.fasta")

  system(paste(MUSCLE,"-in Feature_Seqs.fasta -out Feature_Seqs.mfa"))
  system(paste("source /macqiime/configs/bash_profile.txt;",MAKEPHY,"--tree_method fasttree --root_method=midpoint -i Feature_Seqs.mfa -o Feature_Seqs.tree"))

  TREE<-ape::read.tree("Feature_Seqs.tree")

if (PRINTTREE==TRUE){
  plot.phylo(TREE, cex=0.3, direction="downwards")
}

  return(TREE)

}


Make.Tree.R<-function(FEATURENAMES,FEATURESEQS, MUSCLE, MAKEPHY, PRINTTREE){

  fasta<-DNAStringSet(FEATURESEQS)
  names(fasta)<-FEATURENAMES

  alignment <- AlignSeqs(fasta, anchor=NA)

  phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phang.align)
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  fit = pml(treeNJ, data=phang.align)
  ## negative edges length changed to 0!
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
  TREE<-fitGTR$tree

  if (PRINTTREE==TRUE){
    plot.phylo(TREE, cex=0.3, direction="downwards")
  }

  return(TREE)

}


