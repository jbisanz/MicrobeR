#' Make.Tree
#'
#' Generate a phylogenetic tree and optional plot.
#'
#' @description Takes a series of sequences and uses DECIPHER and FastTree (Make.Tree.Fast) OR DECIPHER and phangorn (Make.Tree.R) as per https://f1000research.com/articles/5-1492/v2.
#' @param FEATURENAMES The names of your features/OTUs/SVs. Will be the row names of your input table in most cases.
#' @param FEATURESEQS The corresponding sequences for each feature. If this is a dada2 table, then feature names are also feature seqs.
#' @param FASTTREE (optional) The location of the fast tree binary on your system
#'
#' @return Returns tree object and prints a phylogenetic tree
#' @export

Make.Tree.R<-function(FEATURENAMES,FEATURESEQS){

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

Make.Tree.Fast<-function(FEATURENAMES,FEATURESEQS, FASTTREE){
  
  if(missing(FASTTREE)){FASTTREE="FastTree"}
  
  fasta<-DNAStringSet(FEATURESEQS)
  names(fasta)<-FEATURENAMES
  
  alignment <- AlignSeqs(fasta, anchor=NA)
  randomname<-paste(sample(LETTERS, 20), collapse="")
  writeXStringSet(alignment, paste0("/tmp/",randomname,".mfa"))
  system(paste0(FASTTREE, " -nt ", paste0("/tmp/",randomname,".mfa"), " > ", paste0("/tmp/",randomname,".tree")))
  TREE<-read.tree(paste0("/tmp/",randomname,".tree"))
  TREE<-midpoint(TREE) #### ADD MIDPOINT TO EXPORTS
  return(TREE)
}

