#' \code{PCoA} Plot PCoAs of Bray-Curtis and (Un)Weighted UniFrac using first 3 axes.
#'
#' @description Will create a 2x6 PDF for each condition passed. UniFrac is implimented as per GUniFrac and PCoA from APE. Plotting via ggplot2.
#'
#' @param PDFTITLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param METADATA Metadata Table with variables to color PCoAs by
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param TREE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param CONDS A vector of column names from METADATA to use: ex c("Treatment","Time")
#' @param SUBSAMPLE Should table be subsampled? TRUE/FALSE (default: TRUE)
#'
#' @usage PCoA("Treatment_and_Time.pdf", experiment_metadata, otutable, ggtree, c("Treatment","Time","ExtractionPlate"), FALSE)
#' @export

PCoA<-function(PDFTITLE,METADATA,OTUTABLE, TREE, CONDS, SUBSAMPLE){

  if(missing(SUBSAMPLE) || SUBSAMPLE==TRUE){
      DEPTH<-min(colSums(OTUTABLE))
      sub.OTUlevel<-Subsample.Table(OTUTABLE,DEPTH)
      filtered.TREE<-prune_taxa(rownames(sub.OTUlevel), tree)
  } else {
     sub.OTUlevel<-OTUTABLE
     filtered.TREE<-TREE
  }

  print("Doing UniFrac...")
  unifracs <- GUniFrac::GUniFrac(t(sub.OTUlevel),filtered.TREE, alpha=c(0, 0.5, 1))$unifracs
  weighted.UniFrac <- unifracs[, , "d_1"] # Weighted UniFrac as per ?GUniFrac
  unweighted.UniFrac <- unifracs[, , "d_UW"] # Unweighted UniFrac
  print("Doing Bray Curtis...")
  braycurtis<-vegan::vegdist(t(sub.OTUlevel), method="bray", binary=FALSE, diag=FALSE, upper=FALSE)
  braycurtis<-as.matrix(braycurtis)

  print("Doing PCoA analysis...")
  pco.braycurtis<-ape::pcoa(braycurtis) #using the principal coordinants analysis of the ape package
  pco.weighted.UniFrac<-ape::pcoa(weighted.UniFrac)
  pco.unweighted.UniFrac<-ape::pcoa(unweighted.UniFrac)

  print ("Plotting...")
  plot.pco.braycurtis<-merge(as.data.frame(pco.braycurtis$vectors), METADATA[rownames(pco.braycurtis$vectors),], by="row.names", all.x=T)
  plot.pco.weighted.UniFrac<-merge(as.data.frame(pco.weighted.UniFrac$vectors), METADATA[rownames(pco.weighted.UniFrac$vectors),], by="row.names", all.x=T)
  plot.pco.unweighted.UniFrac<-merge(as.data.frame(pco.unweighted.UniFrac$vectors), METADATA[rownames(pco.unweighted.UniFrac$vectors),], by="row.names", all.x=T)

  print(paste("Making PDF in current working directory called:", PDFTITLE))
  pdf(PDFTITLE, height=10, width=8)
  for (i in 1:length(CONDS)){
    curcond<-CONDS[i]
    print(paste("...Plotting: ", curcond))
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(3, 2)))
    print(ggplot(plot.pco.braycurtis, aes(x=Axis.1, y=Axis.2, color=plot.pco.braycurtis[,curcond])) + geom_point(alpha=0.5) + theme_bw() + theme(plot.title=element_text(size=10), axis.text=element_text(size=8), aspect.ratio = 1, axis.title=element_text(size=8,face="bold")) + ggtitle("Bray Curtis, PCo1 vs PCo2") + xlab(paste("PCo1:", round((100*pco.braycurtis$values$Eigenvalues/sum(pco.braycurtis$values$Eigenvalues))[1],1), "% Variation Explained")) + ylab(paste("PCo2:", round((100*pco.braycurtis$values$Eigenvalues/sum(pco.braycurtis$values$Eigenvalues))[2],1), "% Variation Explained")) + scale_color_discrete(name=curcond), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(ggplot(plot.pco.braycurtis, aes(x=Axis.2, y=Axis.3, color=plot.pco.braycurtis[,curcond])) + geom_point(alpha=0.5) + theme_bw() + theme(plot.title=element_text(size=10), axis.text=element_text(size=8), aspect.ratio = 1, axis.title=element_text(size=8,face="bold")) + ggtitle("Bray Curtis, PCo2 vs PCo3") + xlab(paste("PCo2:", round((100*pco.braycurtis$values$Eigenvalues/sum(pco.braycurtis$values$Eigenvalues))[2],1), "% Variation Explained")) + ylab(paste("PCo3:", round((100*pco.braycurtis$values$Eigenvalues/sum(pco.braycurtis$values$Eigenvalues))[3],1), "% Variation Explained")) + scale_color_discrete(name=curcond), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

    print(ggplot(plot.pco.weighted.UniFrac, aes(x=Axis.1, y=Axis.2, color=plot.pco.weighted.UniFrac[,curcond])) + geom_point(alpha=0.5) + theme_bw() + ggtitle("Weighted UniFrac, PCo1 vs PCo2") + theme(plot.title=element_text(size=10), axis.text=element_text(size=8), aspect.ratio = 1, axis.title=element_text(size=8,face="bold")) + xlab(paste("PCo1:", round((100*pco.weighted.UniFrac$values$Eigenvalues/sum(pco.weighted.UniFrac$values$Eigenvalues))[1],1), "% Variation Explained")) + ylab(paste("PCo2:", round((100*pco.weighted.UniFrac$values$Eigenvalues/sum(pco.weighted.UniFrac$values$Eigenvalues))[2],1), "% Variation Explained")) + scale_color_discrete(name=curcond), vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(ggplot(plot.pco.weighted.UniFrac, aes(x=Axis.2, y=Axis.3, color=plot.pco.weighted.UniFrac[,curcond])) + geom_point(alpha=0.5) + theme_bw() + ggtitle("Weighted UniFrac, PCo2 vs PCo3") + theme(plot.title=element_text(size=10), axis.text=element_text(size=8), aspect.ratio = 1, axis.title=element_text(size=8,face="bold")) + xlab(paste("PCo2:", round((100*pco.weighted.UniFrac$values$Eigenvalues/sum(pco.weighted.UniFrac$values$Eigenvalues))[2],1), "% Variation Explained")) + ylab(paste("PCo3:", round((100*pco.weighted.UniFrac$values$Eigenvalues/sum(pco.weighted.UniFrac$values$Eigenvalues))[3],1), "% Variation Explained")) + scale_color_discrete(name=curcond), vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

    print(ggplot(plot.pco.unweighted.UniFrac, aes(x=Axis.1, y=Axis.2, color=plot.pco.unweighted.UniFrac[,curcond])) + geom_point(alpha=0.5) + theme_bw() + ggtitle("Unweighted UniFrac, PCo1 vs PCo2") + theme(plot.title=element_text(size=10), axis.text=element_text(size=8), aspect.ratio = 1, axis.title=element_text(size=8,face="bold")) + xlab(paste("PCo1:", round((100*pco.unweighted.UniFrac$values$Eigenvalues/sum(pco.unweighted.UniFrac$values$Eigenvalues))[1],1), "% Variation Explained")) + ylab(paste("PCo2:", round((100*pco.unweighted.UniFrac$values$Eigenvalues/sum(pco.unweighted.UniFrac$values$Eigenvalues))[2],1), "% Variation Explained")) + scale_color_discrete(name=curcond), vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
    print(ggplot(plot.pco.unweighted.UniFrac, aes(x=Axis.2, y=Axis.3, color=plot.pco.unweighted.UniFrac[,curcond])) + geom_point(alpha=0.5) + theme_bw() + ggtitle("Unweighted UniFrac, PCo2 vs PCo3") + theme(plot.title=element_text(size=10), axis.text=element_text(size=8), aspect.ratio = 1, axis.title=element_text(size=8,face="bold")) + xlab(paste("PCo2:", round((100*pco.unweighted.UniFrac$values$Eigenvalues/sum(pco.unweighted.UniFrac$values$Eigenvalues))[2],1), "% Variation Explained")) + ylab(paste("PCo3:", round((100*pco.unweighted.UniFrac$values$Eigenvalues/sum(pco.unweighted.UniFrac$values$Eigenvalues))[3],1), "% Variation Explained")) + scale_color_discrete(name=curcond), vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
  }
  dev.off()

  print("--->Done")

}
