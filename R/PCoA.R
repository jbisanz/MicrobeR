#' \code{PCoA} Generate Distance/Dissimilarities and plot PCoA
#'
#' @description  UniFrac is implimented as per GUniFrac, Bray Curtis from Vegan, Jensen-Shannon divergence from textmineR, and PCoA from APE. Plotting via ggplot2. Will also carry out ADONIS and return p and r2 in plot title.
#'
#' @param METRIC Desired beta-diversity metric, options are Bray Curtis (braycurtis), Weighted UniFrac (weightedunifrac), UnWeighted UniFrac (unweightedunifrac), Jensen-Shannon diversgence (jsd).
#' @param METADATA Metadata Table with variables to color PCoAs by
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param TREE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names, only required for UniFrac Metrics
#' @param COLOR Metadata column to color samples by. Defaults to sample names from metadata row names.
#' @param SHAPE Metadata column to shape samples by. Defaults to none.
#' @param SUBSAMPLE Should table be subsampled? TRUE/FALSE (default: TRUE)
#' @param AXIS Which Axis should be plotted? Expects a numeric vector of length 2. defaults to Pco1 and Pco2 AXIS=c(1,2)
#'
#' @usage PCoA(METRIC="braycurtis", METADATA=experiment_metadata, OTUTABLE=otutable, TREE=ggtree, COLOR="Time", SHAPE="Group", SUBSAMPLE=TRUE, AXIS=c(1,2))
#' @return Returns a ggplot
#' @export

PCoA<-function(METRIC,METADATA,OTUTABLE, TREE, COLOR, SHAPE, CONDS, SUBSAMPLE, AXIS){

  if(missing(TREE) & (METRIC=="weightedunifrac" | METRIC=="unweightedunifrac")){stop("ERROR: UniFrac requested, but tree not supplied...")}


  if(missing(SUBSAMPLE) || SUBSAMPLE==TRUE){
      DEPTH<-min(colSums(OTUTABLE))
      sub.OTUlevel<-Subsample.Table(OTUTABLE,DEPTH)
      filtered.TREE<-prune_taxa(rownames(sub.OTUlevel), tree)
  } else {
     sub.OTUlevel<-OTUTABLE
     filtered.TREE<-TREE
  }

if(missing(AXIS)){AXIS=c(1,2)}

  if(METRIC=="unweightedunifrac"){
    print("Doing Unweighted UniFrac...")
    unifracs <- GUniFrac::GUniFrac(t(sub.OTUlevel),filtered.TREE, alpha=c(0, 0.5, 1))$unifracs
    DISTMATRIX <- unifracs[, , "d_UW"] # Unweighted UniFrac
    METRIC="Unweighted UniFrac"
  } else if(METRIC=="weightedunifrac") {
    print("Doing Weighted UniFrac...")
    unifracs <- GUniFrac::GUniFrac(t(sub.OTUlevel),filtered.TREE, alpha=c(0, 0.5, 1))$unifracs
    DISTMATRIX  <- unifracs[, , "d_1"] # Weighted UniFrac as per ?GUniFrac
    METRIC="Weighted UniFrac"
  } else if(METRIC=="braycurtis"){
    print("Doing Bray Curtis...")
    DISTMATRIX<-as.matrix(vegan::vegdist(t(sub.OTUlevel), method="bray", binary=FALSE, diag=FALSE, upper=FALSE))
    METRIC="Bray Curtis"
  } else if(METRIC=="jsd"){
    DISTMATRIX<-CalcJSDivergence(sub.OTUlevel, by_rows=FALSE)
    METRIC="Jensen-Shannon Divergence"
  }

  PCO<-ape::pcoa(DISTMATRIX)

  PLOT<-merge(as.data.frame(PCO$vectors), METADATA[rownames(PCO$vectors),], by="row.names", all=T )


  if(missing(COLOR)){COLORS=rownames(METADATA)} else{
    ADONIS<-vegan::adonis(DISTMATRIX ~ METADATA[rownames(DISTMATRIX),COLOR], permutations=999)
    print(paste0(COLOR, "-> ADONIS P=", ADONIS$aov.tab$`Pr(>F)`[1], " R2=", signif(ADONIS$aov.tab$R2[1],3)))
    COLORS<-METADATA[[COLOR]]

  }
  if(missing(SHAPE)){SHAPES=rep("Sample", nrow(METADATA))} else{
    ADONIS<-vegan::adonis(DISTMATRIX ~ METADATA[rownames(DISTMATRIX),SHAPE], permutations=999)
    print(paste0(SHAPE, "-> ADONIS P=", ADONIS$aov.tab$`Pr(>F)`[1], " R2=", signif(ADONIS$aov.tab$R2[1],3)))
    SHAPES<-METADATA[[SHAPE]]
  }

    FINALPLOT<-(ggplot(PLOT, aes(x=Axis.1, y=Axis.2, color=COLORS, shape=SHAPES))
          + geom_point(alpha=0.5) + theme_bw()
          + labs(color=COLOR, shape=SHAPE)
          + theme(plot.title=element_text(size=10), axis.text=element_text(size=8), aspect.ratio = 1, axis.title=element_text(size=8,face="bold"))
          + ggtitle(METRIC)
          + xlab(paste0("PCo", AXIS[1],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[1],1), "% Variation Explained"))
          + ylab(paste0("PCo", AXIS[2],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[2],1), "% Variation Explained"))
          )

return(FINALPLOT)
}
