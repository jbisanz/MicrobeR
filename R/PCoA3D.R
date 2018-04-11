#' PCoA3D
#'
#' @description  UniFrac is implimented as per Phyloseq, Bray Curtis from Vegan, Jensen-Shannon divergence from Phyloseq, and PCoA from APE. Plotting via plotly's 3d plotting.
#'
#' @param METRIC Desired beta-diversity metric, options are Bray Curtis (braycurtis), Weighted UniFrac (weightedunifrac), UnWeighted UniFrac (unweightedunifrac), Jensen-Shannon diversgence (jsd).
#' @param METADATA Metadata table with variables to color or shape points by where sample is row name.
#' @param FEATURES Table of feature counts where Samples are columns, and IDs are row names.
#' @param TREE A phylogenetic tree of phylo class *only required for UniFrac Metrics
#' @param COLOR Metadata column to color samples by. Defaults to none.
#' @param SHAPE Metadata column to shape samples by. Defaults to none.
#' @param SUBSAMPLE Should table be subsampled? TRUE/FALSE (default: TRUE)
#' @param AXIS Which PCs should be plotted? Expects a numeric vector of length 3. defaults to  AXIS=c(1,2,3)
#' @param PALETTE a vector of COLORS to be used? defaults to the rainbow palette.
#'
#' @usage PCoA(METRIC="braycurtis", METADATA=experiment_metadata, FEATURES=otutable, TREE=ggtree, COLOR="Time", SHAPE="Group", SUBSAMPLE=TRUE, AXIS=c(1,2,3), PALETTE=c("red","green"))
#' @return Returns plot.
#' @export

PCoA3D<-function(METRIC,METADATA,FEATURES, TREE, COLOR, SHAPE, SUBSAMPLE, AXIS, PALETTE){

  if(missing(METRIC)){METRIC="braycurtis"}
  if(missing(PALETTE) & !missing(COLOR)){PALETTE=rainbow(length(unique(METADATA[[COLOR]])))}

  if(missing(TREE) & (METRIC=="weightedunifrac" | METRIC=="unweightedunifrac")){stop("ERROR: UniFrac requested, but tree not supplied...")}

  if(missing(SUBSAMPLE) || SUBSAMPLE==TRUE){
    DEPTH<-min(colSums(FEATURES))
    sub.OTUlevel<-as.data.frame(Subsample.Table(FEATURES,DEPTH))
    sub.OTUlevel<-Make.Proportion(sub.OTUlevel)
    if(!missing(TREE)){filtered.TREE<-prune_taxa(rownames(sub.OTUlevel), TREE)}
  } else {
    sub.OTUlevel<-Make.Proportion(FEATURES)
    if(!missing(TREE)){filtered.TREE<-TREE}
  }

  if(missing(AXIS)){AXIS=c(1,2,3)}

  if(METRIC=="unweightedunifrac"){
    message("Doing Unweighted UniFrac...")
    DISTMATRIX <- UniFrac(phyloseq(otu_table(sub.OTUlevel, taxa_are_rows = T), filtered.TREE), weighted=F)
    METRIC="Unweighted UniFrac"
  } else if(METRIC=="weightedunifrac") {
    message("Doing Weighted UniFrac...")
    DISTMATRIX <- UniFrac(phyloseq(otu_table(sub.OTUlevel, taxa_are_rows = T), filtered.TREE), weighted=T)
    METRIC="Weighted UniFrac"
  } else if(METRIC=="braycurtis"){
    message("Doing Bray Curtis...")
    DISTMATRIX<-vegan::vegdist(t(sub.OTUlevel), method="bray", binary=FALSE, diag=FALSE, upper=FALSE)
    METRIC="Bray Curtis"
  } else if(METRIC=="jsd"){
    DISTMATRIX<-distance(otu_table(sub.OTUlevel, taxa_are_rows=T), method="jsd", type="samples")
    METRIC="Jensen-Shannon Divergence"
  } else {
    stop("Metric choice not recognized")
  }

  PCO<-ape::pcoa(DISTMATRIX)

  PLOT<-merge(as.data.frame(PCO$vectors), METADATA[rownames(PCO$vectors),], by="row.names", all=T )



  if(missing(COLOR) & missing(SHAPE)){
    FINALPLOT<-(plot_ly(data = PLOT, type = "scatter3d", x = ~get(paste0("Axis.",AXIS[1])), y = ~get(paste0("Axis.",AXIS[2])), z = ~get(paste0("Axis.",AXIS[3])),
                        mode = "markers") %>%
                  layout(scene = list(
                    xaxis = list(title = paste0("PCo", AXIS[1],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[1]],1), "%")),
                    yaxis = list(title = paste0("PCo", AXIS[2],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[2]],1), "%")),
                    zaxis = list(title = paste0("PCo", AXIS[3],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[3]],1), "%"))
                  )))
  } else if(!missing(COLOR) & missing(SHAPE)){
    FINALPLOT<-(plot_ly(data = PLOT, type = "scatter3d", x = ~get(paste0("Axis.",AXIS[1])), y = ~get(paste0("Axis.",AXIS[2])), z = ~get(paste0("Axis.",AXIS[3])),
                        color = ~get(COLOR), colors = PALETTE, mode = "markers") %>%
                  layout(scene = list(
                    xaxis = list(title = paste0("PCo", AXIS[1],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[1]],1), "%")),
                    yaxis = list(title = paste0("PCo", AXIS[2],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[2]],1), "%")),
                    zaxis = list(title = paste0("PCo", AXIS[3],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[3]],1), "%"))
                  )))
  }else if(missing(COLOR) & !missing(SHAPE)){
    FINALPLOT<-(plot_ly(data = PLOT, type = "scatter3d", x = ~get(paste0("Axis.",AXIS[1])), y = ~get(paste0("Axis.",AXIS[2])), z = ~get(paste0("Axis.",AXIS[3])),
                        symbol=~get(SHAPE), symbols=seq(15,20,1), mode = "markers") %>%
                  layout(scene = list(
                    xaxis = list(title = paste0("PCo", AXIS[1],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[1]],1), "%")),
                    yaxis = list(title = paste0("PCo", AXIS[2],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[2]],1), "%")),
                    zaxis = list(title = paste0("PCo", AXIS[3],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[3]],1), "%"))
                  )))
  } else if(!missing(COLOR) & !missing(SHAPE)){
    FINALPLOT<-(plot_ly(data = PLOT, type = "scatter3d", x = ~get(paste0("Axis.",AXIS[1])), y = ~get(paste0("Axis.",AXIS[2])), z = ~get(paste0("Axis.",AXIS[3])),
                        color = ~get(COLOR), colors = PALETTE, symbol=~get(SHAPE), symbols=seq(15,20,1), mode = "markers") %>%
                  layout(scene = list(
                    xaxis = list(title = paste0("PCo", AXIS[1],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[1]],1), "%")),
                    yaxis = list(title = paste0("PCo", AXIS[2],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[2]],1), "%")),
                    zaxis = list(title = paste0("PCo", AXIS[3],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[3]],1), "%"))
                  )))
  }

  return(FINALPLOT)
}
