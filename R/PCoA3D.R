#' \code{PCoA3D} Generate Distance/Dissimilarities and plot PCoA in interactive 3D plot
#'
#' @description  UniFrac is implimented as per GUniFrac, Bray Curtis from Vegan, Jensen-Shannon divergence from textmineR, and PCoA from APE. Plotting via plotly's 3d plotting.
#'
#' @param METRIC Desired beta-diversity metric, options are Bray Curtis (braycurtis), Weighted UniFrac (weightedunifrac), UnWeighted UniFrac (unweightedunifrac), Jensen-Shannon diversgence (jsd).
#' @param METADATA Metadata table with variables to color or shape points by where sample is row name.
#' @param OTUTABLE Table of feature counts where Samples are columns, and IDs are row names.
#' @param TREE A phylogenetic tree of phylo class *only required for UniFrac Metrics
#' @param COLOR Metadata column to color samples by. Defaults to none.
#' @param SHAPE Metadata column to shape samples by. Defaults to none.
#' @param SUBSAMPLE Should table be subsampled? TRUE/FALSE (default: TRUE)
#' @param AXIS Which Axes should be plotted? Expects a numeric vector of length 3. defaults to  AXIS=c(1,2,3)
#' @param PALETTE a vector of COLORS to be used? defaults to the rainbow palette.

#'
#' @usage PCoA(METRIC="braycurtis", METADATA=experiment_metadata, OTUTABLE=otutable, TREE=ggtree, COLOR="Time", SHAPE="Group", SUBSAMPLE=TRUE, AXIS=c(1,2,3), PALETTE=c("red","green"))
#' @return Returns plot.
#' @export

PCoA3D<-function(METRIC,METADATA,OTUTABLE, TREE, COLOR, SHAPE, SUBSAMPLE, AXIS, PALETTE){

  if(missing(METRIC)){METRIC="braycurtis"}
  if(missing(PALETTE) & !missing(COLOR)){PALETTE=rainbow(length(unique(METADATA[[COLOR]])))}
  
  if(missing(TREE) & (METRIC=="weightedunifrac" | METRIC=="unweightedunifrac")){stop("ERROR: UniFrac requested, but tree not supplied...")}

  if(missing(SUBSAMPLE) || SUBSAMPLE==TRUE){
      DEPTH<-min(colSums(OTUTABLE))
      sub.OTUlevel<-Subsample.Table(OTUTABLE,DEPTH)
      if(!missing(TREE)){filtered.TREE<-prune_taxa(rownames(sub.OTUlevel), tree)}
  } else {
     sub.OTUlevel<-OTUTABLE
     if(!missing(TREE)){filtered.TREE<-TREE}
  }

if(missing(AXIS)){AXIS=c(1,2,3)}

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
    DISTMATRIX<-textmineR::CalcJSDivergence(sub.OTUlevel, by_rows=FALSE)
    METRIC="Jensen-Shannon Divergence"
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
