#' PCoA
#'
#' UniFrac is implimented as per Phyloseq, Bray Curtis from Vegan, Jensen-Shannon divergence from Phyloseq, and PCoA from APE. Plotting via ggplot2. Will also carry out ADONIS and return p and r2 in plot title. Distances calculated from proportions.
#'
#' @param METRIC Desired beta-diversity metric, options are Bray Curtis (braycurtis), Weighted UniFrac (weightedunifrac), UnWeighted UniFrac (unweightedunifrac), Jensen-Shannon diversgence (jsd).
#' @param METADATA Metadata Table with variables to color PCoAs by
#' @param FEATURES Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param TREE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names, only required for UniFrac Metrics
#' @param COLOR Metadata column to color samples by. Defaults to sample names from metadata row names.
#' @param SHAPE Metadata column to shape samples by. Defaults to none.
#' @param SUBSAMPLE Should table be subsampled? TRUE/FALSE (default: TRUE)
#' @param AXIS Which Axis should be plotted? Expects a numeric vector of length 2. defaults to Pco1 and Pco2 AXIS=c(1,2)
#' @param ADONIS Should ADONIS test be applied ? TRUE/FALSE (default TRUE)

#'
#' @usage PCoA(METRIC="braycurtis", METADATA=experiment_metadata, FEATURES=otutable, TREE=ggtree, COLOR="Time", SHAPE="Group", SUBSAMPLE=TRUE, AXIS=c(1,2))
#' @return Returns a ggplot
#' @export

PCoA<-function(METRIC,METADATA,FEATURES, TREE, COLOR, SHAPE, SUBSAMPLE, AXIS, ADONIS){

  if(missing(METRIC)){METRIC="braycurtis"}

  if(missing(TREE) & (METRIC=="weightedunifrac" | METRIC=="unweightedunifrac")){stop("ERROR: UniFrac requested, but tree not supplied...")}

  if(missing(ADONIS)){ADONIS=TRUE}

  if(missing(SUBSAMPLE) || SUBSAMPLE==TRUE){
    DEPTH<-min(colSums(FEATURES))
    sub.OTUlevel<-as.data.frame(Subsample.Table(FEATURES,DEPTH))
    sub.OTUlevel<-Make.Proportion(sub.OTUlevel)
    if(!missing(TREE)){filtered.TREE<-prune_taxa(rownames(sub.OTUlevel), TREE)}
  } else {
    sub.OTUlevel<-Make.Proportion(FEATURES)
    if(!missing(TREE)){filtered.TREE<-TREE}
  }

  if(missing(AXIS)){AXIS=c(1,2)}


  if(sum(!colnames(FEATURES) %in% rownames(METADATA))>0){stop("Metadata not available for all samples. Check metadata and feature table match.")}

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


  if (ADONIS==TRUE & !missing(COLOR)){
    if(sum(is.na(METADATA[,COLOR]))>0){stop(paste0("ADONIS could not be calculated with NA variable in: ", COLOR))}
    message(paste("Calculating ADONIS for", COLOR))
    RESULT<-vegan::adonis(DISTMATRIX ~ get(COLOR), METADATA, permutations=999)
    print(paste0(COLOR, "-> ADONIS P=", RESULT$aov.tab$`Pr(>F)`[1], " R2=", signif(RESULT$aov.tab$R2[1],3)))
  }


  if (ADONIS==TRUE & !missing(SHAPE)){
    if(sum(is.na(METADATA[,SHAPE]))>0){stop(paste0("ADONIS could not be calculated with NA variable in: ", SHAPE))}
    message(paste("Calculating ADONIS for", SHAPE))
    RESULT<-vegan::adonis(DISTMATRIX ~ get(SHAPE), METADATA, permutations=999)
    message(paste0(SHAPE, "-> ADONIS P=", RESULT$aov.tab$`Pr(>F)`[1], " R2=", signif(RESULT$aov.tab$R2[1],3)))
  }


  FINALPLOT<-PLOT %>% mutate(SampleID=Row.names) %>%
    mutate(X=get(paste0("Axis.", AXIS[1]))) %>%
    mutate(Y=get(paste0("Axis.", AXIS[2]))) %>%
    ggplot(aes(label=SampleID, x=X, y=Y)) +
    theme_bw() +
    theme(plot.title=element_text(size=10), axis.text=element_text(size=8), aspect.ratio = 1, axis.title=element_text(size=8,face="bold")) +
    ggtitle(METRIC) +
    xlab(paste0("PCo", AXIS[1],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[1]],1), "% Variation Explained")) +
    ylab(paste0("PCo", AXIS[2],": ", round((100*PCO$values$Eigenvalues/sum(PCO$values$Eigenvalues))[AXIS[2]],1), "% Variation Explained"))

  if(missing(COLOR) & missing(SHAPE)){FINALPLOT<-FINALPLOT +geom_point(alpha=0.5)+theme(legend.position="none")}
  else if (!missing(SHAPE) & !missing(COLOR)){FINALPLOT<-FINALPLOT+ geom_point(aes(shape=get(SHAPE), color=get(COLOR)), alpha=0.5) + labs(color=paste0(COLOR), shape=paste0(SHAPE)) }
  else if (!missing(SHAPE)){FINALPLOT<-FINALPLOT+ geom_point(aes(shape=get(SHAPE)), alpha=0.5) + labs(shape=paste0(SHAPE))}
  else if (!missing(COLOR)){FINALPLOT<-FINALPLOT+ geom_point(aes(color=get(COLOR)), alpha=0.5) + labs(color=paste0(COLOR))}

  return(FINALPLOT)
}
