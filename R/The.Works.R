#' The.Works
#'
#' Takes a feature/OTU/SV table and calculates an extensive set of metrics for downstream use. A number of versions of the normalized data are available from the returned object. See below for more information. For extremely large tables (>500x1000) this function is not reccomended as it will make numerous copies and have a significant memory footprint.
#'
#' @param FEATURES Table of feature/OTU/SV counts where Samples are columns, and IDs are row names.
#' @param TREE A phylogenetic tree for UniFrac and phylogenetic distances.
#' @param VERBOSE Should progress and metrics be printed to screen via message()? Default=TRUE
#' @usage Values<-TheWorks(MicrobeR.Demo.SVtable,MicrobeR.Demo.Tree)
#' @return A named list containing the following:
#' \cr.$Features$RawCounts- The original input file of raw counts
#' \cr.$Features$subsampled- A table of counts sampled to the lowest read count
#' \cr.$Features$Filtered- Counts where features present in less than 3 samples with less than 10 reads were removed
#' \cr.$Features$CLR- A CLR-normalized table (see Make.CLR)
#' \cr.$Features$Percent- The features as percent abundance (sum to 100)
#' \cr.$Features$Proportion- The feature as straight proportions (sum to 1)
#' \cr.$Features$PhILR- The Phylogenetic isometric log ratio transform (see package information for philr)
#' \cr.$Features$Tree$Raw- The input Tree
#' \cr.$Tree$... The corresponding pruned trees for the feature tables aboeve
#' \cr.$AlphaDiversity- A number of alpha diversity metrics in a data frame
#' \cr.$DistanceMatrices- A named list of distance matrices for UniFrac (w/uw), bray curtis, jsd, jaccard, PhILR euclidian, and CLR euclidian
#' \cr.$PCoA$Vectors- A named list numeric coordinates for plotting PCoA analysis
#' \cr.$PCoA$Scree- A table containing the percent variation explained for all calculated distance matrices to help in selecting the most appropriate metric.
#' @export

The.Works<-function(FEATURES,TREE, VERBOSE){
  if(missing(VERBOSE)){VERBOSE=TRUE}

  if(VERBOSE==T){message("Starting calculations at ", date())}
  WRX<-list()
  WRX$Features<-list()

  WRX$Features$RawCounts<-FEATURES
  WRX$Features$Subsampled<-Subsample.Table(FEATURES, VERBOSE=VERBOSE)
  WRX$Features$Filtered<-Confidence.Filter(FEATURES, MINSAMPS=3, MINREADS=10, VERBOSE=VERBOSE)
  WRX$Features$CLR<-Make.CLR(FEATURES, CZM=T)
  WRX$Features$Percent<-Make.Percent(FEATURES)
  WRX$Features$Proportion<-Make.Proportion(FEATURES)


  WRX$Tree<-list()
  WRX$Tree$Raw<-TREE
  WRX$Tree$Subsampled<-drop.tip(TREE, TREE$tip.label[!TREE$tip.label %in% rownames(WRX$Features$Subsampled)])
  WRX$Tree$PhILR<-makeNodeLabel(drop.tip(TREE, TREE$tip.label[!TREE$tip.label %in% rownames(WRX$Features$Filtered)]), method="number", prefix='n')
  WRX$Features$PhILR<-suppressMessages(
                                          t(
                                            philr(
                                            t(WRX$Features$Filtered+0.5),
                                             WRX$Tree$PhILR,
                                             part.weights='enorm.x.gm.counts',
                                             ilr.weights='blw.sqrt')
                                          )
                                        )

if(VERBOSE==T){message("Starting Alpha Diversity Calculations at ", date())}
  WRX$AlphaDiversity<-suppressMessages(
    data.frame(Shannon=diversity(WRX$Features$Subsampled, index="shannon", MARGIN=2)) %>% rownames_to_column("Sample") %>%
      left_join(data.frame(InverseSimpson=diversity(WRX$Features$Subsampled, index="invsimpson", MARGIN=2)) %>% rownames_to_column("Sample")) %>%
      left_join(data.frame(SpecRichness=specnumber(WRX$Features$Subsampled, MARGIN=2)) %>% rownames_to_column("Sample")) %>%
      left_join(pd(t(WRX$Features$Subsampled), WRX$Tree$Subsampled, include.root=T) %>% rownames_to_column("Sample")) %>%
      left_join(t(estimateR(t(WRX$Features$Subsampled))) %>% as.data.frame() %>% rownames_to_column("Sample")) %>%
      column_to_rownames("Sample") %>%
      select(-SR, S.obs) %>%
      rename(Chao1=S.chao1, Chao1.SE=se.chao1, ACE=S.ACE, ACE.SE=se.ACE)
    )


  if(VERBOSE==T){message("Starting distance calculations at ", date())}
  WRX$DistanceMatrices<-list()
  suppressMessages(WRX$DistanceMatrices$WeightedUniFrac<-UniFrac(phyloseq(otu_table(Make.Proportion(WRX$Features$Subsampled), taxa_are_rows = T), WRX$Tree$Subsampled), weighted=T))
  suppressMessages(WRX$DistanceMatrices$UnweightedUniFrac<-UniFrac(phyloseq(otu_table(Make.Proportion(WRX$Features$Subsampled), taxa_are_rows = T), WRX$Tree$Subsampled), weighted=F))
  WRX$DistanceMatrices$BrayCurtis<-distance(otu_table(Make.Proportion(WRX$Features$Subsampled), taxa_are_rows=T), method="bray", type="samples")
  WRX$DistanceMatrices$JSD<-distance(otu_table(Make.Proportion(WRX$Features$Subsampled), taxa_are_rows=T), method="jsd", type="samples")
  WRX$DistanceMatrices$Jaccard<-distance(otu_table(Make.Proportion(WRX$Features$Subsampled), taxa_are_rows=T), method="jaccard", type="samples")
  WRX$DistanceMatrices$CLREuclidian<-dist(t(WRX$Features$CLR), method="euclidian")
  WRX$DistanceMatrices$PhILREuclidian<-dist(t(WRX$Features$PhILR), method="euclidian")

  if(VERBOSE==T){message("Starting principal coordinate analyses calculations at ", date())}
  WRX$PCoA<-list()
  tmp<-lapply(WRX$DistanceMatrices, pcoa)
  WRX$PCoA$Vectors<-lapply(names(tmp), function(x) tmp[[x]]$vectors)
  names(WRX$PCoA$Vectors)<-names(tmp)
  WRX$PCoA$Scree<-do.call(rbind, lapply(names(tmp), function(x) data.frame(Metric=x, Axis=1:nrow(tmp[[x]]$values), PercentVariation=100*(tmp[[x]]$values$Relative_eig)) ) )
  rm(tmp)

  if(VERBOSE==T){message("Finished calculations at ", date())}
  return(WRX)
}
