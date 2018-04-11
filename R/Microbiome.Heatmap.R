#' Microbiome.Heatmap.R
#'
#' Creates a heat map based based on user provided table. Transforms to log2(percent), log10(percent), or CLR as requested. For plotting purposes, a prior of 0.01% is added to the log percent abundances.
#'
#' @param FEATURES Table of feature/OTU/SV counts where Samples are columns, and IDs are row names.
#' @param METADATA Metadata file to be used for blocking.
#' @param TRANSFORM Method to transform data with for plotting, valid options are log2, log10, clr, percent, zscore or none (defaults to log10). None would be ideal if for example a list of fold changes was supplied. Zscore is calculated on clr.
#' @param NTOPFEATURES The N most abundance features to plot (defaults to all). Calculated by taking largest row sums of percentage data
#' @param CATEGORY (optional) Category to create separate blocks (default: order of samples in otutable)
#' @param ROWCLUSTER (optional) How to order rows. Valid options are: UPGMA or abundance, default is UPGMA which is UPGMA clustering of euclidean distance of CLR-normalized counts
#' @return Prints a ggplot2 heatmap
#' @export

Microbiome.Heatmap<-function(FEATURES,METADATA, NTOPFEATURES, TRANSFORM, CATEGORY, ROWCLUSTER){

  if(missing(TRANSFORM)){TRANSFORM="log10"; message("No transform specified, using log10")}
  if(missing(METADATA)){ message("No metadata given, using heatmap.2 row+column clustering")}
  if(missing(ROWCLUSTER)){ROWCLUSTER="UPGMA"}

  if(ROWCLUSTER=="UPGMA"){
    message("Row clustering with UPGMA clustering.")
    suppressMessages(roworder<-hclust(dist(Make.CLR(FEATURES)), method="average"))
    roworder<-roworder$labels[roworder$order]
  } else if(ROWCLUSTER=="abundance"){
    message("Rows ranked on average abundance high to low.")
    roworder<-rowMeans(Make.Percent(FEATURES))
    roworder<-names(roworder[order(roworder, decreasing=TRUE)])
  }


  if(TRANSFORM=="log10"){
    FEATURES<-log10(Make.Percent(FEATURES)+0.01)
  } else if(TRANSFORM=="log2"){
    FEATURES<-log2(Make.Percent(FEATURES)+0.01)
  } else if(TRANSFORM=="clr"){
    FEATURES<-Make.CLR(FEATURES)
  } else if(TRANSFORM=="percent"){
    FEATURES<-Make.Percent(FEATURES)
  } else if(TRANSFORM=="none"){
    stop("No transformation applied, this is currently causing a bug, will be corrected in future.")
  } else if(TRANSFORM=="zscore"){
    FEATURES<-t(apply(Make.CLR(FEATURES), 1, function(x){(x-mean(x))/sd(x)}))
  }

  #print(paste("Minimum abundance:", min(FEATURES)))
  #print(paste("Maximum abundance:", max(FEATURES)))

  if(!missing(NTOPFEATURES)){
    FEATURES<-FEATURES[order(rowMeans(FEATURES), decreasing=TRUE)[1:NTOPFEATURES],]
  }

if(TidyConvert.WhatAmI(METADATA)=="data.frame" | TidyConvert.WhatAmI(METADATA)=="matrix") {METADATA<-TidyConvert.ToTibble(METADATA, "SampleID")}


  suppressMessages(
  FEATURES %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    gather(-Taxa, key="SampleID", value="Abundance") %>%
    left_join(METADATA) %>%
    mutate(Taxa=factor(Taxa, levels=rev(roworder))) -> FEATURES
  )

  P<-(ggplot(FEATURES, aes(x=SampleID,y=Taxa, fill=Abundance)) +
      geom_tile() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), axis.text.y = element_text(size=6)) +
      labs(x="Sample", y="Taxa") +
      scale_fill_gradientn(colors=c("black","blue", "cyan","green","yellow","red"), name=paste0(TRANSFORM, " Abundance"))
  )

  if(!missing(CATEGORY)){
    P<-P + facet_grid(~get(CATEGORY), margins=FALSE, drop=TRUE, scales = "free", space = "free")
  }

  return(P)

}
