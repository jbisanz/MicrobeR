#' Microbiome.Barplot.R
#'
#' @description Uses ggplot2 to create a stacked barplot, for example on phylum level abundances. The most abundant features (defaults to 10, based on rowMeans) will be plotted unless user specified. Anything of over 10 features will use default coloring which may be very difficult to interpret.
#' @param FEATURES Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param METADATA A Table of metadata where sample names are row names.
#' @param NTOPLOT A number of features to plot.
#' @param CATEGORY A Metadata category to block samples by (faceting via ggplot2)
#' @return Barplot
#' @usage Microbiome.Barplot(table,metadata,10, "group")
#' @export

Microbiome.Barplot<-function(FEATURES,METADATA, NTOPLOT, CATEGORY){

  if(missing(NTOPLOT) & nrow(FEATURES)>10){NTOPLOT=10}
  else if (missing(NTOPLOT)) {NTOPLOT=nrow(FEATURES)}

  FEATURES<-Make.Percent(FEATURES)
  FEATURES<-FEATURES[order(rowMeans(FEATURES), decreasing = T),]
 if(NTOPLOT<nrow(FEATURES)){ #if left over, is added to remainder
    Remainder<-colSums(FEATURES[(NTOPLOT+1):nrow(FEATURES),])
    FEATURES<-rbind(FEATURES[1:NTOPLOT,], Remainder)
 }

 forplot<-TidyConvert.ToTibble(FEATURES, "Taxa") %>% gather(-Taxa, key="SampleID", value="Abundance")
 forplot$Taxa<-factor(forplot$Taxa,levels=rev(unique(forplot$Taxa)))
 if(!missing(METADATA) & !missing(CATEGORY)){
   if(TidyConvert.WhatAmI(METADATA)=="data.frame" | TidyConvert.WhatAmI(METADATA)=="matrix") {METADATA<-TidyConvert.ToTibble(METADATA, "SampleID")}
  forplot<-left_join(forplot, METADATA, by="SampleID")
 }


 PLOT<-(ggplot(forplot, aes(x=SampleID, y=Abundance, fill=Taxa))
        + geom_bar(stat="identity")
        + theme_classic()
        + ylab("% Abundance")
        + xlab("SampleID")
        + theme(axis.text.x =element_text(angle=45, hjust=1, size=6))
        + theme(legend.text = element_text(size=6))
        )


if(NTOPLOT<=10){
  COLORS<-rev(c(
    "blue4",
    "olivedrab",
    "firebrick",
    "gold",
    "darkorchid",
    "steelblue2",
    "chartreuse1",
    "aquamarine",
    "yellow3",
    "coral",
    "grey"
    ))

 PLOT<-(PLOT + scale_fill_manual(values=COLORS, name=" "))
}

 if(!missing(CATEGORY)){
  PLOT<-( PLOT + facet_grid(~get(CATEGORY), margins=FALSE, drop=TRUE, scales = "free", space = "free") )

 }
  return(PLOT)
}
