#' \code{Microbiome.Barplot.R} Create a stacked barplot showing composition of samples.
#' @description Uses ggplot2 to create a stacked barplot, for example on phylum level abundances. The most abundant features (defaults to 10, based on rowMeans) will be plotted unless user specified. Anything of over 10 features will use default coloring which may be very difficult to interpret.
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param METADATA A Table of metadata where sample names are row names.
#' @param NTOPLOT A number of features to plot.
#' @param CATEGORY A Metadata category to block samples by (faceting via ggplot2)
#' @return Barplot
#' @export

Microbiome.Barplot<-function(OTUTABLE,METADATA, NTOPLOT, CATEGORY){

  if(missing(NTOPLOT) & nrow(OTUTABLE)>10){NTOPLOT=10}
  else if (missing(NTOPLOT)) {NTOPLOT=nrow(OTUTABLE)}

 OTUTABLE<-Make.Percent(OTUTABLE)
 OTUTABLE<-OTUTABLE[order(rowMeans(OTUTABLE), decreasing = T),]
 Remainder<-colSums(OTUTABLE[(NTOPLOT+1):nrow(OTUTABLE),])
 OTUTABLE<-rbind(OTUTABLE[1:NTOPLOT,], Remainder)

 forplot<-reshape2::melt(OTUTABLE)

 if(!missing(METADATA) & !missing(CATEGORY)){
  forplot<-merge(forplot, METADATA, by.x="Var2", by.y="row.names")
 }

 PLOT<-(ggplot(forplot, aes(x=Var2, y=value, fill=factor(Var1, levels = rev(levels(Var1)))))
        + geom_bar(stat="identity")
        + theme_classic()
        + ylab("% Abundance")
        + xlab("SampleID")
        + theme(axis.text.x =element_text(angle=45, hjust=1, size=6))
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
