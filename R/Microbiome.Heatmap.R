#' \code{Microbiome.Heatmap.R} Plot Heat map of microbiome data
#'
#' @description Creates a heat map based based on user provided table. Transforms to log2(percent), log10(percent), or CLR as requested. For plotting purposes, a prior of 0.01% is added to the log percent abundances. Row clustering of ggplot2 output is by upgma of euclidean dist. If no metadata passed, column clustering will be applied in the heatmap.2 output.
#'
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param METADATA Metadata file to be used for blocking.
#' @param TRANSFORM Method to transform data with for plotting, valid options are log2, log10, clr, percent, zscore or none (defaults to log10). None would be ideal if for example a list of fold changes was supplied. Zscore is calculated on clr.
#' @param NTOPFEATURES The N most abundance features to plot (defaults to all). Calculated by taking largest row sums of percentage data
#' @param CATEGORY (optional) Category to create separate blocks (optional, defaults to order of samples in otutable)
#' @param ORDER (optional) An order of levels in Category to use, (optional, defaults to order of unique(Category))
#' @param USERORDER (optional) A vector with user specified sample order for plot with NA for separators, overides metadata, category and order
#' @param PLOTMETHOD (optional) Can be plotted using heatmap.2 or ggplot2, defaults to ggplot2
#' @return Prints a heat map
#' @export

Microbiome.Heatmap<-function(OTUTABLE,METADATA, NTOPFEATURES, TRANSFORM, CATEGORY,ORDER,USERORDER, PLOTMETHOD){

if(missing(TRANSFORM)){TRANSFORM="log10"; print("No transform specified, using log10")}
if(missing(PLOTMETHOD)){PLOTMETHOD="ggplot2"; print("No plotting method specified, using ggplot2")}
if(missing(METADATA)){ print("No metadata given, using heatmap.2 row+column clustering")}

if(missing(METADATA)){
  DEND="both"
  COLCLUST=TRUE
} else {
  DEND="row"
  COLCLUST=FALSE
  }

if(!missing(CATEGORY)){METADATA[[CATEGORY]]<-factor(METADATA[[CATEGORY]])}
  tmp<-OTUTABLE

  if(TRANSFORM=="log10"){
    OTUTABLE<-log10(Make.Percent(OTUTABLE)+0.01)
  } else if(TRANSFORM=="log2"){
    OTUTABLE<-log2(Make.Percent(OTUTABLE)+0.01)
  } else if(TRANSFORM=="clr"){
    OTUTABLE<-Make.CLR(OTUTABLE)
  } else if(TRANSFORM=="percent"){
    OTUTABLE<-Make.Percent(OTUTABLE)
  } else if(TRANSFORM=="none"){
    print("No transformation applied to data")
  } else if(TRANSFORM=="zscore"){
    OTUTABLE<-t(apply(Make.CLR(OTUTABLE), 1, function(x){(x-mean(x))/sd(x)}))
  }

  print(paste("Minimum abundance:", min(OTUTABLE)))
  print(paste("Maximum abundance:", max(OTUTABLE)))

  if(!missing(NTOPFEATURES)){
    tmp<-cbind(rowSums(Make.Percent(tmp)), tmp)
    OTUTABLE<-OTUTABLE[rownames(tmp[order(tmp[,1], decreasing = T),])[1:NTOPFEATURES],]
    remove(tmp)
  }

roworder<-hclust(dist(OTUTABLE), method="average")
roworder<-roworder$labels[roworder$order]

if (!missing(USERORDER)){
    OTUTABLE<-OTUTABLE[,USERORDER]
  } else if(!missing(CATEGORY) && !missing(ORDER)){
    METADATA<-METADATA[which(METADATA[,CATEGORY] %in% ORDER),]
    NEWTABLE<-rep(NA, nrow(OTUTABLE))
    for (level in unique(ORDER)){
      NEWTABLE<-cbind(NEWTABLE, OTUTABLE[,rownames(METADATA[METADATA[[CATEGORY]]==level,])], rep(NA, nrow(OTUTABLE)), rep(NA, nrow(OTUTABLE)), rep(NA, nrow(OTUTABLE)))
    }
    OTUTABLE<-NEWTABLE
  } else if (!missing(CATEGORY)) {
    NEWTABLE<-rep(NA, nrow(OTUTABLE))
    for (level in levels(METADATA[,CATEGORY])){
      NEWTABLE<-cbind(NEWTABLE, OTUTABLE[,rownames(METADATA[METADATA[[CATEGORY]]==level,])], rep(NA, nrow(OTUTABLE)), rep(NA, nrow(OTUTABLE)), rep(NA, nrow(OTUTABLE)))
    }
    OTUTABLE<-NEWTABLE
  }



HTMP2<-heatmap.2(as.matrix(OTUTABLE),
            dendrogram=DEND,
            Rowv=TRUE,
            Colv=COLCLUST,
            revC=FALSE,
            trace="none",
            key.xlab=paste0(TRANSFORM, " relative abundance"),
            keysize="1",
            density.info="none",
            sepwidth=c(0.05,0.05),
            col=colorRampPalette(c("black","blue", "cyan","green","yellow","red"))(100),
            margins=c(12,20),
            symbreaks=FALSE,
            symkey=FALSE
  )



m.OTUTABLE<-melt(OTUTABLE, na.rm = T)
colnames(m.OTUTABLE)<-c("Taxa","Sample","Abundance")
P<-(ggplot(m.OTUTABLE, aes(x=Sample,y=factor(Taxa, levels=rev(roworder)), fill=Abundance, include.lowest = T))
    + geom_tile()
    + theme_bw()
    + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), axis.text.y = element_text(size=6))
    + labs(x="Sample", y="Taxa")
    + scale_fill_gradientn(colors=c("black","blue", "cyan","green","yellow","red"), name=paste0(TRANSFORM, " Abundance"))
)

if(!missing(CATEGORY)){for(i in 1:(length(unique(METADATA[[CATEGORY]]))-1)){
  level=levels(METADATA[[CATEGORY]])[i]
  xloc<-sum(METADATA[[CATEGORY]]==level)
  xloc<-xloc+sum(METADATA[[CATEGORY]] %in% levels(METADATA[[CATEGORY]])[0:(match(level, levels(METADATA[[CATEGORY]]))-1)])
    P<-P + geom_vline(xintercept=xloc+0.5, colour="white", linetype="dashed", size=1)
}

for (level in unique(METADATA[[CATEGORY]])){
  xloc<-sum(METADATA[[CATEGORY]]==level)/2
  xloc<-xloc+sum(METADATA[[CATEGORY]] %in% levels(METADATA[[CATEGORY]])[0:(match(level, levels(METADATA[[CATEGORY]]))-1)])
      P<-P + annotate("text", x=xloc, y=nrow(OTUTABLE)+1, label=level, size=4, vjust=1, hjust=0.5)
}}

if(PLOTMETHOD=="ggplot2"){ print(P)}
else if (PLOTMETHOD=="heatmap.2") {print(HTMP2)}

}
