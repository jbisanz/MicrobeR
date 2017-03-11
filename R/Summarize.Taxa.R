#' \code{Summarize.Taxa} Make taxa summary tables similar to QIIME summarize_taxa.py
#' @description  Summarizes through the first 7 levels of taxonomy of OTUs/SVs where in the levels are labeled as:
#' 1-Kingdom
#' 2-Phylum
#' 3-Class
#' 4-Order
#' 5-Family
#' 6-Genus
#' 7-Species
#'
#' @param OTUTABLE Table of feature/OTU/SV COUNTS where Samples are columns, and IDs are row names
#' @param TAXONOMY Table of taxonomies to aggregate on as identified above wherein the column names correspond to the taxonomic level
#' @return List of data frames of each level. Access by TABLE[[i]] where i is taxonomic level
#' @export


Summarize.Taxa<-function(OTUTABLE, TAXONOMY){
  SUMMARIZED.TAXA<-list()
  for (i in 1:7){
    print(paste("Now summarizing taxonomic level", i, "....."))
    if(i==1){TaxonomicLevel<-cbind(rownames(TAXONOMY),TAXONOMY$Kingdom)}
    if(i==2){TaxonomicLevel<-cbind(rownames(TAXONOMY),paste(TAXONOMY$Kingdom,TAXONOMY$Phylum,sep=";"))}
    if(i==3){TaxonomicLevel<-cbind(rownames(TAXONOMY),paste(TAXONOMY$Kingdom,TAXONOMY$Phylum,TAXONOMY$Class,sep=";"))}
    if(i==4){TaxonomicLevel<-cbind(rownames(TAXONOMY),paste(TAXONOMY$Kingdom,TAXONOMY$Phylum,TAXONOMY$Class,TAXONOMY$Order,sep=";"))}
    if(i==5){TaxonomicLevel<-cbind(rownames(TAXONOMY),paste(TAXONOMY$Kingdom,TAXONOMY$Phylum,TAXONOMY$Class,TAXONOMY$Order,TAXONOMY$Family ,sep=";"))}
    if(i==6){TaxonomicLevel<-cbind(rownames(TAXONOMY),paste(TAXONOMY$Kingdom,TAXONOMY$Phylum,TAXONOMY$Class,TAXONOMY$Order,TAXONOMY$Family,TAXONOMY$Genus,sep=";"))}
    if(i==7){TaxonomicLevel<-cbind(rownames(TAXONOMY),paste(TAXONOMY$Kingdom,TAXONOMY$Phylum,TAXONOMY$Class,TAXONOMY$Order,TAXONOMY$Family,TAXONOMY$Genus,TAXONOMY$Species,sep=";"))}
    rownames(TaxonomicLevel)<-rownames(TAXONOMY)
    colnames(TaxonomicLevel)<-c("temp","MergedTaxonomy")

    merged.OTUTABLE<-merge(TaxonomicLevel, OTUTABLE, all.x=FALSE, all.y=TRUE, by.x="temp", by.y="row.names")
    merged.OTUTABLE$temp<-NULL
    Summarized.OTUTABLE<-aggregate(. ~ MergedTaxonomy, data=merged.OTUTABLE, FUN=sum)
    rownames(Summarized.OTUTABLE)<-Summarized.OTUTABLE$MergedTaxonomy
    Summarized.OTUTABLE$MergedTaxonomy<-NULL
    SUMMARIZED.TAXA[[length(SUMMARIZED.TAXA)+1]] <- Summarized.OTUTABLE
  }
  taxonomy.OTUTABLE<-OTUTABLE #with taxonomy merged to OTUIDS and all items
  rownames(taxonomy.OTUTABLE)<-paste(row.names(taxonomy.OTUTABLE), TAXONOMY[row.names(taxonomy.OTUTABLE),]$mergedtaxonomy,sep="|")
  SUMMARIZED.TAXA[[length(SUMMARIZED.TAXA)+1]]<-taxonomy.OTUTABLE

  if((sum(SUMMARIZED.TAXA[[8]])!=sum(SUMMARIZED.TAXA[[7]])) &
     (sum(SUMMARIZED.TAXA[[8]])!=sum(SUMMARIZED.TAXA[[6]])) &
     (sum(SUMMARIZED.TAXA[[8]])!=sum(SUMMARIZED.TAXA[[5]])) &
     (sum(SUMMARIZED.TAXA[[8]])!=sum(SUMMARIZED.TAXA[[4]])) &
     (sum(SUMMARIZED.TAXA[[8]])!=sum(SUMMARIZED.TAXA[[3]])) &
     (sum(SUMMARIZED.TAXA[[8]])!=sum(SUMMARIZED.TAXA[[2]])) &
     (sum(SUMMARIZED.TAXA[[8]])!=sum(SUMMARIZED.TAXA[[1]]))){
    message("ERROR!!!!!!!!!!!!!!!, Reads went missing during taxa summary")
  }
  else{
    message("No reads went missing during summarizing, returning summarized info")
    return(SUMMARIZED.TAXA)
  }
} #Takes an OTU table and a master data frame of OTUs to taxonomy, returns a list of dataframes where the list index is the taxonomic level, the 8th level is OTU
