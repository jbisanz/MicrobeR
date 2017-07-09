#' \code{Summarize.Taxa} Make taxa summary tables similar to QIIME summarize_taxa.py
#' @description  Summarizes through the first 7 levels of taxonomy of OTUs/SVs wherein the levels are labeled as:
#' 1-Kingdom
#' 2-Phylum
#' 3-Class
#' 4-Order
#' 5-Family
#' 6-Genus
#' 7-Species
#' 8- Features with full taxonomy strings
#'
#' @param OTUTABLE Table of feature/OTU/SV COUNTS where Samples are columns, and IDs are row names
#' @param TAXONOMY Table of taxonomies to aggregate on as identified above wherein the column names correspond to the taxonomic level
#' @return A Named List of data frames of each level. Access by like TABLE$Genus or TABLE[[i]] where i is taxonomic level above.
#' @usage SummarizedTaxa<-Summarize.Taxa(SVtable, TaxTable)
#' @export


Summarize.Taxa<-function(OTUTABLE, TAXONOMY){

  
  if(sum(!rownames(OTUTABLE) %in% rownames(TAXONOMY))>0){stop("Mismatch between feature table and taxonomy table detected. There are features which are not in the taxonomy table. Check dimensions and row names.")}#check that all features in feature table are also in taxonomy table  
  
  TAXONOMY<-TAXONOMY[rownames(OTUTABLE),] #order the taxonomy table in the order as the taxonomy table
  
  taxastrings<-data.table(
  FeatureID=rownames(TAXONOMY),
  Kingdom=paste0(TAXONOMY$Kingdom),
  Phylum=paste(TAXONOMY$Kingdom, TAXONOMY$Phylum, sep=';'),
  Class=paste(TAXONOMY$Kingdom, TAXONOMY$Phylum, TAXONOMY$Class, sep=';'),
  Order=paste(TAXONOMY$Kingdom, TAXONOMY$Phylum, TAXONOMY$Class, TAXONOMY$Order, sep=';'),
  Family=paste(TAXONOMY$Kingdom, TAXONOMY$Phylum, TAXONOMY$Class, TAXONOMY$Order, TAXONOMY$Family, sep=';'),
  Genus=paste(TAXONOMY$Kingdom, TAXONOMY$Phylum, TAXONOMY$Class, TAXONOMY$Order, TAXONOMY$Family, TAXONOMY$Genus, sep=';'),
  Species=paste(TAXONOMY$Kingdom, TAXONOMY$Phylum, TAXONOMY$Class, TAXONOMY$Order, TAXONOMY$Family, TAXONOMY$Genus, TAXONOMY$Species, sep=';'),
  Feature=paste(rownames(TAXONOMY), TAXONOMY$Kingdom, TAXONOMY$Phylum, TAXONOMY$Class, TAXONOMY$Order, TAXONOMY$Family, TAXONOMY$Genus, TAXONOMY$Species, sep=';'),
  stringsAsFactors = F)
  taxastrings[is.na(taxastrings)]<-"NotAssigned"
  
  SUMMARIZED.TAXA=list()
  for (level in colnames(taxastrings)[2:8]){
    SUMMARIZED.TAXA[[level]]<-data.table(OTUTABLE)
    SUMMARIZED.TAXA[[level]]$TaxaString<-taxastrings[,get(level)] # was important that this was ordered above
    SUMMARIZED.TAXA[[level]]<-as.data.frame(SUMMARIZED.TAXA[[level]][, lapply(.SD, sum), by=TaxaString]) #from data.table vignette, SD is a data.table short hand for the subset of data
    rownames(SUMMARIZED.TAXA[[level]])<-SUMMARIZED.TAXA[[level]]$TaxaString
    SUMMARIZED.TAXA[[level]]$TaxaString<-NULL
  }
  SUMMARIZED.TAXA$Feature<-OTUTABLE
  rownames(SUMMARIZED.TAXA$Feature)<-taxastrings$Feature
  
  if(diff(range(lapply(SUMMARIZED.TAXA, sum)))!=0){stop("Reads went missing during summary. Debugging will be necesssary.")}
  
  return(SUMMARIZED.TAXA)
}
