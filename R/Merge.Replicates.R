#' Merge.Replicates
#'
#' Sum counts together for replicates (ex. pcr/extraction/sequencing run)
#'
#' @description Collapses technical replicates by summing replicate counts.
#'
#' @param FEATURES Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param METADATA Table of metadata with sample name as row name
#' @param CATEGORY A category used to dereplicate sequences , ex. "RepID"
#' @param VERBOSE Should report on number of replicates be given (T/F)
#' @return A named list of [[1]] new Feature/OTU/SV table and [[2]] collapsed metadata file
#' @usage x<-Merge.Replicates(svtable, metadata, "UniqueId")
#' newtable<-x$Features
#' newmeta<-x$Metadata
#'
#' @export

Merge.Replicates<-function(FEATURES, METADATA, CATEGORY, VERBOSE){
  if(missing(VERBOSE)){VERBOSE=T}
  if(VERBOSE==T){
    print(paste("Found",length(METADATA[,CATEGORY])-length(unique(METADATA[,CATEGORY])), "samples with replicates. Suming reads."))
  }

suppressMessages(
    TidyConvert.ToTibble(FEATURES, "Feature") %>%
    gather(-Feature, key="MBRSAMPLEIDS", value="Count") %>%
    left_join(TidyConvert.ToTibble(METADATA, "MBRSAMPLEIDS") %>%
    select(MBRSAMPLEIDS, CATEGORY)) %>%
    group_by_(as.name(CATEGORY), "Feature" ) %>%
    summarize(Count=sum(Count)) %>%
    spread(key=CATEGORY, value=Count, fill=0) %>%
    as.data.frame() %>%
    remove_rownames() %>%
    column_to_rownames("Feature") %>%
    as.matrix() -> FEATURES
    )

TidyConvert.ToTibble(METADATA, "MBRSAMPLEIDS") %>%
  filter(!duplicated(get(CATEGORY))) %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames(CATEGORY) %>%
  select(-MBRSAMPLEIDS) -> METADATA

return(list(Features=FEATURES,Metadata=METADATA))
}
