#' \code{Generate.Test.Data} Generate random OTU table and metadata for testing purposes
#'
#' @description Samples random numbers between 1 and 1000 for read counts and creates a categorical variable for metadata into group A or B. Also gives phylogenetic tree by dropping tips from the gg 97% tree.
#'

#' @return Returns list of mock OTU table [[1]] and metadata [[2]]
#' @export
#'
#'

Generate.Test.Data<-function(){

Test.OTUs<-data.frame(row.names=sample(), SampleA=sample(1:1000,100),SampleB=sample(1:1000,100),SampleC=sample(1:1000,100),SampleD=sample(1:1000,100))
Test.Meta<-data.frame(row.names=c("SampleA","SampleB","SampleC","SampleD"), TestGroup=sample(c("GroupA","GroupB"), 4, replace=T), TestContinuous=sample(1:100, 4))


return(list(Test.OTUs, Test.Meta))
}
