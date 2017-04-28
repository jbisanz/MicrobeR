#' \code{Generate.Test.Data} Generate random OTU table and metadata for testing purposes.
#'
#' @description Creates a mock dataset based on a user defined number of samples and matching metadata for testing purposes. Returns a categorical variable called TestGroup and a continuous number called TestContinuous.
#'
#' @return Returns list of mock OTU table [[1]] and metadata [[2]]
#' @export
#'
#'

Generate.Test.Data<-function(NSAMP, NFEATURE){
Test.OTUs<-data.frame(matrix(ncol = NSAMP, nrow = NFEATURE))
Test.OTUs<-apply(Test.OTUs, 2, function(x){sample(1000, NFEATURE)})
colnames(Test.OTUs)=paste0("Sample_", seq(1,NSAMP,1))
rownames(Test.OTUs)<-paste0("OTU_", seq(1,NFEATURE,1))


Test.Meta<-data.frame(row.names=paste0("Sample_", seq(1,NSAMP,1)), TestGroup=sample(c("GroupA","GroupB"), NSAMP, replace=T), TestContinuous=sample(1:100, NSAMP))


return(list(Test.OTUs, Test.Meta))
}
