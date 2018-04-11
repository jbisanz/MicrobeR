#' Nice.Table
#'
#' @description Function to simplify the following command DT::datatable(TABLE, extensions='Buttons', options=list(pageLength=LENGTH, dom='bfrtip', buttons=c('copy','csv','excel', 'pdf') ))
#'
#' @param TABLE Any table to be printed
#' @param NROW (optional) number of rows to display per page, defaults to 10
#'
#' @export
#'
#'

Nice.Table<-function(TABLE, NROW){
if(missing(NROW))(NROW=10)
DTABLE<-DT::datatable(TABLE, extensions='Buttons', filter="top", options=list(pageLength=NROW, dom='Bfrtip', buttons=c('copy','csv','excel', 'pdf') ))
return(DTABLE)
}
