#' Tidy.Convert
#'
#' Convert between table formats as needed (used behind the scenes) and is not fully implemented. TidyConvert Detects formatting of tables and interconverts between data.frame/data.table/matrix/tibble as needed.
#' When a format is encountered that does not use row names (data.table, tibble), the first column is always expected to be the key although right now a user supplied key is required.
#'
#' @param TABLE Table of feature/OTU/SV counts where Samples are columns, and feature IDs are row names.
#' @param KEY The column name for a unique identifier.
#' @return Converted table
#' @export
#'

TidyConvert.WhatAmI<-function(TABLE){
class(TABLE)[1]
}

TidyConvert.ToMatrix<-function(TABLE, KEY){
  if(TidyConvert.WhatAmI(TABLE)=="matrix"){return(TABLE)}
  else if(TidyConvert.WhatAmI(TABLE)=="data.frame"){return(TABLE %>% as.matrix())}
  else if(TidyConvert.WhatAmI(TABLE)=="data.table"){return(TABLE %>% as.data.frame() %>% column_to_rownames(KEY) %>% as.matrix())}
  else if(TidyConvert.WhatAmI(TABLE)=="tibble"){return(TABLE %>% as.data.frame() %>% column_to_rownames(KEY) %>% as.matrix())}
}

TidyConvert.ToDataFrame<-function(TABLE, KEY){
  if(TidyConvert.WhatAmI(TABLE)=="matrix"){return(TABLE %>% as.data.frame())}
  else if(TidyConvert.WhatAmI(TABLE)=="data.frame"){return(TABLE)}
  else if(TidyConvert.WhatAmI(TABLE)=="data.table"){return(TABLE %>% as.data.frame() %>% column_to_rownames(KEY))}
  else if(TidyConvert.WhatAmI(TABLE)=="tibble"){return(TABLE %>% as.data.frame() %>% column_to_rownames(KEY))}
}



TidyConvert.ToTibble<-function(TABLE, KEY){
  if(TidyConvert.WhatAmI(TABLE)=="matrix"){return(TABLE %>% as.data.frame %>% rownames_to_column(KEY) %>% as.tibble)}
  else if(TidyConvert.WhatAmI(TABLE)=="data.frame"){return(TABLE %>% rownames_to_column(KEY) %>% as.tibble)}
  else if(TidyConvert.WhatAmI(TABLE)=="data.table"){return(TABLE %>% as.tibble)}
  else if(TidyConvert.WhatAmI(TABLE)=="tibble"){return(TABLE)}
}

TidyConvert.ToDataTable<-function(TABLE, KEY){
  if(TidyConvert.WhatAmI(TABLE)=="matrix"){return(TABLE %>% as.data.frame %>% rownames_to_column(KEY) %>% as.data.table())}
  else if(TidyConvert.WhatAmI(TABLE)=="data.frame"){return(TABLE %>% rownames_to_column(KEY) %>% as.data.table)}
  else if(TidyConvert.WhatAmI(TABLE)=="data.table"){return(TABLE)}
  else if(TidyConvert.WhatAmI(TABLE)=="tibble"){return(TABLE %>% as.data.table())}
}
