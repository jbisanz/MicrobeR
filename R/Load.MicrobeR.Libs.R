#' \code{Load.MicrobeR.Libs} Import handy libraries
#' @description dada2
#' ggplot2
#' reshape2
#' vegan
#' ZIBR
#' ALDEx2
#' seqinr
#' phangorn
#' gplots
#' phyloseq
#' GUniFrac
#' grid
#' doBy
#' RColorBrewer
#' randomForest
#' ggord
#' plyr
#' dplyr
#' lme4
#' lmerTest
#' multcomp
#' DT
#' @export

Load.MicrobeR.Libs<-function(){
require(dada2, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(reshape2, quietly = TRUE)
require(vegan, quietly = TRUE)
require(ZIBR, quietly = TRUE)
require(ALDEx2, quietly = TRUE)
require(seqinr, quietly = TRUE)
require(phangorn, quietly = TRUE)
require(gplots, quietly = TRUE)
require(phyloseq, quietly = TRUE)
require(GUniFrac, quietly = TRUE)
require(grid, quietly = TRUE)
require(doBy, quietly = TRUE)
require(RColorBrewer, quietly = TRUE)
require(randomForest, quietly = TRUE)
require(ggord, quietly = TRUE)
require(plyr, quietly = TRUE)
require(dplyr, quietly = TRUE)
require(lme4, quietly = TRUE)
require(lmerTest, quietly = TRUE)
require(multcomp, quietly = TRUE)
require(DT, quietly = TRUE)
}
