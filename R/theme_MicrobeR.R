#' \code{theme_MicrobeR} A ggplot2 theme
#'
#' @export
#' 
 
theme_MicrobeR <- function () { 
  theme_classic(base_size=9, base_family="Helvetica") +
    theme(panel.border = element_rect(color="black", size=1, fill=NA)) +
    theme(axis.line = element_blank(), strip.background = element_blank())
}
mean_sd=function(x){data.frame(y=mean(x), ymin=mean(x)-sd(x), ymax=mean(x)+sd(x))}