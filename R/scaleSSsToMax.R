#'
#' @title Scale sample sizes, with a maximum allowed value
#'
#' @description Function to scale sample sizes, with a maximum allowed value.
#'
#' @param dfrSS - dataframe with sample sizes
#' @param colSS - name of column with sample sizes
#' @param sclSS - sample size to scale by
#' @param maxSS - max output sample size allowed
#'
#' @return copy of input dataframe, but with scaled sampled sizes in column \code{colSS}
#'
#' @details The sample sizes in column \code{colSS} are scaled such that
#'     dfrSS[[colSS]] = min(maxSS,maxSS/sclSS * dfrSS[[colSS]])
#'
#'@export
#'
scaleSSsToMax<-function(dfrSS,
                        colSS="ss",
                        sclSS=1,
                        maxSS=200){
  dfr<-dfrSS;
  scl<-maxSS/sclSS;
  dfr[[colSS]]<-scl*dfr[[colSS]];#scale observed sample sizes
  dfr[[colSS]][dfr[[colSS]]>maxSS]<-maxSS;#limit scaled ss's to maxSS
  return(dfr)
}
