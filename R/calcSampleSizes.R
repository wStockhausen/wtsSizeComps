#'
#' @title Calculate sample sizes from count or "size" composition data across a set of factors
#'
#' @description Function to sample sizes from count or "size" composition data across a set of factors.
#'
#' @param dfr - input dataframe with count/abundance data by size and other factors
#' @param id.value - name of column in dfr with count (or abundance) values
#' @param id.facs - character vector of column names used as factors in calculating the sample sizes
#' @param verbose - flag to print diagnostic output to console
#'
#' @return dataframe with sample sizes across a set of factors.
#'
#' @details Uses \code{sqldf::sqldf}. Note that sample sizes can be calculated from unnormalized size compositions,
#' as well as "raw" count data, using this function.
#'
#' @export
#'
calcSampleSizes<-function(dfr,
                          id.value="count",
                          id.facs="",
                          verbose=FALSE){
  #--normalize size compositions by norm.facs so
  #--they sum to 1 across remaining id.vars.
  if (verbose) cat("\n\n#----Starting calcSampleSizes\n");
  tmp<-dfr;
  if (verbose) cat("#--names(tmp): ",names(tmp),"\n");

  #----Step 1: backquote column names
  if (verbose) cat("#--Step 1-------------------------------------\n");
  idq.value<-paste0("`",id.value,"`");
  if (verbose) cat("#--id.faqs: ",id.facs,"\n");
  if (wtsUtilities::isBlankString(id.facs)) {
    idx<-!(names(tmp) %in% c(id.value));
    id.facs<-names(tmp)[idx];
    if (verbose) {
      cat("#--factor selection: ",idx,"\n");
      cat("#--id.faqs: ",id.facs,"\n");
    }
  }
  idq.facs<-paste0("`",id.facs,"`");
  if (verbose) cat("#--idq.facs: ",idq.facs,"\n");

  #----Step 2: calculate sample sizes
  if (verbose) cat("#--Step 2-------------------------------------\n");
  #----Example query:
  # qry<-"select
  #         year, fishery, area,
  #         sum(count) as ss
  #       from tmp
  #       group by year, fishery, area;";
  qry<-"select
          &&facs
          sum(&&idq.value) as ss
        from tmp
        &&groupby;";
  str.groupby  <-"";
  str.facs<-"";
  if (!wtsUtilities::isBlankString(id.facs)) {
    str<-paste(idq.facs,collapse=",");
    str.facs<-paste0(str,",");
    str.groupby  <-paste0("group by ",str);
  }
  qry<-gsub("&&facs", str.facs,qry,fixed=TRUE);
  qry<-gsub("&&idq.value",idq.value,   qry,fixed=TRUE);
  qry<-gsub("&&groupby",  str.groupby, qry,fixed=TRUE);
  if (verbose){cat("#--str.facs: '",str.facs,"'\n",sep="");}
  if (verbose){cat("#--str.groupby  : '",str.groupby,  "'\n",sep="");}
  if (verbose) cat("#--Query to calculate sample sizes:\n",qry,"\n");
  tmp1<-sqldf::sqldf(qry);

  if (verbose) cat("#--nrow(final) =",nrow(tmp2),"\n");

  if (verbose) cat("\n\n#----Finished calcSampleSizes\n");
  return(tmp1);
}

