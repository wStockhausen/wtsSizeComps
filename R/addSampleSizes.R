#'
#' @title Add sample sizes to "size" compositions in a dataframe
#'
#' @description Function to add sample sizes to "size" compositions in a dataframe.
#'
#' @param dfr - input dataframe with composition data by "size" and other factors
#' @param id.facs - character vector of column names used as factors in calculating the sample sizes
#' @param id.size - name of column in dfr with "size" data
#' @param id.value - name of column in dfr with count data
#' @param dfrSS - previously calculated dataframe with sample sizes by "other factors" (or NULL)
#' @param id.ss - name of column in dfrSS with sample sizes
#' @param verbose - flag to print diagnostic output to console
#'
#' @return dataframe with original size compositions and sample sizes in column "ss".
#'
#' @details Non-standard column names SHOULD NOT be back-quoted. If \code{dfrSS} is
#' provided, it will be used as the source of the sample size information (useful
#' for appending modified sample sizes to size compositions); otherwise, "raw"
#' sample sizes are calculated directly from \code{dfr}.
#'
#' @importFrom sqldf sqldf
#'
#' @export
#'
addSampleSizes<-function(dfr,
                         id.facs="",
                         id.size="size",
                         id.value="N",
                         dfrSS=NULL,
                         id.ss="ss",
                         verbose=FALSE){
  if (verbose) cat("\n\n#----Starting addSampleSizes\n");
  tmp<-dfr;
  #--Step 1: determine sample sizes
  if (is.null(dfrSS)){
    #--calculate sample sizes from size composition values
    #--Step 1a: get factors (assumed to be all columns other than size and value if id.facs is NULL)
    if (verbose) cat("\n#--Step 1-------------------------------------\n");
    if (is.null(id.facs)) id.facs<-names(tmp)[!(names(tmp) %in% c(id.size,id.value))];
    if (verbose) cat("#--id.facs =",paste(paste0("'",id.facs,"'",collapse=","),"\n"));

    #--Step 1b: calculate sample sizes
    tmp1<-calcSampleSizes(tmp,id.value=id.value,id.facs=id.facs,verbose=verbose);
  } else {
    #--copy sample sizes from input dataframe
    tmp1<-dfrSS;
    i<-which(names(tmp1)==id.ss);
    names(tmp1)[i]<-"ss";#--make sure sample size column is named 'ss'
  }

  #--Step 2: backquote column names
  idq.size<-paste0("`",id.size,"`");
  idq.value<-paste0("`",id.value,"`");
  idq.facs<-paste0("`",id.facs,"`");

  #--Step 3: combine size compositions and sample size information
  if (verbose) cat("\n#--Step 4-------------------------------------\n");
  #--Example query:
  # qry<-"select
  #         u.fishery,u.area,u.`shell condition`,u.year,
  #         u.ss,t.size,t.abundance
  #       from tmp1 as u left join tmp1 as t
  #       on
  #         t.fishery = u.fishery and
  #         t.area    = u.area and
  #         t.`shell condition` = u.`shell condition` and
  #         t.year    = u.year
  #       order by u.fishery,u.area,u.`shell condition`,u.year,u.size;";
  qry<-"select
          &&factors
          u.ss,t.&&idq.size,t.&&idq.value
        from tmp1 as u left join tmp as t
        &&onCond
        order by &&factors t.&&idq.size;";
  #determine strings for query
  str.factors<-"";
  str.onCond <-"";
  if (!is.null(id.facs)) {
    str.factors<-paste0("u.",idq.facs[1],",");
    str.onCond<-paste0("on ","t.",idq.facs[1],"=u.",idq.facs[1]);
    if (length(idq.facs)>1){
      for (i in 2:length(idq.facs)) {
        str.factors<-paste0(str.factors,paste0("u.",idq.facs[i],","));
        str.onCond<-paste0(str.onCond,"and \n",
                           paste0("t.",idq.facs[i],"=u.",idq.facs[i]));
      }
    }
  }
  #substitute
  qry<-gsub("&&factors",   str.factors, qry, fixed=TRUE);
  qry<-gsub("&&idq.size",  idq.size,    qry, fixed=TRUE);
  qry<-gsub("&&idq.value", idq.value,   qry, fixed=TRUE);
  qry<-gsub("&&onCond",    str.onCond,  qry, fixed=TRUE);
  if (verbose) cat("Query to join ss and size composition tables:\n",qry,"\n");
  tmp2<-sqldf::sqldf(qry);
  if (verbose) cat("#--nrow(final) =",nrow(tmp2),"\n");

  if (verbose) cat("#----Finished addSampleSizes\n");
  return(tmp2);
}
