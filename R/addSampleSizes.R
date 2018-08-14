#'
#' @title Add sample sizes to "size" compositions in a dataframe
#'
#' @description Function to add sample sizes to "size" compositions in a dataframe.
#'
#' @param dfr - input dataframe with composition data by "size" and other factors
#' @param id.size - name of column in dfr with "size" data
#' @param id.value - name of column in dfr with count (or abundance) values
#' @param verbose - flag to print diagnostic output to console
#'
#' @return dataframe with original size compositions and sample sizes in column "ss".
#'
#' @details Non-standard column names SHOULD NOT be back-quoted.
#'
#' @export
#'
addSampleSizes<-function(dfr,
                         id.size="size",
                         id.value=NULL,
                         verbose=FALSE){
  if (verbose) cat("\n\n#----Starting addSampleSizes\n");
  tmp<-dfr;
  #--Step 1: get factors (assumed to be all columns other than size and value)
  if (verbose) cat("\n#--Step 1-------------------------------------\n");
  id.facs<-names(tmp)[!(names(tmp) %in% c(id.size,id.value))];
  if (verbose) cat("#--id.facs =",paste(paste0("'",id.facs,"'",collapse=","),"\n"));

  #--Step 2: back-quote column names
  if (verbose) cat("\n#--Step 2-------------------------------------\n");
  idq.size<-paste0("`",id.size,"`");
  idq.value<-paste0("`",id.value,"`");
  if (!is.null(id.facs)) idq.facs<-paste0("`",id.facs,"`");

  #--Step 3: calculate sample sizes by factor combination
  cat("\n#--Step 3-------------------------------------\n");
  #----Example query:
  # qry<-"select
  #         fishery,area,`shell condition`,year,
  #         sum(count) as ss
  #       from tmp
  #       group by fishery,area,`shell condition`,year
  #       order by fishery,area,`shell condition`,year;";
  qry<-"select
          &&factors
          sum(&&idq.value) as ss
        from tmp
        &&groupby
        &&orderby";
  #determine strings for query
  str.factors<-"";
  str.groupby<-"";
  str.orderby<-"";
  if (!is.null(id.facs)) {
    str.factors<-paste0(idq.facs[1],",");
    str.groupby<-paste0("group by ",idq.facs[1]);
    str.orderby<-paste0("order by ",idq.facs[1]);
    if (length(idq.facs)>1){
      for (i in 2:length(idq.facs)) {
        str.factors<-paste0(str.factors,paste0(idq.facs[i],","));
        str.groupby<-paste0(str.groupby,", ",idq.facs[i]);
        str.orderby<-paste0(str.orderby,", ",idq.facs[i]);
      }
    }
  }
  #substitute
  qry<-gsub("&&factors",   str.factors, qry, fixed=TRUE);
  qry<-gsub("&&idq.value", idq.value,   qry, fixed=TRUE);
  qry<-gsub("&&groupby",   str.groupby, qry, fixed=TRUE);
  qry<-gsub("&&orderby",   str.orderby, qry, fixed=TRUE);
  if (verbose) cat("Query to calculate sample sizes:\n",qry,"\n");
  tmp1<-sqldf::sqldf(qry);

  #--Step 4: combine size compositions and sample size information
  cat("\n#--Step 4-------------------------------------\n");
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

  if (verbose) cat("#----Finished addSampleSizes\n");
  return(tmp2);
}
