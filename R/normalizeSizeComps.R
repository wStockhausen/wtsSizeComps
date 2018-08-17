#'
#' @title Normalize "size" compositions across a set of factors
#'
#' @description Function to create a dataframe of "size" compositions normalized (summing to 1) across a set of factors.
#'
#' @param dfrZ - input dataframe with count/abundance data by size and other factors
#' @param id.size - name of column in dfr with "size" data
#' @param id.value - name of column in dfr with count (or abundance) values
#' @param id.facs - character vector of column names used as factors in calculating the size compositions
#' @param id.normfacs - character vector of column names of factors to normalize the size compositions across (a subset of id.facs)
#' @param verbose - flag to print diagnostic output to console
#'
#' @return dataframe with "size" compositions normalized to sum to 1 across a set of factors.
#'
#' @details Uses \code{sqldf::sqldf}.
#'
#' @export
#'
normalizeSizeComps<-function(dfr,
                             id.size="size",
                             id.value="count",
                             id.facs="",
                             id.normfacs="",
                             verbose=FALSE){
  #--normalize size compositions by norm.vars so
  #--they sum to 1 across remaining id.vars.
  if (verbose) cat("\n\n#----Starting normalizeSizeComps\n");
  tmp<-dfr;
  if (verbose) cat("#--names(tmp): ",names(tmp),"\n");

  #----Step 1: backquote column names
  if (verbose) cat("#--Step 1-------------------------------------\n");
  idq.size<-paste0("`",id.size,"`");
  idq.value<-paste0("`",id.value,"`");
  if (verbose) cat("#--id.faqs: ",id.facs,"\n");
  if (wtsUtilities::isBlankString(id.facs)) {
    idx<-!(names(tmp) %in% c(id.size,id.value));
    id.facs<-names(tmp)[idx];
    if (verbose) {
      cat("#--factor selection: ",idx,"\n");
      cat("#--id.faqs: ",id.facs,"\n");
    }
  }
  idq.facs<-paste0("`",id.facs,"`");
  if (verbose) cat("#--idq.facs: ",idq.facs,"\n");
  if (!wtsUtilities::isBlankString(id.normfacs)) {
    idq.normfacs<-paste0("`",id.normfacs,"`");
    if (verbose) cat("#--idq.normfacs: ",idq.normfacs,"\n");
  }

  #----Step 2: calculate totals to normalize by
  if (verbose) cat("#--Step 2-------------------------------------\n");
  #----Example query:
  # qry<-"select
  #         year, fishery, area,
  #         sum(count) as totcount
  #       from tmp
  #       group by year, fishery, area;";
  qry<-"select
          &&normfacs
          sum(&&idq.value) as total
        from tmp
        &&groupby;";
  str.groupby  <-"";
  str.norm.vars<-"";
  if (!wtsUtilities::isBlankString(id.normfacs)) {
    str<-paste(idq.normfacs,collapse=",");
    str.normfacs<-paste0(str,",");
    str.groupby  <-paste0("group by ",str);
  }
  qry<-gsub("&&normfacs", str.normfacs,qry,fixed=TRUE);
  qry<-gsub("&&idq.value",idq.value,   qry,fixed=TRUE);
  qry<-gsub("&&groupby",  str.groupby, qry,fixed=TRUE);
  if (verbose){cat("#--str.normfacs: '",str.normfacs,"'\n",sep="");}
  if (verbose){cat("#--str.groupby  : '",str.groupby,  "'\n",sep="");}
  if (verbose) cat("#--Query to calculate total counts:\n",qry,"\n");
  tmp1<-sqldf::sqldf(qry);

  #----Step 3: normalize compositions
  if (verbose) cat("#--Step 3-------------------------------------\n");
  #----Example query:
  # qry<-"select
  #         t.year, t.fishery, t.area,t.`shell condition`,
  #         t1.total as ss,
  #         t.size,
  #         t.count/t1.total as fraction
  #       from tmp as t, tmp1 as t1
  #       where
  #         t.year    = t1.year and
  #         t.fishery = t1.fishery and
  #         t.area    = t1.area;";
  qry<-"select
          &&idq.vars
          t.&&idq.size,
          t.&&idq.value/t1.total as fraction
        from tmp as t, tmp1 as t1
        &&where
        &&orderby;";
  str.id.vars<-"";
  if (!wtsUtilities::isBlankString(id.facs)) {
    str.id.vars<-paste0(paste("t",idq.facs,sep=".",collapse=","),",");
    str.orderby<-paste0("order by\n",paste("t",idq.facs,sep=".",collapse=","));
  }
  str.where<-"";
  if (!wtsUtilities::isBlankString(idq.normfacs)) {
    str.where<-paste0("where\n",paste0("t.",idq.normfacs[1],"=t1.",idq.normfacs[1]));
    nNVs<-length(idq.normfacs);
    if (nNVs>1){ for (i in 2:nNVs) {str.where<-paste0(str.where," and \n",paste0("t.",idq.normfacs[i],"=t1.",idq.normfacs[i]));}}
  }
  qry<-gsub("&&idq.vars", str.id.vars, qry, fixed=TRUE);
  qry<-gsub("&&idq.size", idq.size,    qry,fixed=TRUE);
  qry<-gsub("&&idq.value",idq.value,   qry, fixed=TRUE);
  qry<-gsub("&&where",    str.where,   qry, fixed=TRUE);
  qry<-gsub("&&orderby",  str.orderby, qry, fixed=TRUE);
  if (verbose){cat("#--str.id.vars: '",str.id.vars,"'\n",sep="");}
  if (verbose){cat("#--str.where  : '",str.where,  "'\n",sep="");}
  if (verbose){cat("#--str.orderby: '",str.orderby,"'\n",sep="");}
  if (verbose) cat("Query to calculate normalized compositions:\n",qry,"\n");
  tmp2<-sqldf::sqldf(qry);
  #--set NA's to zeros
  idx<-is.na(tmp2[["fraction"]]);
  tmp2[idx,"fraction"]<-0;
  if (verbose) cat("#--setting",sum(idx),"NA values to 0\n");
  if (verbose) cat("#--nrow(final) =",nrow(tmp2),"\n");

  if (verbose) cat("\n\n#----Finished normalizeSizeComps\n");
  return(tmp2);
}

