
#'
#' @title Determine the unique factor and size combinations from a dataframe of "size" compositions
#'
#' @description Function to determine the unique factor and size combinations from a dataframe of "size" compositions.
#'
#' @param dfr - input dataframe with composition data by "size" and other factors
#' @param id.size - name of column in dfr with "size" data
#' @param id.value - name of column in dfr with count (or abundance) values
#' @param id.facs - character vector of factor column names
#' @param expandToAllFactorCombos - flag to expand to all factor combinations (TRUE) or only observed combinations (FALSE)
#' @param verbose - flag to print diagnostic output to console
#'
#' @return dataframe with compositions by size (and possibly other factors) expanded to all sizes.
#'
#' @details If id.facs is blank, all columns other than the size and value columns
#' are assumed to be factor columns.
#'
#' Non-standard column names SHOULD NOT be backquoted.
#'
#' @export
#'
getUniqueFZs<-function(dfr,
                       id.size="size",
                       id.value=NULL,
                       id.facs="",
                       expandToAllFactorCombos=FALSE,
                       verbose=verbose){
  tmp<-dfr;
  #--Step 1: determine unique factor x size combinations
  if (verbose) cat("\n#--Step 1-------------------------------------\n");
  if (wtsUtilities::isBlankString(id.facs)) id.facs<-names(tmp)[!(names(tmp) %in% c(id.size,id.value))];
  if (verbose) cat("#--id.facs =",paste(paste0("'",id.facs,"'",collapse=","),"\n"));
  uZs<-data.frame(size=cutpts[1:(nCPs-1)]);
  tmp2<-unique(tmp[,id.facs]);
  uFs<-tmp2;
  #back-quote column names
  idq.size<-paste0("`",id.size,"`");
  idq.value<-paste0("`",id.value,"`");
  idq.facs<-paste0("`",id.facs,"`");
  if (expandToAllFactorCombos) {
    #--Example query:
    # qry<-"select * from
    #         (select distinct fishery from tmp2),
    #         (select distinct area from tmp2),
    #         (select distinct `shell condition` from tmp2),
    #         (select distinct year from tmp2);";
    qry<-"select * from
            &&str.factors;";
    str.factors<-paste0("(select distinct ",idq.facs[1]," from tmp2)");
    if (length(idq.facs)>1){
      for (i in 2:length(idq.facs)) {
        str.factors<-paste0(str.factors,
                            paste0(",\n(select distinct ",idq.facs[i]," from tmp2)"));
      }
    }
    qry<-gsub("&&str.factors",str.factors,qry, fixed=TRUE);
    if (verbose) cat("Query to calculate factor combinations:\n",qry,"\n");
    uFs <- sqldf::sqldf(qry);
  }
  uFZs<-sqldf::sqldf("select * from uFs,uZs;");
  return(uFZs);
}

#'
#' @title Re-bin "size" compositions in a dataframe to all given "size" values (and possibly all factor level combinations)
#'
#' @description Function to re-bin "size" compositions in a dataframe to all given "size" values (and possibly all factor level combinations).
#'
#' @param dfr - input dataframe with composition data by "size" and other factors
#' @param id.size - name of column in dfr with "size" data
#' @param id.value - name of column in dfr with count (or abundance) values
#' @param id.ss - name of column in dfr with sample sizes
#' @param id.facs - character vector of factor column names
#' @param cutpts - vector of cutpoints to use for re-binning size compositions
#' @param truncate.low -flag to truncate size compositions below the minimum size (TRUE) or include in smallest size bin (FALSE)
#' @param truncate.high - flag to truncate size compositions above the maximum size (TRUE) or include in largest size bin (FALSE)
#' @param expandToAllFactorCombos - flag to expand to all factor combinations (TRUE) or only observed combinations (FALSE)
#' @param verbose - flag to print diagnostic output to console
#'
#' @return dataframe with compositions by size (and possibly other factors) expanded to all sizes by either unique or all factor combinations.
#'
#' @details Non-standard column names SHOULD NOT be back-quoted.
#'
#' @export
#'
rebinSizeComps<-function(dfr,
                         id.size="size",
                         id.value=NULL,
                         id.ss="ss",
                         id.facs="",
                         cutpts=NULL,
                         truncate.low=TRUE,
                         truncate.high=FALSE,
                         expandToAllFactorCombos=TRUE,
                         verbose=FALSE){
  tmp<-dfr;
  #--Step 1: rebin compositions (if necessary) to new size bins
  if (verbose) cat("#--Step 1-------------------------------------\n");
  if (is.null(cutpts)) {
    cutpts<-data.frame(size=min(tmp$size,na.RM=TRUE):(max(tmp$size,na.rm=TRUE)+1));
    newcutpts<-cutpts;
  } else {
    nCPs<-length(cutpts);
    newcutpts<-cutpts;#make copy to apply truncation correctly
    if (!truncate.low ) newcutpts[1]<-0;
    if (!truncate.high) newcutpts[nCPs]<-Inf;
  }
  if (verbose) cat("#--cutpoints:",cutpts,"\n");
  cuts<-cut(tmp[[id.size]],newcutpts,right=FALSE,labels=FALSE);#make cuts based on new bins adjusted for truncation
  tmp[[id.size]]<-cutpts[cuts];                                #assign to original bins using cutpoints
  tmp<-tmp[!is.na(tmp[[id.size]]),];#drop truncated data
  #----Example query:
  # qry<-"select
  #         fishery,area,`shell condition`,year,ss,size,
  #         sum(abundance) as abundance
  #       from tmp
  #       group by fishery,area,`shell condition`,year,ss,size
  #       order by fishery,area,`shell condition`,year,ss,size;";
  qry<-"select
          &&factors &&id.size,
          sum(&&idq.value) as &&idq.value
        from tmp
        group by &&factors &&idq.size
        order by &&factors &&idq.size;";
  #back-quote column names
  idq.size<-paste0("`",id.size,"`");
  idq.value<-paste0("`",id.value,"`");
  if (wtsUtilities::isBlankString(id.facs)) id.facs<-names(tmp)[!(names(tmp) %in% c(id.size,id.value))];
  idq.facs<-paste0("`",id.facs,"`");
  #determine factor string for query
  str.factors<-"";
  if (length(idq.facs)>0){
    str.factors<-paste0(idq.facs[1],",");
    if (length(idq.facs)>1){
      for (i in 2:length(idq.facs)) {
        str.factors<-paste0(str.factors,paste0(idq.facs[i],","));
      }
    }
  }
  #substitute
  qry<-gsub("&&factors",   str.factors, qry, fixed=TRUE);
  qry<-gsub("&&idq.size",  idq.size,    qry, fixed=TRUE);
  qry<-gsub("&&idq.value", idq.value,   qry, fixed=TRUE);
  if (verbose) cat("Query to rebin size comps:\n",qry,"\n");
  tmp1<-sqldf::sqldf(qry);

  #--Step 2: determine unique factor x size combinations
  if (verbose) cat("\n#--Step 2-------------------------------------\n");
    uFZs<-getUniqueFZs(tmp,
                       id.size=id.size,
                       id.value=id.value,
                       id.facs=id.facs,
                       expandToAllFactorCombos=expandToAllFactorCombos,
                       verbose=verbose)

  #--Step 3 expand to all sizes
  if (verbose) cat("\n#--Step 3-------------------------------------\n");
  #--Example query:
  # qry<-"select
  #         u.fishery,u.area,u.`shell condition`,u.year,
  #         u.size,t.abundance
  #       from uFZs as u left join tmp as t
  #       on
  #         t.fishery = u.fishery and
  #         t.area    = u.area and
  #         t.`shell condition` = u.`shell condition` and
  #         t.year    = u.year and
  #         t.size    = u.size
  #       order by u.fishery,u.area,u.`shell condition`,u.year,u.size;";
  qry<-"select
          &&factors
          u.&&id.size,t.&&idq.value
        from uFZs as u left join tmp as t
        on
          &&onCond
          t.&&idq.size    = u.&&idq.size
        order by &&factors u.&&idq.size;";
  str.factors<-"";
  str.onCond <-"";
  if (length(idq.facs)>0){
    str.factors<-paste0("u.",idq.facs[1],",");
    str.onCond <-paste0("t.",idq.facs[1],"=u.",idq.facs[1]," and \n");
    if (length(idq.facs)>1){
      for (i in 2:length(idq.facs)) {
        str.factors<-paste0(str.factors,paste0("u.",idq.facs[i],","));
        str.onCond <-paste0(str.onCond, paste0("t.",idq.facs[i],"=u.",idq.facs[i]," and \n"));
      }
    }
  }
  cat("#--str.factors string = ",str.factors,"\n")
  qry<-gsub("&&factors",  str.factors, qry, fixed=TRUE);
  qry<-gsub("&&idq.size",  idq.size,     qry, fixed=TRUE);
  qry<-gsub("&&idq.value", idq.value,    qry, fixed=TRUE);
  qry<-gsub("&&onCond",   str.onCond,  qry, fixed=TRUE);
  if (verbose) cat("Query to expand abundance to sizes:\n",qry,"\n");
  tmp3<-sqldf::sqldf(qry);
  tmp3[[id.value]][is.na(tmp3[[id.value]])]<-0;

  return(tmp3);
  #
  # # #--Step 4: expand sample sizes to all sizes
  # cat("\n#--Step 4-------------------------------------\n");
  # #--Example query:
  # # qry<-"select
  # #         u.fishery,u.area,u.`shell condition`,u.year,
  # #         u.size,t.ss
  # #       from uFZs as u left join tmp as t
  # #       on
  # #         t.fishery = u.fishery and
  # #         t.area    = u.area and
  # #         t.`shell condition` = u.`shell condition` and
  # #         t.year    = u.year
  # #       order by u.fishery,u.area,u.`shell condition`,u.year,u.size;";
  # qry<-"select
  #         &&factors
  #         u.&&id.size,t.&&id.ss
  #       from uFZs as u left join tmp as t
  #       on
  #         &&onCond
  #         t.&&id.size    = u.&&id.size
  #       order by &&factors u.&&id.size;";
  # qry<-gsub("&&factors", str.factors, qry, fixed=TRUE);
  # qry<-gsub("&&id.size", id.size,     qry, fixed=TRUE);
  # qry<-gsub("&&id.ss",   id.ss,       qry, fixed=TRUE);
  # qry<-gsub("&&onCond",  str.onCond,  qry, fixed=TRUE);
  # if (verbose) cat("Query to expand ss to sizes:\n",qry,"\n");
  # tmp4<-sqldf::sqldf(qry);
  # tmp4[[id.ss]][is.na(tmp4[[id.ss]])]<-0;
  #
  # #--Step 5: combine abundance and sample size information
  # cat("\n#--Step 5-------------------------------------\n");
  # #--Example query:
  # # qry<-"select
  # #         u.fishery,u.area,u.`shell condition`,u.year,
  # #         u.ss,u.size,t.abundance
  # #       from tmp4 as u, tmp3 as t
  # #       on
  # #         t.fishery = u.fishery and
  # #         t.area    = u.area and
  # #         t.`shell condition` = u.`shell condition` and
  # #         t.year    = u.year and
  # #         t.size    = u.size
  # #       order by u.fishery,u.area,u.`shell condition`,u.year,u.size;";
  # qry<-"select
  #         &&factors
  #         u.&&id.ss,u.&&id.size,t.&&id.value
  #       from tmp4 as u, tmp4 as t
  #       on
  #         &&onCond
  #         t.size = u.size
  #       order by &&factors u.&&id.size;";
  # qry<-gsub("&&factors",  str.factors, qry, fixed=TRUE);
  # qry<-gsub("&&id.ss",    id.ss,       qry, fixed=TRUE);
  # qry<-gsub("&&id.size",  id.size,     qry, fixed=TRUE);
  # qry<-gsub("&&id.value", id.value,    qry, fixed=TRUE);
  # qry<-gsub("&&onCond",   str.onCond,  qry, fixed=TRUE);
  # if (verbose) cat("Query to join ss and size composition tables:\n",qry,"\n");
  # tmp5<-sqldf::sqldf(qry);
  #
  # return(tmp5);
}


#'
#' @title Normalize "size" compositions across a set of factors
#'
#' @description Function to create a dataframe of "size" compositions normalized (summing to 1) across a set of factors.
#'
#' @param dfrZ - input dataframe with count/abundance data by size and other factors
#' @param id.size - name of column in dfr with "size" data
#' @param id.value - name of column in dfr with count (or abundance) values
#' @param id.vars - character vector of column names used as factors in calculating the size compositions
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
  tmp<-dfr;
  #--normalize size compositions by norm.vars so
  #--they sum to 1 across remaining id.vars.
  #back-quote column names
  #----Step 1: backquote column names
  if (verbose) cat("#--Step 1-------------------------------------\n");
  idq.size<-paste0("`",id.size,"`");
  idq.value<-paste0("`",id.value,"`");
  if (wtsUtilities::isBlankString(id.facs)) id.facs<-names(tmp)[!(names(tmp) %in% c(id.size,id.value))];
  idq.facs<-paste0("`",id.facs,"`");
  if (wtsUtilities::isBlankString(id.normfacs)) idq.normfacs<-paste0("`",id.normfacs,"`");
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
  if (verbose){cat("#--str.norm.vars: '",str.norm.vars,"'\n",sep="");}
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

  return(tmp2);
}

#'
#' @title Scale "size" compositions in a dataframe by factor level-specific scalars
#'
#' @description Function to scale "size" compositions in a dataframe by factor level-specific scalars.
#'
#' @param dfrZCs - input dataframe with size composition data by size and other factors
#' @param dfrScalars - input dataframe with scalars for a subset of the factors in dfrZCs
#' @param id.value - name of column in dfrZCs with (possibly normalized) count/abundance values
#' @param id.scalevars - character vector of column names of factors for the values in dfrScalars to scale the size compositions in dfrZCs by
#' @param id.scalevalue - name of column in dfrScalars with the values used to scale the size compositions in dfrZCs
#' @param verbose - flag to print diagnostic output to console
#'
#' @return dataframe with scaled compositions by size (and possibly other factors).
#'
#' @details The names of columns in dfrScalars used as factors to scale the
#' composition data in dfrZCs should be the same in both datasets.
#'
#' @export
#'
scaleSizeComps<-function(dfrZCs,
                         dfrScalars,
                         id.value="fraction",
                         id.scalefacs="",
                         id.scalevalue="abundance",
                         verbose=FALSE){
  tmp<-dfrZCs;
  #----Step 1: backquote column names
  if (verbose) cat("#--Step 1-------------------------------------\n");
  idq.value<-paste0("`",id.value,"`");
  id.facs<-paste0("`",names(tmp)[!(names(tmp)==id.value)],"`");
  idq.facs<-paste0("`",id.facs,"`");
  if (!wtsUtilities::isBlankString(id.scalefacs)) idq.scalefacs<-paste0("`",id.scalefacs,"`");
  idq.scalevalue<-paste0("`",id.scalevalue,"`");

  #--Step 2: scale "normalized" size compositions by scalar values
  if (verbose) cat("#--Step 2-------------------------------------\n");
  #--Example query:
  # qry<-"select
  #         t.fishery, t.area, t.`shell condition`,
  #         t.year, t.ss, t.size,
  #         t.fraction*s.abundance as abundance
  #       from tmp as t, dfrScalars as s
  #       where
  #         t.year    = s.year and
  #         t.fishery = s.fishery and
  #         t.area    = s.area;";
  qry<-"select
          &&idq.facs,
          t.&&idq.value*s.&&idq.scalevalue as &&idq.scalevalue
        from tmp as t, dfrScalars as s
        &&whereCond;";
  str.idq.facs<-paste("where t",idq.facs,sep=".",collapse=",");
  str.where<-"";
  if (!wtsUtilities::isBlankString(id.scalefacs)){
    str.where<-paste0("where t.",idq.scalefacs[1],"=s.",idq.scalefacs[1]);
    if (length(id.scalefacs)>1){
      for (i in 2:length(id.scalefacs)) str.where<-paste0(str.where," and \n",
                                                          paste0("t.",idq.scalefacs[i],"=s.",idq.scalefacs[i]));
    }
  }
  qry<-gsub("&&idq.facs",       str.idq.facs,   qry, fixed=TRUE);
  qry<-gsub("&&idq.value",      idq.value,      qry, fixed=TRUE);
  qry<-gsub("&&idq.scalevalue", idq.scalevalue, qry, fixed=TRUE);
  qry<-gsub("&&whereCond",      str.where,      qry, fixed=TRUE);
  if (verbose){cat("#--idq.facs    : '",idq.facs,    "'\n",sep="");}
  if (verbose){cat("#--str.idq.facs: '",str.idq.facs,"'\n",sep="");}
  if (verbose){cat("#--str.where   : '",str.where,   "'\n",sep="");}
  if (verbose) cat("Query to calculate normalized compositions:\n",qry,"\n");
  tmp<-sqldf::sqldf(qry);

  return(tmp);
}

