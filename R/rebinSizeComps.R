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
                         id.facs="",
                         cutpts=NULL,
                         truncate.low=TRUE,
                         truncate.high=FALSE,
                         expandToAllFactorCombos=TRUE,
                         verbose=FALSE){
  if (verbose) cat("\n\n#----Starting rebinSizeComps\n");
  tmp<-dfr;
  #--Step 0: back-quote column names
  if (verbose) cat("#--Step 0-------------------------------------\n");
  idq.size<-paste0("`",id.size,"`");
  idq.value<-paste0("`",id.value,"`");
  if (wtsUtilities::isBlankString(id.facs)) id.facs<-names(tmp)[!(names(tmp) %in% c(id.size,id.value))];
  idq.facs<-paste0("`",id.facs,"`");

  #--Step 1: rebin compositions (if necessary) to new size bins
  if (verbose) cat("#--Step 1-------------------------------------\n");
  if (verbose) cat("#----nrow(dfr) =",nrow(tmp),"\n");
  if (is.null(cutpts)) {
    cutpts<-data.frame(size=min(tmp$size,na.RM=TRUE):(max(tmp$size,na.rm=TRUE)+1));
    newcutpts<-cutpts;
  } else {
    nCPs<-length(cutpts);
    newcutpts<-cutpts;#make copy to apply truncation correctly
    if (!truncate.low ) newcutpts[1]<-0;
    if (!truncate.high) newcutpts[nCPs]<-Inf;
  }
  if (verbose) {
    cat("#--no. cutpoints :",length(cutpts),"\n");
    cat("#--cutpoints:",cutpts,"\n");
  }
  cuts<-cut(tmp[[id.size]],newcutpts,right=FALSE,labels=FALSE);#make cuts based on new bins adjusted for truncation
  tmp[[id.size]]<-cutpts[cuts];                                #assign to original bins using cutpoints
  tmp<-tmp[!is.na(tmp[[id.size]]),];#drop truncated data
  #----Example query:
  # qry<-"select
  #         fishery,area,`shell condition`,year,size,
  #         sum(abundance) as abundance
  #       from tmp
  #       group by fishery,area,`shell condition`,year,size
  #       order by fishery,area,`shell condition`,year,size;";
  qry<-"select
          &&factors &&idq.size,
          sum(&&idq.value) as &&idq.value
        from tmp
        group by &&factors &&idq.size
        order by &&factors &&idq.size;";
  #determine factor string for query
  str.factors<-"";
  if (!wtsUtilities::isBlankString(id.facs)) {
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
  if (verbose) cat("nrow(rebinned) =",nrow(tmp1),"\n");

  #--Step 2: determine unique factor x size combinations
  if (verbose) cat("\n#--Step 2-------------------------------------\n");
  uFZs<-getUniqueFZs(tmp,
                     id.size=id.size,
                     id.value=id.value,
                     id.facs=id.facs,
                     cutpts=cutpts,
                     expandToAllFactorCombos=expandToAllFactorCombos,
                     verbose=verbose);

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
          u.&&idq.size,t.&&idq.value
        from uFZs as u left join tmp1 as t
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
  if (verbose) cat("#--str.factors string = ",str.factors,"\n")
  qry<-gsub("&&factors",  str.factors, qry, fixed=TRUE);
  qry<-gsub("&&idq.size",  idq.size,     qry, fixed=TRUE);
  qry<-gsub("&&idq.value", idq.value,    qry, fixed=TRUE);
  qry<-gsub("&&onCond",   str.onCond,  qry, fixed=TRUE);
  if (verbose) cat("Query to expand abundance to sizes:\n",qry,"\n");
  tmp3<-sqldf::sqldf(qry);
  idx<-is.na(tmp3[[id.value]]);
  tmp3[idx,id.value]<-0;
  if (verbose) {
    cat("#--setting",sum(idx),"NA values to 0\n");
    cat("#--nrow(final) =",nrow(tmp3),"\n");
  }

  if (verbose) cat("#----Finished rebinSizeComps()\n");
  return(tmp3);
}
