#'
#' @title Determine the unique factor and size combinations from a dataframe of "size" compositions
#'
#' @description Function to determine the unique factor and size combinations from a dataframe of "size" compositions.
#'
#' @param dfr - input dataframe with composition data by "size" and other factors
#' @param id.size - name of column in dfr with "size" data
#' @param id.value - name of column in dfr with count (or abundance) values
#' @param id.facs - character vector of factor column names
#' @param cutpts - vector of cutpoints for size compositions (or NULL to extract from dfr)
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
                       cutpts=NULL,
                       expandToAllFactorCombos=FALSE,
                       verbose=verbose){
  if (verbose) cat("\n\n#----Starting getUniqueFZs\n");
  tmp<-dfr;
  #--Step 1: back-quote column names
  if (verbose) cat("\n#--Step 1-------------------------------------\n");
  idq.size<-paste0("`",id.size,"`");
  idq.value<-paste0("`",id.value,"`");
  idq.facs<-paste0("`",id.facs,"`");

  #--Step 2: determine unique factor combinations
  if (verbose) cat("\n#--Step 2-------------------------------------\n");
  if (wtsUtilities::isBlankString(id.facs)) id.facs<-names(tmp)[!(names(tmp) %in% c(id.size,id.value))];
  if (verbose) cat("#--id.facs =",paste(paste0("'",id.facs,"'",collapse=","),"\n"));
  tmp2<-unique(tmp[,id.facs]);
  uFs<-tmp2;

  #--Step 3: determine unique sizes
  if (verbose) cat("\n#--Step 3-------------------------------------\n");
  if (is.null(cutpts)){
    qry<-"select distinct &&idq.size from tmp order by &&idq.size";
    qry<-gsub("&&idq.size",idq.size,qry, fixed=TRUE);
    uZs<-sqldf::sqldf(qry);
  } else {
    uZs <- list();
    uZs[[id.size]]<-cutpts[1:(length(cutpts)-1)];
    uZs <- as.data.frame(uZs,optional=TRUE);
  }
  if (verbose) cat("#--unique Zs:",uZs[[id.size]]);

  #--Step 4: determine unique factor x size combinations
  if (verbose) cat("\n#--Step 4-------------------------------------\n");
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

  if (verbose) cat("#----Finished getUniqueFZs\n");
  return(uFZs);
}