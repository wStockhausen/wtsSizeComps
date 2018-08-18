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
  if (verbose) cat("\n\n#----Starting scaleSizeComps\n");
  tmp<-dfrZCs;
  #----Step 1: backquote column names
  if (verbose) cat("#--Step 1-------------------------------------\n");
  idq.value<-paste0("`",id.value,"`");
  if (verbose) cat("#--idq.value: ",idq.value,"\n");
  if (verbose) cat("#--names(tmp): ",names(tmp),"\n");
  id.facs<-names(tmp)[!(names(tmp)==id.value)];
  idq.facs<-paste0("`",id.facs,"`");
  if (verbose) cat("#--idq.facs: ",idq.facs,"\n");
  if (!wtsUtilities::isBlankString(id.scalefacs)) {
    idq.scalefacs<-paste0("`",id.scalefacs,"`");
    if (verbose) cat("#--idq.scalefacs: ",idq.scalefacs,"\n");
  }
  idq.scalevalue<-paste0("`",id.scalevalue,"`");
  if (verbose) cat("#--idq.scalevalue: ",idq.scalevalue,"\n");

  #--Step 2: scale size compositions by scalar values
  if (verbose) cat("#--Step 2-------------------------------------\n");
  #--Example query:
  # qry<-"select
  #         t.fishery, t.area, t.`shell condition`,
  #         t.year, t.ss, t.size,
  #         t.fraction*s.abundance as abundance
  #       from tmp as t left join dfrScalars as s
  #       on
  #         t.year    = s.year and
  #         t.fishery = s.fishery and
  #         t.area    = s.area;";
  qry<-"select
          &&idq.facs,
          t.&&idq.value*s.&&idq.scalevalue as &&idq.scalevalue
        from tmp as t left join dfrScalars as s
        &&onCond;";
  str.idq.facs<-paste("t",idq.facs,sep=".",collapse=",");
  str.on<-"";
  if (!wtsUtilities::isBlankString(id.scalefacs)){
    str.on<-paste0("on t.",idq.scalefacs[1],"=s.",idq.scalefacs[1]);
    if (length(id.scalefacs)>1){
      for (i in 2:length(id.scalefacs)) str.on<-paste0(str.on," and \n",
                                                       paste0("t.",idq.scalefacs[i],"=s.",idq.scalefacs[i]));
    }
  }
  qry<-gsub("&&idq.facs",       str.idq.facs,   qry, fixed=TRUE);
  qry<-gsub("&&idq.value",      idq.value,      qry, fixed=TRUE);
  qry<-gsub("&&idq.scalevalue", idq.scalevalue, qry, fixed=TRUE);
  qry<-gsub("&&onCond",      str.on,      qry, fixed=TRUE);
  if (verbose){cat("#--idq.facs    : '",idq.facs,    "'\n",sep="");}
  if (verbose){cat("#--str.idq.facs: '",str.idq.facs,"'\n",sep="");}
  if (verbose){cat("#--str.on      : '",str.on,      "'\n",sep="");}
  if (verbose) cat("Query to calculate scaled compositions:\n",qry,"\n");
  tmp<-sqldf::sqldf(qry);
  idx<-is.na(tmp[[id.scalevalue]]);
  tmp[idx,id.scalevalue]<-0;
  if (verbose) {
    cat("#--setting",sum(idx),"NA values to 0\n");
    cat("#--nrow(final) =",nrow(tmp),"\n");
  }

  if (verbose) cat("#----Finished scaleSizeComps()\n");
  return(tmp);
}

