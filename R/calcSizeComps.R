#'
#' @title Create "size" compositions from a dataframe with count data
#'
#' @description Function to create a dataframe of "size" compositions from a dataframe with count data.
#'
#' @param dfrZ - input dataframe with count/abundance data by size and other factors
#' @param id.size - name of column in dfr with "size" data
#' @param id.value - name of column in dfrZ with count (or abundance) values
#' @param id.facs - character vector of column names to use as factors in calculating size compositions
#' @param cutpts - cut points for size bins
#' @param truncate.low -flag to truncate counts/abundance below the minimum size or include in first size bin (FALSE)
#' @param truncate.high - flag to truncate counts/abundance above the maximum size (TRUE) or include in final size bin (FALSE)
#' @param expandToAllFactorCombos - flag (T/F) to expand size comps to include all factor combinations
#' @param normalize - flag (T/F) to normalize size comps (i.e., sum to 1 across normalizing factors)
#' @param id.normfacs - character vector of column names of factors to normalize the size compositions across when scaling to total abundance
#' @param rescale - flag to rescale size compositions by factor combination using scalars in dfrScalars
#' @param dfrScalars - [optional] input dataframe with scalars to apply by a subset of the factors in dfrZ
#' @param id.scalefacs - character vector of column names of factors in dfrScalars to scale the size compositions across
#' @param id.scalevalue - column name for scale value in dfrScalars
#' @param verbose - flag to print diagnostic output to console
#'
#' @return dataframe with summed counts (or abundance) as compositions by size (and possibly other factors).
#'
#' @details Non-standard column names should be back-quoted. The names of columns (id.normfacs) used to normalize the composition data in dfrZ
#' prior to expanding to total abundance using dfrScalars should represent the same factors (id.scalefacs) and levels in both datasets.
#'
#' @importFrom wtsUtilities isBlankString
#'
#' @export
#'
calcSizeComps<-function(dfrZ,
                        id.size="size",
                        id.value="abundance",
                        id.facs="",
                        cutpts=NULL,
                        truncate.low=FALSE,
                        truncate.high=FALSE,
                        expandToAllFactorCombos=FALSE,
                        normalize=FALSE,
                        id.normfacs="",
                        rescale=FALSE,
                        dfrScalars=NULL,
                        id.scalevalue,
                        id.scalefacs=id.normfacs,
                        verbose=FALSE){
  tmp0<-dfrZ;
  if (verbose) cat("#-In calcSizeComps(): names(tmp0) = ",names(tmp0),"\n");
  #--Step 1: bin counts/abundance to size bins and expand to all sizes
  tmp1<-rebinSizeComps(tmp0,
                       id.size=id.size,
                       id.value=id.value,
                       id.facs=id.facs,
                       cutpts=cutpts,
                       truncate.low=truncate.low,
                       truncate.high=truncate.high,
                       expandToAllFactorCombos=expandToAllFactorCombos,
                       verbose=verbose);
  if (verbose) cat("#-In calcSizeComps(): number of NAs after rebinning =",sum(is.na(tmp1[[id.value]])),"\n")
  if (verbose) cat("#-In calcSizeComps(): names(tmp1) = ",names(tmp1),"\n");

  #--Step 2: add sample sizes in column "ss"
  tmp2 <- addSampleSizes(tmp1,
                         id.size=id.size,
                         id.value=id.value,
                         verbose=verbose);
  if (verbose) cat("#-In calcSizeComps(): number of NAs after adding sample sizes =",sum(is.na(tmp2[[id.value]])),"\n")
  if (verbose) cat("#-In calcSizeComps(): names(tmp2) = ",names(tmp2),"\n");

  if (normalize){
    #--Step 3: normalize size compositions to sum to 1 across a subset of factors
    idp.facs<-"ss"
    if (!wtsUtilities::isBlankString(id.facs)){
      idp.facs<-c(id.facs,"ss");
    }
    tmp3<-normalizeSizeComps(tmp2,
                             id.size=id.size,
                             id.value=id.value,
                             id.facs=idp.facs,
                             id.normfacs=id.normfacs,
                             verbose=verbose);
  } else {
    tmp3 <- tmp2;
  }
  idp.value <- id.value;
  if (normalize) idp.value <- "fraction";
  if (verbose) cat("#-In calcSizeComps(): number of NAs after normalizing =",sum(is.na(tmp3[[idp.value]])),"\n")
  if (verbose) cat("#-In calcSizeComps(): names(tmp3) = ",names(tmp3),"\n");

  if (rescale){
    #--step 4: re-scale size compositions across a subset of factors
    tmp4<-scaleSizeComps(tmp3,
                         dfrScalars,
                         id.value=idp.value,
                         id.scalefacs=id.scalefacs,
                         id.scalevalue=id.scalevalue,
                         verbose=verbose);
  } else {
    tmp4 <- tmp3;
  }
  idp.value <- id.value;
  if (rescale) idp.value <- id.scalevalue;
  if (verbose) cat("#-In calcSizeComps(): number of NAs after rescaling =",sum(is.na(tmp4[[idp.value]])),"\n")
  if (verbose) cat("#-In calcSizeComps(): names(tmp4) = ",names(tmp4),"\n");

  return(tmp4);
}
