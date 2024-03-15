#'
#' Plot annual size comps (or other similar) using \pkg{ggridges}
#'
#' @description Function to plot annual size comps (or other similar) using \pkg{ggridges}
#'
#' @param dfr_ : input dataframe
#' @param sizes_in : name of column with "size" information (default = "size")
#' @param values_in : name of column with "abundance" information (default = "abundance")
#' @param y_positions : name of column with vertical axis location information (default = "year")
#' @param removeZeros : remove comps with zero abundance in all bins(default = TRUE)
#' @param reverseY : reverse y scale (default = FALSE)
#' @param colour : column name colour variable (default = NULL)
#' @param fill : column name for fill variable (default = NULL)
#' @param group_by : string vector with column names for grouping before summarizing (default = value of y_positions)
#' @param normalize_by : string vector with column names for grouping before normalizing (default = NULL)
#' @param xlim : x axis limits (default = c(0,200))
#' @param ylim : y axis limits (default = c(1980,2022))
#' @param x_breaks : x axis breaks (default = seq(5,200,5))
#' @param y_breaks : y axis breaks (default = seq(1900,2100,5))
#' @param legend.position : legend position (as in [ggplot2::theme()], default = c(0.98,0.98))
#' @param legend.justification : legend justification on plot (as in [ggplot2::theme()], default = c(1,1))
#'
#' @return \pkg{ggplot2} plot object
#'
#' @details
#' Uses [ggridges::geom_density_ridges()] to plot the size compositions.
#' The returned plot object can be faceted, etc., for more complex plotting situations.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges
#' @importFrom rlang sym syms
#' @importFrom scales squish
#' @importFrom wtsPlots getStdTheme
#' @importFrom wtsUtilities Sum
#'
#' @export
#'
plotSizeCompsAsRidges<-function(dfr_,
                          sizes_in="size",
                          values_in="abundance",
                          y_positions="year",
                          removeZeros=TRUE,
                          reverseY=FALSE,
                          colour=NULL,
                          fill=NULL,
                          group_by=y_positions,
                          normalize_by=NULL,
                          xlim=c(100,200),
                          ylim=c(1980,2022),
                          x_breaks=seq(5,200,5),
                          y_breaks=seq(1900,2100,5),
                          legend.position=c(0.98,0.98),
                          legend.justification = c(1,1)){
  #--keep only relevant columns
  keep_cols = names(dfr_) %in% c(group_by,normalize_by,y_positions,sizes_in,values_in);
  dfr = dfr_[,keep_cols];
  #--create syms for size, value and y_position columns
  sym_siz = rlang::sym(sizes_in);
  sym_val = rlang::sym(values_in);
  sym_yps = rlang::sym(y_positions);
  sym_clr = rlang::sym(colour);
  sym_fil = rlang::sym(fill);
  #--aggregate size comps by size across grouping variables (if requested)
  if (!is.null(group_by)){
    syms_vrs = rlang::syms(c(group_by,sizes_in));
    dfr = dfr |>
              dplyr::group_by(!!!syms_vrs) |>
              dplyr::summarize(tot_=wtsUtilities::Sum(!!sym_val)) |>
              dplyr::ungroup();
    rm(syms_vrs);
  } else {
    dfr[["tot_"]] = dfr[[values_in]];#--copy values to column named "tot_"
  }
  #--remove size comps with zero "mass" at all sizes (if requested)
  if (removeZeros){
    #--calculate years with non-zero size comps
    syms_vrs = rlang::syms(c(group_by));
    tmp = dfr |> dplyr::group_by(!!!syms_vrs) |>
                    dplyr::summarize(tot_=wtsUtilities::Sum(tot_)) |>
                    dplyr::ungroup() |>
                    dplyr::filter(tot_>0) |>
                    dplyr::select(!!!syms_vrs);
    #--keep only those
    dfr = dfr |> dplyr::inner_join(tmp);
    rm(syms_vrs,tmp);
  }
  #--normalize size comps over normalizing variables (if requested)
  if (!is.null(normalize_by)){
    syms_vrs = rlang::data_syms(normalize_by);
    dfr = dfr |>
              dplyr::group_by(!!!syms_vrs) |>
              dplyr::mutate(tot_=tot_/wtsUtilities::Sum(tot_)) |>
              dplyr::ungroup();
    rm(syms_vrs);
  }
  #------plot size comps as ridges
  if (!is.null(group_by)){
    syms_vrs = rlang::syms(group_by);
  } else {
    rlang::syms(y_positions);
  }
  p = ggplot(data=dfr,mapping=aes(x=!!sym_siz,y=!!sym_yps,colour=!!sym_clr,fill=!!sym_fil)) +
       ggridges::geom_density_ridges(mapping=aes(group=paste0(!!!syms_vrs),height=tot_),stat="identity",alpha=0.4)+
       guides(fill=guide_legend(override.aes=list(alpha=1)))+
       # geom_segment(data=tmp4,mapping=aes(x=medianSize,y=y0,xend=medianSize,yend=y1,colour=area)) +
       # geom_point(data=tmp4,mapping=aes(x=medianSize,y=y0,colour=area,fill=area)) +
       scale_x_continuous(breaks=x_breaks,limits=xlim) +
       scale_y_continuous(breaks=y_breaks,limits=ylim,oob=scales::squish) +
       labs(x="size (mm CW)",y="crab year") +
       wtsPlots::getStdTheme() +
       theme(legend.position=legend.position,
             legend.justification=legend.justification);
  if (reverseY) p = p + scale_y_reverse();
  return(p)
}

#--function testing
# p = plotZCsAsRidges(dfrRC_ZCs_ByFAYXMS,
#                 sizes_in="size",
#                 values_in="abundance",
#                 y_positions="year",
#                 removeZeros=TRUE,
#                 reverseY=TRUE,
#                 colour="area",
#                 fill="area",
#                 group_by=c("year","area"),
#                 normalize_by=c("year","area"),
#                 xlim=c(100,200),
#                 ylim=c(1980,2022),
#                 x_breaks=seq(5,200,5),
#                 y_breaks=seq(1900,2100,5),
#                 legend.position=c(0.98,0.98),
#                 legend.justification = c(1,1));
# print(p);

#--for "internal" testing
# dfr_ = dfrRC_ZCs_ByFAYXMS |> dplyr::filter(fishery=="TCF",year>=1980,size>100);
# sizes_in="size";
# values_in="abundance";
# y_positions="year";
# removeZeros=TRUE;
# colour="area";
# fill="area";
# group_by=c("year","area");
# normalize_by=c("year","area");
# x_breaks=seq(5,200,5);
# y_breaks=seq(1900,2100,5);
# legend.position=c(0.98,0.98);
# legend.justification=c(1,1);
# #--clean up "internal" assignments
# rm(dfr_,sizes_in,values_in,y_positions,removeZeros,colour,fill,
#    group_by,normalize_by,x_breaks,y_breaks,legend.position,legend.justification);



