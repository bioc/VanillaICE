#' Constructor for HmmTrellisParam class
#'
#' @param ylimits length-two list of the y-axis limits for B allele
#' frequencies and log R ratios, respectively
#' @param expandfun a function that takes a length-one \code{GRanges}
#' object as an argument and computes a width relative to the width of
#' the \code{GRanges} object
#' @export
HmmTrellisParam <- function(ylimits=list(c(0,1), c(-3,1)),
                            expandfun=function(g){  width(g) * 50}){
  new("HmmTrellisParam", ylimits=ylimits, expandfun=expandfun)
}

expandFun <- function(object) object@expandfun
yLimits <- function(object) object@ylimits

#' Default viewports for plotting CNV data with lattice-style graphics
#'
#' @examples
#' vps <- viewports()
#' @seealso \code{\link{xyplotList}} \code{\link{xygrid}}
#' @return list
#' @export
viewports <- function(){
  data_vp <- viewport(x=0, y=0.01, width=unit(0.99, "npc"), height=unit(0.9, "npc"),
                      name="vp1", just=c("left", "bottom"))
  idiogram_vp <- viewport(x=unit(0.23,"npc"), y=unit(0.85, "npc"), width=unit(0.72, "npc"),
                          height=unit(0.15, "npc"), name="vp1", just=c("left", "bottom"))
  title_left <- viewport(x=0.05, y=0.85, width=unit(0.2, "npc"),
                         height=unit(0.15, "npc"), name="vp2",just=c("left", "bottom"))
##  title_right <- viewport(x=unit(0.8, "npc"), y=unit(0.93, "npc"), width=unit(0.2, "npc"),
##                          height=unit(0.05, "npc"), name="vp2",just=c("left", "bottom"))
  list(data=data_vp, idiogram=idiogram_vp, title_left=title_left)
}

.find_xlim_percent <- function(g, percent=0.05){
  wd <- width(g)
  w <- wd/percent
  d <- (w-wd)*1/2
  st <- max(start(g)[1]-d, 1)
  en <- min(end(g)[1]+d, seqlengths(g))
  lim <- as.integer(c(st, en))
  ##ILimit(start=lim[1], end=lim[2])
  lim
}


#' @param trellis_plot  an object of class \code{trellis}
#' @param viewports a list of viewports as provided by the \code{viewports} function
#' @rdname plotting
#' @seealso \code{\link{viewports}}
#' @export
xygrid <- function(trellis_plot, viewports, granges){
  vp1 <- viewports[["data"]]
  vp2 <- viewports[["title_left"]]
##  vp3 <- viewports[["title_right"]]
  vp4 <- viewports[["idiogram"]]
  grid.newpage()
  pushViewport(vp1)
  print(trellis_plot, vp=1, newpage=FALSE)
  upViewport()
  pushViewport(vp2)
  locs <- prettyNum(round(c(start(granges), end(granges))/1000, 1), big.mark=",")
  label <- paste0(seqnames(granges), ": ", locs[1], "-", locs[2], "kb")
  label <- paste0(granges$id, "\n", label)
  state <- granges$state
  state <- as.integer(factor(state, levels=c(1,2,5,6)))
  CN <- c("0", "1", "3", "4+")[state]
  if("prCall" %in% colnames(mcols(granges))){
    prob <- round(granges$prCall, 3)
    label <- paste0(label, "\nCN=", CN, "\nPr(CN|data)=", prob, "\n")
  }
  grid.text(label, y=unit(0, "npc"), x=unit(0, "npc"),
            just=c("left", "bottom"), gp=gpar(cex=0.7))
  seqlevels(granges, force=TRUE) <- chromosome(granges)
  xlim <- .find_xlim_percent(granges, 0.05)
  iparams <- IdiogramParams(seqnames=chromosome(granges),
                            genome=genome(granges)[[1]],
                            seqlengths=seqlengths(granges),
                            box=list(xlim=xlim, color="orange"))
  idiogram <- plotIdiogram(iparams)
  upViewport(0)
  pushViewport(vp4)
  print(idiogram, vp=vp4, newpage=FALSE)
}
