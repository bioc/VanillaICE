#' Constructor for IdiogramParam objects
#'
#' Parameters for plotting idiograms
#' @param seqnames length-one character vector providing chromosome name
#' @param seqlengths length-one numeric vector indicating size of chromosome
#' @param unit character string indicating unit for genomic position
#' @param genome character string indicating genome build
#' @param box a list of parameters for plotting the box around the part of the idiogram that is plotted
#' @return IdiogramParam object
#' @rdname IdiogramParam-class
#' @export
IdiogramParams <- function(seqnames=character(),
                           seqlengths=numeric(), unit="kb", genome="hg19",
                           box=list(xlim=ILimit(), color="blue", lwd=1)){
  new("IdiogramParams", seqnames=seqnames, seqlengths=seqlengths, unit=unit,
      genome=genome, box=box)
}

setGeneric("boxIdiogram", function(object) standardGeneric("boxIdiogram"))
setMethod("chromosome", "IdiogramParams", function(object) object@seqnames)
setMethod("seqnames", "IdiogramParams", function(x) x@seqnames)
setMethod("seqlengths", "IdiogramParams", function(x) x@seqlengths)
setMethod("genome", "IdiogramParams", function(x) x@genome)
setMethod("boxIdiogram", "IdiogramParams", function(object) object@box)

#' @param object an IdiogramParam object
#' @aliases show,IdiogramParams-method
#' @rdname IdiogramParams-class
setMethod("show", "IdiogramParams", function(object){
  cat("'IdiogramParams' class \n")
  cat("   seqnames  : ", seqnames(object), "\n")
  cat("   genome    : ", genome(object), "\n")
  cat("   seqlengths: ", seqlengths(object), "\n")
})


panelIdiogram <- function(params){
  panel.xyplot(0, 0, type="n", col="transparent")
  build <- genome(params)
  cytoband <- getCytoband(params)
  cytoband_p <- cytoband[grep("^p", rownames(cytoband), value=TRUE), ]
  cytoband_q <- cytoband[grep("^q", rownames(cytoband), value=TRUE), ]
  N <- nrow(cytoband)
  p.bands <- nrow(cytoband_p)

  left_cuts <- c(1, p.bands+1)
  right_cuts <- c(N, p.bands)
  cut.right <- cut.left <- rep(FALSE, N)
  cut.left[left_cuts] <- TRUE
  cut.right[right_cuts] <- TRUE
  ## this is a "stalk", do not draw box. Draw two vertical lines instead
  is_stalk <- cytoband[, "gieStain"] == "stalk"
  index_stalk <- which(is_stalk)
  cut.right[index_stalk - 1] <- TRUE
  cut.left[index_stalk  + 1] <- TRUE

  colors <- .cytobandColors(cytoband[, "gieStain"])
  xx <- c(0, cytoband[nrow(cytoband), "end"])
  yy <- c(0,1)

  starts <- cytoband[, "start"]
  ends <- cytoband[, "end"]
  if(any(is_stalk)) .drawStalk(starts[is_stalk], ends[is_stalk])
  taper_right <- taper_left <- rep(0, length(starts))
  taper_left[cut.left] <- 0.15
  taper_right[cut.right] <- 0.15
  .drawPolygon(starts, ends, colors, taper_left=taper_left, taper_right=taper_right)

  box_param <- boxIdiogram(params)
  xlim <- box_param$xlim
  col <- box_param$color
  lwd <- box_param$lwd
  panel.rect(xleft=xlim[1], xright=xlim[2], ybottom=0, ytop=1,
             border=col,
             lwd=lwd)
}



getCytoband <- function(params){
  path <- system.file("extdata", package="SNPchip", mustWork=TRUE)
  cytoband <- read.table(file.path(path, paste("cytoBand_", genome(params), ".txt", sep="")), as.is=TRUE, header=FALSE)
  colnames(cytoband) <- c("chrom", "start", "end", "name", "gieStain")
  cytoband <- cytoband[cytoband[, "chrom"] %in% seqnames(params), ]
  rownames(cytoband) <- as.character(cytoband[, "name"])
  return(cytoband)
}

plotIdiogram <- function(object){
  dat <- data.frame(x=c(1, seqlengths(object)),
                    y=c(0, 1))
  fig <- xyplot(y~x, dat, xlim=c(range(dat$x)),
                ylim=c(0,1), type="n", xlab="", ylab="",
                scales=list(draw=FALSE),
                par.settings=list(axis.line=list(col="transparent")),
                panel=panelIdiogram, params=object)
}

#' @param x an IdiogramParam object
#' @param y ignored
#' @param ... ignored
#' @export
#' @aliases plot,IdiogramParams-method
#' @rdname IdiogramParam-class
setMethod("plot", "IdiogramParams", function(x, y, ...)  plotIdiogram(x))


.cytobandColors <- function(stains){
  colors <- rep("white", length(stains))
  colors[stains == "gneg" | stains == "gvar"] <- "gray100"
  colors[stains=="gpos25"] <- "gray90"
  colors[stains=="gpos50"] <- "gray70"
  colors[stains=="gpos75"] <- "gray40"
  colors[stains=="gpos100"] <- "gray0"
  colors[stains=="stalk"] <- "brown3"
  colors[stains=="acen"] <- "brown4"
  colors
}

.drawStalk <- function(starts, ends){
  deltas <- (ends-starts)/3
  y0 <- rep(0.2, length(starts))
  y1 <- rep(0.8, length(starts))
  lsegments(starts+deltas, y0, starts+deltas, y1)
  lsegments(ends+deltas, y0, ends+deltas, y1)
}

.drawPolygon <- function(starts, ends, colors, taper_left=0, taper_right=0){
  top <- 0.8
  bot <- 0.2
  h <- top-bot
  deltas <- (ends-starts)/4
  ## how to vectorize?
  for(i in seq_along(starts)){
    yy <- c(bot + taper_left[i]*h, bot, bot, bot + taper_right[i]*h, top - taper_right[i]*h,
            top, top, top - taper_left[i]*h)
    lpolygon(x=c(starts[i], starts[i]+deltas[i],
               ends[i]-deltas[i],
               rep(ends[i], 2),
               ends[i]-deltas[i],
               starts[i]+deltas[i],
               starts[i]),
             y=yy, col=colors[i])
  }
}
