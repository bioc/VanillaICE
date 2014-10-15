test_sweepMode.R <- function(){
  data(snp_exp)
  x <- lrr(snp_exp)
  stat <- function(x, MARGIN=c(1,2)) {
    MARGIN <- match.arg(MARGIN)
    if(MARGIN==1){
      y <- rowModes(x)
    } else y <- colModes(x)
    y
  }
  xrows <- sweep(x, 1, rowModes(x))
  xcols <- sweep(x, 2, colModes(x))
  xrows2 <- sweepMode(snp_exp, 1)
  xcols2 <- sweepMode(snp_exp, 2)
  checkTrue(all.equal(xrows, lrr(xrows2)))
  checkTrue(all.equal(xcols, lrr(xcols2)))
}

test_duplicatedMapLocs <- function(){
  data(snp_exp)
  g <- rowData(snp_exp)
  ## duplicate the first row
  g[length(g)] <- g[1]
  rowData(snp_exp) <- g
  snp_exp2 <- dropDuplicatedMapLocs(snp_exp)
  checkIdentical(nrow(snp_exp2), nrow(snp_exp)-1L)
}
