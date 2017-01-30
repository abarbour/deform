#' Deform a convex hull
#' @description
#' Deform the coordinates of a convex hull
#' @details
#' This effectively applies a 2D plane-strain tensor to the coordinates where the 
#' diagonal components equal \code{out.by} and the off diagonal components equal \code{shear.by}.
#' @export
#' @param hull.coords a matrix or data.frame with x and y coordinates of the convex hull
#' @param out.by numeric; the factor to expand the area of the hull by
#' @param shear.by numeric; the factor to shear the area of the hull by
#' @param plot logical; should the results be plotted?
#' @author A.J. Barbour
#' @examples
#' set.seed(1234)
#' X <- matrix(stats::rnorm(2000), ncol = 2)
#' hull.inds <- chull(X)
#' hull.inds <- c(hull.inds, hull.inds[1])
#' hull <- X[hull.inds, ]
#' 
#' deform_hull(hull, plot=TRUE)
#' deform_hull(hull, out.by=2, plot=TRUE)
#'
deform_hull <- function(hull.coords, out.by=1, shear.by=0, plot=FALSE){
  M <- as.matrix(hull.coords)
  stopifnot(ncol(M)>=2)
  cn <- names(hull.coords)
  cms <- colMeans(M)
  M.dem <- sweep(M, 2, cms)
  Affine <- matrix(c(out.by, shear.by, shear.by, out.by), ncol=2)
  Trans <- matrix(cms, ncol=2, nrow=nrow(hull.coords), byrow=TRUE)
  M.sc <- M.dem %*% Affine + Trans
  if (plot){
    plot(M.sc, col=NA)
    segments(Trans[,1], Trans[,2], M.sc[,1], M.sc[,2], col='grey', lty=2)
    points(cms[1], cms[2], pch="+", font=2, cex=1.5)
    lines(M, type='b', pch=16, cex=0.5)
    lines(M.sc, type='b', pch=16, col='red')
  }
  M.sc <- data.frame(M.sc)
  names(M.sc) <- cn
  M.sc
}