#' Surface deformation associated with fluid withdrawl
#' @name segall85
#' @aliases segall
#' @export
#' @param help logical; load documentation for \code{\link{segall85}}
#' @seealso \code{\link{Simple-deformation}}
#' @examples
#' \dontrun{
#' x. <- c(-10:-2, seq(-1.9,1.9,by=0.1), 2:10)
#' su <- surface_displacement(x.*1e3, z_src=1e3)
#' plot(uz ~ x, su, type="b")
#' }
segall85 <- function(help=FALSE){
  cat("
This is just a placeholder function. See documentation (`?segall85`).
      ")
  if (help) ?segall85
}

#  spatial coordinates in segall
#   x is positive downwards (here z)
#   y is positive away from source (here x)
#   z is along line source (here y)
#  source coordinates:
#   z_src is zeta (depth)
#   y_src is xi (distance away)
#  parameters:
#   mu shear modulus
#   C units of force, proportional to the source strength
#' @param x numeric; spatial coordinate relative to extraction point
#' @param C. numeric; units of force, proportional to source strength (e.g., extraction rate)
#' @param mu. numeric; the shear modulus in Pascals
#' @param ... additional arguments to \code{\link{.surface_g}}
#' @rdname segall85
#' @export
surface_displacement <- function(x, C.=1, mu.=1e9, ...){
  sg <- .surface_g(x, ...)
  c. <- C./mu.
  sg <- mutate(sg,
               ux = gx * c.,
               uz = gz * c.,
               uxz.mag = sqrt(ux^2 + uz^2),
               uxz.ang = atan2(uz,ux) * 180/pi )
  return(sg)
}

#' @rdname segall85
#' @export
#' @param x_src numeric; the horizontal distance from the source
#' @param z_src numeric; the depth of the source below the surface
#' @param nuu numeric; the 'undrained' Poisson's ratio (typically 1/3)
.surface_g <- function(x=0, x_src=0, z_src=0, nuu=1/3){
  # segall85 C9
  gx <-  2*(1 - nuu)*(x - x_src)/(z_src^2 + (x - x_src)^2)
  gz <- -2*(1 - nuu)*z_src/(z_src^2 + (x - x_src)^2)
  data.frame(x, gx, gz, xz=x/z_src)
}

#' @rdname segall85
#' @export
#' @param Time numeric; the time from xxx
#' @param Vdot. numeric; the volumetric flow rate xxx
#' @param L. numeric; the xxx
#' @param t. numeric; the xxx
#' @param HD. numeric; the xxx
#' @param phi. numeric; the xxx
timevarying_fluidmass <- function(x, Time, Vdot., L., t., HD., phi.){
  #
  # [ ] account for multiple Vdot.
  #
  FUN <- function(ti){
    -1 * Vdot. * sqrt(ti / HD.) * ierfc2(sqrt(x^2 / (4 * HD. * ti))) / (L. * t. * phi.)
  }
  apply(matrix(Time), 1, FUN)
}

#' @rdname segall85
#' @export
#' @param B. numeric; the xxx
#' @param D. numeric; the xxx
#' @param Pt.Sources.x numeric; a vector of point-source locations in the x direction
#' @param x.lim numeric; the limit of integration in both the positive and negative directions; if 
#' missing this is based on the absolute maximum of \code{x}
timevarying_surface_displacement <- function(x, Time, Vdot., B., L., D., HD., nuu.=1/3, Pt.Sources.x=0, x.lim){
  # Time varying deformation associated with fluid extraction
  # segall85 eq 26
  #
  x.lim <- if (missing(x.lim)){
    2 * round(max(abs(x), na.rm=TRUE))
  } else {
    as.numeric(x.lim)
  }
  #
  Sources <- matrix(Pt.Sources.x, ncol=1)
  #
  .FUN <- function(ti, Xi.sources){
    # sum over source contributions
    .fun <- function(Xi., dXi.=0, X.=x, ti.=ti){
      # from the spatio-temporal component
      ti.d <- Xi.^2 / (4 * HD. * ti.)
      ierfc2(sqrt(ti.d))/((D.)^2 + (X. - Xi. - dXi.)^2)
    }
    #
    # integrate over Xi
    .xi.integ <- function(dxi){
      apply(matrix(x), 1, FUN=function(xx){
        integrate(.fun, -1*x.lim, x.lim, X.=xx, dXi.=dxi, 
                  subdivisions = 50L, 
                  rel.tol=.Machine$double.eps^0.35)$value
      })
    }
    srci <- seq_along(Xi.sources)
    xvar.t2 <- apply(Xi.sources, MARGIN=1, FUN=.xi.integ)
    # time varying component (scaling)
    tvar <- -2 * B. * (1 + nuu.) * Vdot. * D. * sqrt(ti / HD.) / (3 * pi * L.)
    #
    nctv2 <- ncol(xvar.t2)
    sc <- if (length(tvar) == nctv2){
      # accounts for multiple values of Vdot.
      matrix(rep(tvar, nrow(xvar.t2)), ncol=nctv2, byrow = TRUE)
    } else {
      if (length(tvar)==1){
        tvar
      } else {
        warning("Number of values for Vdot is not 1, or not equal to the number of sources.\nUsed first value.")
        tvar[1]
      }
    }
    # combine and return
    res <- rowSums(xvar.t2 * sc)
    return(res)
  }
  # apply FUN through time
  apply(X=matrix(Time), MARGIN=1, FUN=.FUN, Xi.sources=Sources)
}

#' @rdname segall85
#' @export
#' @param nu. numeric; the drained Poisson's ratio (typically 1/4)
#' @param mu.gpa. numeric; the shear modulus in giga-Pascals (GPa)
timevarying_porepressure <- function(x, z, Time, Vdot., B., L., D., HD., t., mu.gpa., nu.=1/4, nuu.=1/3, Pt.Sources.x=0, x.lim){
  #
  # Time varying p.p. associated with fluid extraction
  # segall85 eq 28
  #
  x.lim <- if (missing(x.lim)){
    round(max(abs(x), na.rm=TRUE))
  } else {
    as.numeric(x.lim)
  }
  #
  sc <- 1e9
  mu. <- sc * mu.gpa.
  #
  Sources <- matrix(Pt.Sources.x, ncol=1)
  #
  .FUN <- function(ti, zi=0, Xi.sources){
    message(paste("Time:", ti/365/86400, "years"))
    #
    # Function that returns the integral value at
    # a position Xi away from the source at dXi
    .fun <- function(Xi., dXi.=0, X.=x, Z.=0, ti.=ti){
      # from the spatio-temporal component
      ti.d <- Xi.^2 / (4 * HD. * ti.)
      C1 <- ierfc2(sqrt(ti.d))
      #
      xmxi <- (X. - Xi. - dXi.)^2
      zpd <- (Z. + D.)
      C2A <- zpd / (zpd^2 + xmxi^2)
      zpdt <- zpd + t.
      C2B <- -1 * zpdt / (zpdt^2 + xmxi^2)
      sc * (C2A + C2B) * C1
    }
    #
    # function used to integrate over Xi at different sources
    .xi.integ2D <- function(dxi, ti.=ti){
      message(paste("\t2D integration of source at:", dxi/1e3, "km"))
      grd <- as.matrix(expand.grid(VarX=x, VarZ=z))
      ires <- apply(grd, 1, FUN=function(xxzz){
        #message(paste(xxzz, collapse="//"))
        zi <- try(integrate(.fun, -1*x.lim, x.lim, X.=xxzz[1], Z.=xxzz[2], dXi.=dxi, subdivisions = 150L, rel.tol=.Machine$double.eps^0.5))
        if (inherits(zi,"try-error")){
          #message(paste("\t==== ==== ====>  shrink integration limits to below",x.lim))
          message(paste(xxzz, collapse=" x // z "))
          NA
        } else {
          #message(paste("integration",zi$message))
          ti.d <- xxzz[1]^2 / (4 * HD. * ti.)
          C1 <- (1 - 2*nu.) * ierfc2(sqrt(ti.d)) / ((nuu. - nu.) * (1 - 2*nuu.))
          C2 <- -2/(pi * (1 - nuu.))
          C1 + C2 * zi$value
        }
      })
      matrix(ires, ncol=length(z), byrow = FALSE)
    }
    # Array of depth slices for all sources
    sourcePP <- abind(lapply(X = Sources, FUN = .xi.integ2D), along=3)
    # time varying component (scaling)
    tvar <- -2 * mu. * (1 + nuu.)^2 * B.^2 * Vdot. * sqrt(ti / HD.) / (9 * L. * t.)
    #
    nsrc <- length(Sources)
    vsc <- if (length(tvar) == nsrc){
      # accounts for multiple values of Vdot.
      dsp <- dim(sourcePP)
      #array(data = tvar, dim = dim(sourcePP))
      abind(lapply(tvar, function(tv) array(tv, c(dsp[1:2],1))))
    } else {
      if (length(tvar)==1){
        tvar
      } else {
        warning("Number of values for Vdot is not 1, or not equal to the number of sources.\nUsed first value.")
        tvar[1]
      }
    }
    # multiple source solution array by tvar (which will be an equal-size array is n_vdot == n_sources)
    sourcePP <- vsc * sourcePP
    # superpostion of all sources
    res <- apply(X = sourcePP, MARGIN = c(1,2), FUN = sum, na.rm = TRUE)
    return(-1*res/sc)
  }
  # apply FUN through time
  abind(lapply(X = Time, FUN = .FUN, Xi.sources=Sources), along=3)
}

#' Simple numerical deformation estimates
#' 
#' Calculate tilts and extensions based on spatially varying displacements
#' 
#' @details
#' \code{\link{Uniaxial_extension}} calculates the component of
#' strain associated with deformation along a single axis.
#' For example, the change in the Eastward displacements (or rates)
#' in the East direction.
#' 
#' \code{\link{Tilt}} calculates the tilt field associates with
#' spatial variations in vertical positions (or rates of change);
#' the sign convention used is such that
#' a ball placed on the x-y plane would roll in the direction of the
#' tilt vector.  Or, in other words, the tilt vector is the direction a
#' plumb bob would move.
#' 
#' Calculations are done with \code{\link{diff}}.
#' 
#' @note \code{\link{Tilt}} has not been well-tested for two-dimensional
#' results!
#' @name Simple-deformation
#' @param x numeric; the spatial vector in real units (not an index!)
#' @param X numeric; the quantity varying in the direction of \code{x}
#' @param y numeric; the spatial vector perpendicular to \code{x} in real units
#' @param z numeric; values which vary in the direction perpendicular to the \code{x}-\code{y} plane
#' @param left logical; should padding be on the "left" side of the vector (i.e., index number 1)?
#' @export
#' @seealso \code{\link{segall85}}
#' @examples
#' # Some x values
#' xval. <- sort(unique(c(-7:-3, seq(-3.90,3.90,by=0.05), 3:7)))
#' set.seed(1221)
#' xanis <- rnorm(length(xval.), sd = 10)
#' 
#' #  surface displacements
#' su <- surface_displacement(xval.*1e3, C.=1e13, z_src=0.7e3)
#' 
#' # Vertical tilt -- assumes axial symmetry
#' sut <- with(su, Tilt(x, z=uz))
#' #               -- including anisotropic effects
#' sut.anis <- with(su, Tilt(x, x = sort(x + xanis), z=uz))
#' 
#' plot(ztilt ~ x, sut.anis, col='blue', pch=16, cex=0.5)
#' lines(ztilt ~ x, sut, lwd=2)
#' 
#' plot(ztilt ~ abs(x), sut.anis, col='blue', pch=ifelse(sign(x)==1,16,1), cex=0.5)
#' lines(ztilt ~ abs(x), sut, lwd=2)
#'  
#' # Uniaxial strain in the 'x' direction
#' sue <- with(su, Uniaxial_extension(x, X=ux))
#' sue.anis <- with(su, Uniaxial_extension(sort(x + xanis), X=ux))
#' 
#' plot(dXdx ~ x, sue.anis, col='blue', pch=16, cex=0.5)
#' lines(dXdx ~ x, sue, lwd=2)
#'  
#' plot(dXdx ~ abs(x), sue.anis, col='blue', pch=ifelse(sign(x)==1,16,1), cex=0.5)
#' lines(dXdx ~ abs(x), sue, lwd=2)
#' 
#' # see how the 'left' argument affects things:
#' .setleft(1,TRUE) # NA  1
#' .setleft(1,FALSE) # 1  NA
Uniaxial_extension <- function(x, X, left=TRUE){
  deriv <- diff(X)/diff(x)
  dXdx <- .setleft(deriv, left)
  data.frame(x, dXdx)
}
#' @rdname Simple-deformation
#' @export
Tilt <- function(x, y=NULL, z, left=TRUE){
  deriv1 <- diff(z)/diff(x)
  dzdx <- .setleft(deriv1, left)
  dzdy <- if (is.null(y)){
    dzdx
  } else {
    # this needs to be more intelligent:
    # what if z is a matrix!?
    deriv2 <- diff(z)/diff(x)
    .setleft(deriv2, left)
  }
  # reverse sign so that the sign convention
  # is such that a positive tilt corresponds to
  # the direction a plumb bob would head, or
  # the direction a ball on the surface would roll
  tlt <- -1*(dzdx + dzdy)
  ang <- atan2(dzdy, dzdx) * 180/pi 
  data.frame(x = x, ztilt = tlt, xy.direction = ang)
}
#' @rdname Simple-deformation
#' @export
.setleft <- function(x, left=TRUE){
  x <- as.vector(x)
  if (left){
    c(NA, x)
  } else {
    c(x, NA)
  }
}