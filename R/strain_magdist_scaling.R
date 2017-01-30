#' Predict the level of strain expected for magnitude and distance
#' @name magnitude-distance
#' @param Mw numeric; the moment magnitude of the earthquake
#' @param Distance.km numeric; the epicentral distance, in kilometers
#' @param km2deg numeric; specify the scaling coefficient used to
#' convert the distance (in km) to degrees; this may be useful if,
#' for example, the distance is already in degrees (then set it equal to 1)
#' @param model character; the strain-scaling model to use; currently ignored
#' because there are no alternatives
#' @param strn.type character; the type of strain-scaling relationship to
#' provide a prediction for
NULL

#' @rdname magnitude-distance
#' @export
static <- function(Mw, Distance.km, model="w88"){
  Distance.m <- Distance.km*1e3
  1.5*Mw - 3*log10(Distance.m) - 2.3 - 9
}

#' @rdname magnitude-distance
#' @export
dynamic <- function(Mw, Distance.km, km2deg=NULL, model="aw14",
                    strn.type=c('general','areal','diffext','shear','dil')){
  if (is.null(km2deg)){
    km2deg <- 1/111.3195
  }
  Distance.deg <- Distance.km * km2deg
  ld <- log10(Distance.deg)
  strn.type <- match.arg(strn.type)
  lstrn <- if (strn.type %in% c("areal",'diffext')){
    # areal or differential extension
    Err <- 10**(0.93*Mw - 1.63*ld - 2.66)
    Ett <- 10**(0.86*Mw - 1.73*ld - 2.46)
    strn <- if (strn.type=="areal"){
      Ett + Err
    } else {
      Ett - Err
    }
    log10(strn)
  } else if (strn.type=="shear"){
    # shear
    0.92*Mw - 1.77*ld - 2.64
  } else if (strn.type=="dil"){
    # dilatation
    0.94*Mw - 1.64*ld - 2.88
  } else {
    # general relation
    0.95*Mw - 1.65*ld - 2.8
  }
  return(lstrn - 9)
}