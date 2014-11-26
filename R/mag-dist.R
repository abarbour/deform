#' Predict the level of strain expected for magnitude and distance
#' @name magnitude-distance
NULL

#' @rdname magnitude-distance
static <- function(Mw, Distance.km, model="w88"){
  Distance.m <- Distance.km*1e3
  1.5*Mw - 3*log10(Distance.m) - 2.3 - 9
}

#' @rdname magnitude-distance
dynamic <- function(Mw, Distance.km, km2deg=NULL, model="aw14"){
  if (is.null(km2deg)){
    km2deg <- 1/111.3195
  }
  Distance.deg <- Distance.km * km2deg
  0.95*Mw - 1.65*log10(Distance.deg) - 2.8 - 9
}