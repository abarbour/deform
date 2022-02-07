#' Predict the level of strain expected for magnitude and distance
#' 
#' @details 
#' Currently the options for \code{model} are:
#' \describe{
#'   \item{\code{'static'}}{Static RMS strain from an earthquake (Wyatt, 1988)}
#'   \item{\code{'dynamic'}}{Dynamic RMS strain from a seismic wave:
#'   \code{'aw14'} for Agnew and Wyatt (2014): applicable to teleseisms at thousands of km (inaccurate within 500km or so); or
#'   \code{'bc17'} for Barbour and Crowell (2017): applicable to regional earthquakes at < 500 km (inaccurate
#'   for very short distances); or
#'   \code{'fbl20'} for Farghal, Barbour and Langbein (2020): applicable to regional earthquakes 
#'   at < 500 km (best accuracy at short distances; inaccurate > 500 km);
#'   \code{'blf21'} for Barbour, Langbein, and Farghal (2021): best for regional earthquakes as
#'   it accounts for station and event terms in the fbl20 dataset; the basis of the M_{DS} 
#'   magnitude scale, which is basically equivalent to Mw.
#'   }
#' }
#' @name magnitude-distance
#' @param Mw numeric; the moment magnitude of the earthquake
#' @param Distance.km numeric; the epicentral distance, in kilometers
#' @param km2deg numeric; specify the scaling coefficient used to
#' convert the distance (in km) to degrees; this may be useful if,
#' for example, the distance is already in degrees (then set it equal to 1)
#' @param model character; the strain-scaling model to use (see Details)
#' @param strn.type character; the type of strain-scaling relationship to
#' provide a prediction for
#' 
#' @references Agnew, D. C., and F. K. Wyatt (2014),
#' Dynamic strains from regional and teleseismic earthquakes, 
#' Bull. Seismol. Soc. Am. 104, no. 4, 1846– 1859, 
#' doi: 10.1785/0120140007
#' 
#' @references Barbour, A., and B. Crowell (2017),
#' Dynamic strains for earthquake source characterization, 
#' Seismol. Res. Lett. 88, no. 2A, 354–370, 
#' doi: 10.1785/0220160155
#' 
#' @references Farghal, N., A. Barbour, and J. Langbein (2020),
#' The Potential of Using Dynamic Strains in Earthquake Early Warning Applications, 
#' Seismol. Res. Lett., in press, 1–12, 
#' doi: 10.1785/0220190385
#' 
#' @references Barbour, A., J. Langbein, and N. Farghal (2021),
#' Earthquake Magnitudes from Dynamic Strain,
#' Bull. Seismol. Soc. Am.m, 111 (3): 1325–1346
#' doi: 10.1785/0120200360
#' 
#' @references Wyatt, F. K. (1988),
#' Measurements of coseismic deformation in southern California: 1972–1982, 
#' J. Geophys. Res. 93, no. B7, 7923–7942, 
#' doi: 10.1029/JB093iB07p07923
#' 
#' @examples 
#' static(7.2, 220)
#' dynamic(7.2, 15000, model='aw14')
#' 
#' dynamic(7.2, 150, model='aw14')
#' dynamic(7.2, 150, model='bc17')
#' dynamic(7.2, 150, model='fbl20')
NULL

#' @rdname magnitude-distance
#' @export
static <- function(Mw, Distance.km, model="w88"){
  Distance.m <- Distance.km*1e3
  1.5*Mw - 3*log10(Distance.m) - 2.3 # equation 4
}

# Hanks and Kanamori
#' @export
Mw2M0 <- function(Mw) 10 ** (1.5 * Mw + 9.1) # result in N*m

#' @rdname magnitude-distance
#' @export
dynamic <- function(Mw, Distance.km, model=c("aw14",'bc17','fbl20','blf21'),
                    strn.type=c('general','areal','diffext','shear','dil'),
                    km2deg=NULL){
  
  model <- match.arg(model)
  Distance <- Distance.km
  
  if (is.null(km2deg)){
    km2deg <- 1/111.3195
  }
  if (model %in% c('aw14')){
    Distance <- Distance * km2deg
  }
  
  ld <- log10(Distance)
  
  lstrn <- if (model == 'blf21'){
    # Barbour Langbein and Farghal (2021)
    warning('BLF21 scaling result does not account for site or (longitude < -124) bias', 
            immediate. = FALSE, call. = FALSE)
    0.92*Mw - 1.45 * ld - 0.00072 * Distance + log10(3e-9)
  } else if (model == 'fbl20'){
    # Farghal, Barbour and Langbein (2020) Equation 11
    0.83*Mw - 1.43 * ld - 0.0013 * Distance - 8.051
  } else if (model == 'bc17'){
    # Barbour and Crowell (2017) Equation 8
    1.25*Mw - 1.992*ld - 9.513
  } else if (model == 'aw14'){
    # Agnew and Wyatt (2014)
    strn.type <- match.arg(strn.type)
    LS <- if (strn.type %in% c("areal",'diffext')){
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
    # they report in nanostrain, so correct to unscaled strain
    LS - 9
  } else {
    stop('bad model specification')
  }
  
  return(lstrn)
}