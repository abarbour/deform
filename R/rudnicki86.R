#' Rudnicki's response to point and step injection 
#' 
#' @param r,z numeric; receiver position (radial, depth-positive) [m]
#' @param t numeric; time [s]
#' @param zinj numeric; injection depth [m]
#' @inheritParams effstress
#' @param response character; the type of response to output: 
#' either the points-source impulsive (\code{'impulse'}) or Heaviside (\code{'step'}) response
#' @param ... additional arguments
#' @param x an object with class \code{'rudnicki.pt'}, e.g. from \code{\link{rudnicki86}}
#' 
#' @references 
#' Rudnicki, J. W. (1986), Fluid mass sources and point forces in linear elastic diffusive solids, 
#' \emph{Mechanics of Materials}, 5(4), 383-393, \url{https://doi.org/10.1016/0167-6636(86)90042-6}
#'
#' @export
#' 
#' @examples
#' r <- 1
#' zi <- 1
#' t <- 10^seq(-4, 1, length.out=301)
#' 
#' # Impulse response
#' rt <- rudnicki86(r,1,t,zi)
#' rtdeep <- rudnicki86(r,10,t,zi)
#' 
#' # Step response
#' rt_s <- rudnicki86(r,1,t,zi, response='step')
#' rtdeep_s <- rudnicki86(r,10,t,zi, response='step')
#' 
#' layout(matrix(1:10, 2, byrow=TRUE))
#' plot(rt, col=1)
#' plot(rt, response=TRUE, col=1)
#' 
#' plot(rtdeep, col=2)
#' plot(rtdeep, response=TRUE, col=2)
#' 
#' plot(rt_s, col=3)
#' plot(rt_s, response=TRUE, col=3)
#' 
#' plot(rtdeep_s, col=4)
#' plot(rtdeep_s, response=TRUE, col=4)
#'
rudnicki86 <- function(r, z, t, 
                     zinj=0, mu=1, diffusiv=1, nu=0.25, nuu=0.33, B=1, 
                     response=NULL){
  
  response <- match.arg(response, c('impulse','step'))
  check_range(B, 0, 1)
  check_range(nu, 0, 0.5)
  check_range(nuu, 0, 0.5)
  stopifnot(diffusiv > 0)
  stopifnot(mu > 0)
  
  parms <- sprintf("nu=%s, nu_u=%s, B=%s, mu=%s, D=%s", nu, nuu, B, mu, diffusiv)
  message("calculating ", response, " response for ", parms)
  
  alpha <- effstress(nu, nuu, B)
  beta <- biot_compressibility(nu, nuu, B, mu)
  chi <- darcy_conductivity(nu, nuu, B, mu, diffusiv)
  la <- lame_first(nu, mu)

  ptres <- .rudnicki_pt(Z=z, Zinj=zinj, R=r, Time=`t`, Diffusiv = diffusiv, 
                        Alpha = alpha, Beta = beta, Chi = chi, Lambda = la, Mu = mu, 
                        impulse=response == 'impulse')
  
  ptres[['params']][['input']] <- list(nu, nuu, B, diffusiv)

  class(ptres) <- c('rudnicki.pt','list')
  return(ptres)
}

.rudnicki_pt <- function(Z, Zinj, R, Time,
                         Diffusiv, Alpha, Beta, Chi, Lambda, Mu,
                         impulse = TRUE) {
    Zrel <- Z - Zinj
    p0 <- 1 / 4 / pi / Chi
    u0 <- p0 * Alpha / 2 / (Lambda + 2 * Mu)
    
    Rsq <- R^2
    Dis <- sqrt(Zrel^2 + Rsq)
    Dis3 <- Dis^3
    
    eta <- Dis / 2 / sqrt(Diffusiv * Time)
    dot_eta <- -eta / Time / 2
    etasq <- eta^2
    eta3 <- eta^3
    eta4 <- eta^4
    
    eta_erf <- erf_re(eta)
    eta_erfc <- erfc_re(eta)
    
    sqpi <- sqrt(pi)
    expeta <- exp(-etasq)
    
    fct  <- eta_erfc + eta_erf / etasq / 2 - expeta / eta / sqpi
    pfct <- -eta_erf / eta3 + 2 * expeta / etasq / sqpi
    
    perfc <- -2 * expeta / sqpi
    ppfct <- perfc / eta3 + 3 * eta_erf / eta4 - 4 * expeta * (1 + etasq) / sqpi / eta3
    
    if (impulse) {
      #
      # response to impulsive source
      #
      C <- dot_eta / Dis
      uz <- C * u0 * pfct * Zrel
      ur <- C * u0 * pfct * R
      p <- C * p0 * perfc
      
      ezz <- u0 * (pfct * Rsq + (ppfct * eta + pfct) * Zrel^2) * dot_eta / Dis3
      err <- u0 * ((ppfct * eta + pfct) * Rsq + pfct * Zrel^2) * dot_eta / Dis3
      ett <- C * u0 * pfct
      ezr <- u0 * eta * ppfct * dot_eta * Zrel * R / Dis3
      
      dp <- p0 * dot_eta * (perfc + 2 * (1 - 2 * etasq) * expeta / sqpi)
      
    } else {
      #
      # response to heaviside source
      #
      
      uz <- u0 * fct * Zrel / Dis
      ur <- u0 * fct * R / Dis
      p <- p0 * eta_erf / Dis
      
      ezz <- u0 * (fct * Rsq  +  pfct * eta * Zrel^2) / Dis3
      err <- u0 * (pfct * eta * Rsq  +  fct * Zrel^2) / Dis3
      ett <- u0 * fct / Dis
      ezr <- u0 * (eta * pfct - fct) * Zrel * R / Dis3
      
      dp <- p0 * (eta_erfc + 2 * eta * expeta / sqpi)
      
    }
    
    tlt <- -ezr
    dpz <- -dp * Zrel / Dis3
    dpr <- -dp * R / Dis3
    
    res <- list(
      is.impulse = impulse,
      params = list(
        spatial = list(R, Z, Zinj),
        params = list(Diffusiv, Alpha, Beta, Chi, Lambda, Mu)
      ),
      time = Time,
      displacement = cbind(z = uz, r = ur),
      strain = cbind(
        zz = ezz,
        rr = err,
        tt = ett,
        zr = ezr
      ),
      tilt = cbind(zr = tlt),
      pore.pressure = cbind(p),
      darcy.flux = cbind(p = dp, pz = dpz, pr = dpr)
    )
    return(res)
  }

#' @rdname rudnicki86
#' @export
#' @method plot rudnicki.pt
plot.rudnicki.pt <- function(x, response=FALSE, ...){
  
  typ <- ifelse(x[['is.impulse']], "Impulse", "Step")
  
  xt <- as.vector(x[['time']])
  
  app <- ifelse(response, "\n(response)", "")
  
  .addleg <- function(yn){
    ny <- length(yn)
    if (ny>1) legend('top', yn, ncol = ny, col=seq_len(ny), lty=seq_len(ny), bty='n')
  }
  
  y <- x[['displacement']]
  if (!response) y <- apply(y, 2, cumsum)
  matplot(xt, y, type='l', 
          ylab="", xlab="",
          main=paste0(typ, ': displacements', app), ...)
  .addleg(colnames(y))
  
  y <- x[['strain']]
  if (!response) y <- apply(y, 2, cumsum)
  matplot(xt, y, type='l', 
          ylab="", xlab="",
          main=paste0(typ, ': strain', app), ...)
  .addleg(colnames(y))
  
  y <- x[['tilt']]
  if (!response) y <- apply(y, 2, cumsum)
  matplot(xt, y, type='l', 
          ylab="", xlab="",
          main=paste0(typ, ': tilt', app), ...)
  .addleg(colnames(y))
  
  y <- x[['pore.pressure']]
  if (!response) y <- apply(y, 2, cumsum)
  matplot(xt, y, type='l', 
          ylab="", xlab="",
          main=paste0(typ, ': pore pressure', app), ...)
  .addleg(colnames(y))
  
  y <- x[['darcy.flux']]
  if (!response) y <- apply(y, 2, cumsum)
  matplot(xt, y, type='l', 
          ylab="", xlab="",
          main=paste0(typ, ': Darcy fluid flux', app), ...)
  .addleg(colnames(y))
}
