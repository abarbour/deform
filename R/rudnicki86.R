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
#' t <- seq(1,1000,length.out=101)
#' 
#' # Impulse response
#' rt <- rudnicki86(r,1,t,zi)
#' rtdeep <- rudnicki86(r,10,t,zi)
#' 
#' # Step response
#' rt_s <- rudnicki86(r,1,t,zi, response='step')
#' rtdeep_s <- rudnicki86(r,10,t,zi, response='step')
#' 
#' layout(matrix(1:6,2))
#' plot(rt, col=1)
#' plot.new()
#' 
#' plot(rtdeep, col=2)
#' plot.new()
#' 
#' plot(rt_s, col=3)
#' plot.new()
#' 
#' plot(rtdeep_s, col=4)
#' plot.new()
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
  
  Zrel <- z - zinj
  p0 <- 1/(4*pi*chi)
  u0 <- p0 * alpha/2/(la + 2*mu)
  
  Rsq <- r^2
  Dis <- sqrt(Zrel^2 + Rsq)
  Dis3 <- Dis^3
  
  eta <- Dis/sqrt(diffusiv * `t`)/2
  etasq <- eta^2
  eta3 <- eta^3
  eta4 <- eta^4
  
  eta_erf <- deform::erf(eta)
  eta_erfc <- deform::erfc(eta)
  
  sqpi <- sqrt(pi)
  
  fct  <- eta_erfc + eta_erf/etasq/2 - exp(-etasq)/eta/sqpi
  pfct <- -eta_erf/eta3 + 2*exp(-etasq)/etasq/sqpi
  
  perfc <- -2 * exp(-etasq)/sqpi
  ppfct <- perfc/eta3 + 3*eta_erf/eta4 - 4*exp(-etasq)*(1 + etasq)/sqpi/eta3
  doteta <- -eta/t/2
  
  if (response == 'impulse'){
    #
    # response to impulsive source
    #
    C <- doteta/Dis
    uz <- C*u0*pfct*Zrel
    ur <- C*u0*pfct*r
    p <- C*p0*perfc
    
    ezz <- u0*(pfct*Rsq + (ppfct*eta + pfct)*Zrel^2)*doteta / Dis3
    err <- u0*(pfct*Zrel^2 + (ppfct*eta + pfct)*Rsq)*doteta / Dis3
    ett <- C*u0*pfct
    ezr <- u0*eta*ppfct*doteta*Zrel*r / Dis3
    
    dp <- p0*doteta*(perfc + 2*(1 - 2*etasq)*exp(-etasq)/sqpi)
  
  } else {
    #
    # response to heaviside source
    #
    
    uz <- u0*fct*Zrel/Dis
    ur <- u0*fct*r/Dis
    p <- p0 * eta_erf / Dis
    
    ezz <- u0*(fct*Rsq + pfct*eta*Zrel^2) / Dis3
    err <- u0*(fct*Zrel^2 + pfct*eta*Rsq) / Dis3
    ett <- u0*fct / Dis
    ezr <- u0*(eta*pfct - fct)*Zrel*r / Dis3
    
    dp <- p0*(eta_erfc + 2*eta*exp(-etasq)/sqpi)

  }
  
  tlt <- -ezr
  dpz <- -dp*Zrel / Dis3
  dpr <- -dp*r / Dis3
  
  res <- list(type = response,
              params = list(spatial=list(r, z, zinj),
                            params=list(nu, nuu, B, diffusiv),
                            calc_params=list(alpha, beta, chi, lambda = la)),
              time = `t`,
              displacement = cbind(z=uz, r=ur),
              strain = cbind(zz = ezz, rr = err, tt = ett, zr = ezr),
              tilt = cbind(zr=tlt),
              pore.pressure = cbind(p),
              darcy.flux = cbind(p=dp, pz=dpz, pr=dpr)
  )
  class(res) <- c('rudnicki.pt','list')
  return(res)
}

#' @rdname rudnicki86
#' @export
#' @method plot rudnicki.pt
plot.rudnicki.pt <- function(x, response=FALSE, ...){
  
  typ <- x[['type']]
  
  xt <- as.vector(x[['time']])
  
  app <- ifelse(response, " (response)","")
  
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
