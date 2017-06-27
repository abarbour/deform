biot_willis <- function(nu, nuu, B){
  3*(nuu - nu)/(1 - 2*nu)/(1 + nuu)/B
}
compressibility <- function(nu, nuu, B, mu=1){
  9*(1 - 2*nuu)*(nuu - nu)/2/mu/(1 - 2*nu)/(1 + nuu)^2/B^2
}
darcy_cond <- function(nu, nuu, B, mu=1, D=1){
  9*D*(1 - nuu)*(nuu - nu)/2/mu/(1 - nu)/(1 + nuu)^2/B^2
}
lame <- function(nu, mu=1){
  2*mu*nu/(1 - 2*nu)
}

#' Rudnicki's response to point and step injection 
#' @param r,z numeric; receiver position (radial, depth-positive) [m]
#' @param t numeric; time [s]
#' @param zinj numeric; injection depth [m]
#' @param mu numeric; elastic shear modulus [Pa]
#' @param d numeric; hydraulic diffusivity [m^2/s]
#' @param nu,nuu numeric; drained and undrained Poisson's ratio [0,0.5]
#' @param B numeric; Skempton's coefficient [0,1]
#' @param response character; the type of reponse be the impulsive or heaviside
#' @export
#' @examples
#' r <- 1
#' zi <- 1
#' t <- seq(1,1000,length.out=101)
#' mu <- 15e9
#' 
#' # Impulse response
#' rt <- rudnicki(r,1,t,zi, mu=mu)
#' rtdeep <- rudnicki(r,10,t,zi, mu=mu)
#' 
#' # Step response
#' rt_s <- rudnicki(r,1,t,zi, mu=mu, response='step')
#' rtdeep_s <- rudnicki(r,10,t,zi, mu=mu, response='step')
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
rudnicki <- function(r, z, t, zinj=0, mu=10e9, d=1, nu=0.25, nuu=0.33, B=1, response=c('impulse','step')){
  
  alpha <- biot_willis(nu, nuu, B)
  beta <- compressibility(nu, nuu, B, mu)
  chi <- darcy_cond(nu, nuu, B, mu, d)
  la <- lame(nu, mu)
  
  Zrel <- z - zinj
  p0 <- 1/(4*pi*chi)
  u0 <- p0 * alpha/2/(la + 2*mu)
  
  Rsq <- r^2
  Dis <- sqrt(Zrel^2 + Rsq)
  Dis3 <- Dis^3
  
  eta <- Dis/sqrt(d*`t`)/2
  etasq <- eta^2
  eta3 <- eta^3
  eta4 <- eta^4
  
  eta_erf <- deform::erf(eta)
  eta_erfc <- deform::erfc(eta)
  
  fct  <- eta_erfc + eta_erf/etasq/2 - exp(-etasq)/eta/sqrt(pi)
  pfct <- -eta_erf/eta3 + 2*exp(-etasq)/etasq/sqrt(pi)
  
  perfc <- -2 * exp(-etasq)/sqrt(pi)
  ppfct <- perfc/eta3 + 3*eta_erf/eta4 - 4*exp(-etasq)*(1 + etasq)/sqrt(pi)/eta3
  doteta <- -eta/t/2
  
  response <- match.arg(response)
  parms <- sprintf("nu=%s, nu_u=%s, B=%s, mu=%s, D=%s", nu, nuu, B, mu, d)
  message("calculating ", response, " response for ", parms)
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
    
    dp <- p0*doteta*(perfc + 2*(1 - 2*etasq)*exp(-etasq)/sqrt(pi))
  
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
    
    dp <- p0*(eta_erfc + 2*eta*exp(-etasq)/sqrt(pi))

  }
  
  tlt <- -ezr
  dpz <- -dp*Zrel / Dis3
  dpr <- -dp*r / Dis3
  
  res <- list(type = response,
              params = list(r, z, zinj, 
                            alpha = alpha, beta = beta, chi = chi, lambda = la),
              time = `t`,
              displacement = cbind(z=uz, r=ur),
              strain = cbind(zz = ezz, rr = err, tt = ett, zr = ezr),
              tilt = cbind(zr=tlt),
              pore.pressure = cbind(p),
              darcy.flux = cbind(p=dp, pz=dpz, pr=dpr)
  )
  class(res) <- c('rudnicki','list')
  return(res)
}

#' @rdname rudnicki
#' @export
#' @method plot rudnicki 
plot.rudnicki <- function(x, response=FALSE, ...){
  
  typ <- x[['type']]
  
  xt <- as.vector(x[['time']])
  
  app <- ifelse(response, " (response)","")
  
  y <- x[['displacement']]
  ynms <- paste(colnames(y), collapse=" - ")
  if (!response) y <- apply(y, 2, cumsum)
  matplot(xt, y, type='l', main=paste0(typ, ': displacements', app), ...)
  mtext(ynms, line=-1)
  
  y <- x[['strain']]
  ynms <- paste(colnames(y), collapse=" - ")
  if (!response) y <- apply(y, 2, cumsum)
  matplot(xt, y, type='l', main=paste0(typ, ': strain', app), ...)
  mtext(ynms, line=-1)
  
  y <- x[['tilt']]
  ynms <- paste(colnames(y), collapse=" - ")
  if (!response) y <- apply(y, 2, cumsum)
  matplot(xt, y, type='l', main=paste0(typ, ': tilt', app), ...)
  mtext(ynms, line=-1)
  
  y <- x[['pore.pressure']]
  ynms <- paste(colnames(y), collapse=" - ")
  if (!response) y <- apply(y, 2, cumsum)
  matplot(xt, y, type='l', main=paste0(typ, ': pore pressure', app), ...)
  mtext(ynms, line=-1)
  
  y <- x[['darcy.flux']]
  ynms <- paste(colnames(y), collapse=" - ")
  if (!response) y <- apply(y, 2, cumsum)
  matplot(xt, y, type='l', main=paste0(typ, ': darcy flux', app), ...)
  mtext(ynms, line=-1)
}
