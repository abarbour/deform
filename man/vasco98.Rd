% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vasco98.R
\name{vasco98}
\alias{vasco98}
\alias{.surface_g_pt}
\alias{surface_displacement_point}
\title{Surface deformation associated with fluid withdrawl: an alternative formulation}
\usage{
vasco98(help = FALSE)

.surface_g_pt(x = 0, x_src = 0, z_src = 0, nuu. = 1/3)

surface_displacement_point(x, depth, delV., B. = 1, rho_f. = 1000, tol, ...)
}
\arguments{
\item{help}{logical; load documentation for \code{\link{vasco98}}}

\item{x}{numeric; spatial coordinate relative to extraction point}

\item{x_src}{numeric; the horizontal distance from the source}

\item{z_src}{numeric; the depth of the source below the surface}

\item{nuu.}{numeric; the 'undrained' Poisson's ratio (typically 1/3)}

\item{depth}{numeric; the depth below the surface to the source}

\item{delV.}{numeric; uniform change in pore-fluid volume (positive = loss)}

\item{B.}{numeric; the xxx}

\item{rho_f.}{numeric; the fluid density}

\item{tol}{numeric; the numerical tolerance in the integration; if 
supplied, this should be much smaller than the smallest difference between
any values in \code{x} divided by the depth}

\item{...}{additional arguments passed to \code{\link{.surface_g_pt}}}
}
\description{
Surface deformation associated with fluid withdrawl: an alternative formulation
}
\examples{
# Reproduce...
r <- seq(-10,10,by=0.2) # km
r.m <- r*1e3
xx <- surface_displacement_point(r.m, depth=2e3, delV.=1e7)
plot(uz ~ x, xx, type='l')
}
\references{
Vasco, D., et al. (1998), Monitoring of Fluid Injection and 
Soil Consolidation Using Surface Tilt Measurements,
\url{http://ascelibrary.org/doi/abs/10.1061/(ASCE)1090-0241(1998)124:1(29)}
}
\seealso{
\code{\link{segall85}}, or \code{\link{segall89}}; \code{\link{Simple-deformation}}
}
