#!/usr/bin/env Rscript --no-save

###
#	test_lv_gps.R
#	/Users/abl/survey.processing/development/R/packages/deform/NOBUILD/Sandbox/LV
#	Created by 
#		/Users/abl/bin/ropen ( v. 2.6.9 )
#	on 
#		2019:151 (31-May)
#
#	[ Explain what this script does, broadly ]
#
###

## local functions
try(source('segall92.R'))

## libs

library(tools)

#if (!require("pacman")) install.packages("pacman", dependencies=TRUE)
#pacman::p_load(package1, package2, package_n)

# loads core tidy packages:  ggplot2, tibble, tidyr, readr, purrr, and dplyr
library(tidyverse)
#tidyverse_update(TRUE)

## local/github libs
# devtools::install_github("abarbour/kook")
#library(kook)
#Set1 <- brew.set1()
#Set1l <- brew.set1(TRUE)
#Dark2 <- brew.dark2()
#Dark2l <- brew.dark2(TRUE)
#Dark2ll <- brew.dark2(TRUE,TRUE)

#+++++++++++

shake <- FALSE
redo <- FALSE
inter <- interactive()

if (!exists("gps") | redo){
	gps <- readr::read_table('crossection.out')
}

#+++++++++++

# keep distances in km to stabilize integration
r <- 10**seq(0.9, 2, length.out=91)
Pa <- 80e3 #Pa (N / m^2)
# rescale pressure to maintain unit consistency
pressure <- 2.5 * Pa * (1000^2)
# Dimensions of inflating reservoir
disk.radius <- 10
disk.thickness <- 4 #15

redo.comp <- TRUE
if (!exists('res') | redo.comp) res <- surface_displacement_ringdisk(r, radius = disk.radius, thickness = disk.thickness, depth = disk.thickness / 2, p0 = pressure)

head(res)

#+++++++++++

FIG <- function(...){
	with(res,{
		uz_sc <- -uz / 1 #scaling
		ur_sc <- ur / 1 #scaling
		yl <- range(c(uz_sc, ur_sc), na.rm=TRUE)
		plot(x, uz_sc, col=NA, log='', ylim=yl, ylab="Displacement, mm", xlab="Axial Distance, km", main="Long Valley GPS")
		pu <- par('usr')
		#segments(0,c(0,1),c(disk.radius,disk.thickness), c(0,1), lwd=c(1,2))
		points(c(disk.radius,disk.thickness), pu[3]+1+c(0,1), pch="|",lwd=3)
		text(max(c(disk.radius,disk.thickness)), pu[3]+1+c(0,1), c('Radius','Thickness'), pos=4, font=3, cex=0.9)
		lines(x, uz_sc, col='red', lwd=1.5)
		lines(x, ur_sc, lty=2, col='blue', lwd=1.5)
	})
	legend('topright', c('Vertical','Radial'), lty=c(1,2), col=c('red','blue'))
	with(gps,{
		d <- `Distance(km)`
		x <- d + (disk.radius)
		H <- `Hor (mm)`
		Z <- `Vert (mm)`
		points(x, H, lwd=1.5, col='blue')
		points(x, Z, lwd=1.5, pch=3, col='red')
	})
	axis(1, at=1:20, labels=FALSE, col=NA, col.ticks=1, tcl=-0.2)
}

#+++++++++++

shake.fig <- TRUE
if (shake.fig & inter){
	FIG() 
} else {
#figfi <- "test_lv_gps"
#h <- w <- 7
#niceEPS(figfi, h=h, w=w, toPDF=TRUE)
#try(FIG())
#niceEPS()
#nicePNG(figfi, h=h, w=w)
#try(FIG())
#nicePNG()
}

#+++++++++++

#rdafi <- "test_lv_gps.rda"
#save(, file=rdafi)
#message("Saved ", rdafi)

