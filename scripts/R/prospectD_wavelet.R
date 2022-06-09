
# -------------------------------------------------------------------------------------------------------
#
# PROSPECT-D LUT generation and model inversion using Wavelet transformation
#
# Based on Alan, G., & Garke, J. (2008). Retrieval of chlorophyll concentration from leaf 
# reflectance spectra using wavelet analysis, 112, 1614-1632. http://doi.org/10.1016/j.rse.2007.08.005
#
# Author: Javier Lopatin
# javier.lopatin@uai.cl
# Universidad Adolfo Ibanez
# Santiago
# Chile
#
# -------------------------------------------------------------------------------------------------------


# load libraries
require(wmtsa)

# set variables for modeling and testing
#################################################################################
#no_cores = 4
# Optimal scale = 1,64 and 12 predictors (wavelet_no)
wavelet_yesno = 1
#number of wavelets
wavelet_no = 6
#scale range
wavelet_range = c(1,200) # 0.08163825

LUTsize = 1000 # 1,000 for testing and 10,000 for proper analysis;  Duan et al. 2011 / Atzberger 2012
resolution = 1000

percent = 1 #: percent of the simulated spectra with lowest cost function ( 20 % recommended by Vohland et al. 2010 and Duan et al. 2014)

# set working directory
setwd("~/Documentos/GitHub/RTM")

# load functions
# main.r contains the PROSPECT-D functions
source("~/Documentos/GitHub/RTM/scripts/R/PROSPECT_D/main.r")
# utils.r contains the LUT generation and inversion funcitons
source("~/Documentos/GitHub/RTM/scripts/R/PROSPECT_D/utils.r")

###################################################
# get measured spectra
spectra = as.matrix(read.csv("data/FAB2016_smoothed.csv", sep=",", header=TRUE))

# minband = min(as.numeric(substr(colnames(spectra), 2,5)))
minband = 470
maxband = max(as.numeric(substr(colnames(spectra), 2,5)), na.rm=TRUE)
spectra = spectra[,(minband-400+9):ncol(spectra)]

# test wavelet transform
wavelet_spectra = as.matrix(wmtsa::wavCWT(x=as.numeric(spectra[21,]), n.scale=wavelet_no, scale.range=wavelet_range))
par(mfrow=c(2,3))
plot(wavelet_spectra[,1], type="l")
plot(wavelet_spectra[,2], type="l")
plot(wavelet_spectra[,3], type="l")
plot(wavelet_spectra[,4], type="l")
plot(wavelet_spectra[,5], type="l")
plot(wavelet_spectra[,6], type="l")
par(mfrow=c(1,1))

################################################
#: LUT Inversion
################################################

#define parameters
N_minmax      = c(0.9, 2.4,resolution)   	 #% mesophyll / leaf structure coefficient
Cab_minmax    = c(2, 75, resolution)         #% chlorophyll content (�g.cm-2)
Car_minmax    = c(2,20, resolution)          #% carotenoid content (�g.cm-2)
Ant_minmax    = c(0.3,8, resolution)         #% Anth content (�g.cm-2)
Cbrown_minmax = c(0,0.4, resolution)		 #% brown pigment content (arbitrary units)
Cw_minmax     = c(0.004, 0.055, resolution)  #% EWT (cm) = equivalent water thickness
Cm_minmax     = c(0.0015,0.015, resolution)	 #% LMA (g.cm-2) = dry matter conent

# get list of PROSAIL parameters
parameters = list(LUTsize);
parameters[[1]] = N_minmax
parameters[[2]] = Cab_minmax
parameters[[3]] = Car_minmax
parameters[[4]] = Ant_minmax
parameters[[5]] = Cw_minmax
parameters[[6]] = Cm_minmax
parameters[[7]] = Cbrown_minmax

#: generate LUT
################################################
# list containing the spectra generated using the wavelet transform for scales=3,4,5,6
LUT1 = getProsailLUT(parameters=parameters, LUTsize=LUTsize)
names(LUT1)


# get inverted traits                                           
################################################
inversion = invertProsail(parameters, obtainedLUT=LUT1, measuredSpectra=spectra,
                                 LUTsize=LUTsize, percent = 1)                                       
names(inversion)                                          
inversion$estimates_


################################################
#: export results
################################################

# export plot with histograms for each trait
png(file = "mn_traithistogram.png",width = 800,units="px", height = 1200, res = 210,bg = "white")
par(mfrow=c(3,2))
for (i in 1:7){
  hist(estimates_traits[,i], main=colnames(estimates_traits)[i], xlab="")
}
par(mfrow=c(1,1))
dev.off()

write.table(estimates_traits, file="inversion_PD_estimates_traits.txt")
write.table(estimates_sd, file="inversion_PD_estimates_sd.txt")
write.table(estimates_meanRMSE, file="inversion_PD__estimates_meanRMSE.txt")
write.table(estimates_spectra, file="inversion_PD_estimates_spectra.txt")
