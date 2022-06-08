
# ------------------------------------------------------------------------------
#
#
#
#
#
#
# ------------------------------------------------------------------------------



# load libraries

require(wmtsa)
require(foreach)
require(doParallel)

# set variables for modeling
#################################################################################
no_cores = 4
#Alan, G., & Garke, J. (2008). Retrieval of chlorophyll concentration from leaf reflectance spectra using wavelet analysis, 112, 1614-1632. http://doi.org/10.1016/j.rse.2007.08.005
# Optimal scale = 1,64 and 12 predictors (wavelet_no)
wavelet_yesno = 1
#number of wavelets
wavelet_no = 6
#scale range
wavelet_range = c(1,200) # 0.08163825

LUTsize = 1000 # 1,000 for testing and 100,000 for proper analysis;  Duan et al. 2011 / Atzberger 2012
resolution = 1000

percent = 1 #: percent of the simulated spectra with lowest cost function ( 20 % recommended by Vohland et al. 2010 and Duan et al. 2014)

# set working directory
setwd("~/Documentos/GitHub/RTM")

# load functions
source("scripts/PROSPECT_D/main.r")


spectra = as.matrix(read.csv("data/FAB2016_smoothed.csv", sep=",", header=TRUE))
head(spectra)
#minband = min(as.numeric(substr(colnames(spectra), 2,5)))
minband = 470
maxband = max(as.numeric(substr(colnames(spectra), 2,5)), na.rm=TRUE)
spectra = spectra[,(minband-400+9):ncol(spectra)]
head(spectra)

#test
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
N_minmax      = c(0.9, 2.4,resolution)   	   #% mesophyll / leaf structure coefficient
Cab_minmax    = c(2, 75, resolution)         #% chlorophyll content (�g.cm-2)
Car_minmax    = c(2,20, resolution)          #% carotenoid content (�g.cm-2)
Ant_minmax    = c(0.3,8, resolution)         #% Anth content (�g.cm-2)
Cbrown_minmax	=	c(0,0.4, resolution)		     #% brown pigment content (arbitrary units)
Cw_minmax     = c(0.004, 0.055, resolution)  #% EWT (cm) = equivalent water thickness
Cm_minmax     = c(0.0015,0.015, resolution)	 #% LMA (g.cm-2) = dry matter conent


#: generate LUT
################################################

parameters = list();
parameters[[1]] = N_minmax
parameters[[2]] = Cab_minmax
parameters[[3]] = Car_minmax
parameters[[4]] = Ant_minmax
parameters[[5]] = Cw_minmax
parameters[[6]] = Cm_minmax
parameters[[7]] = Cbrown_minmax

#: create vectors of possible parameter values with uniform distribution and replacement
################################################
N_sim = LUTsize #: number of parameter combinations (1000 for testing, 100,000 for proper analysis, Duan et al. 2011 / Atzberger 2012)
create_val = function(variable){
  variablevec = seq(from = variable[1], to = variable[2], length.out = variable[3])
  variablevec = sample(variablevec, size = N_sim, replace=TRUE)
  return(variablevec)
}

#: fill LUT with sampled parameter values
################################################
LUT = matrix(NA, nrow=N_sim, ncol=length(parameters))
for (kk in 1:length(parameters)){
  LUT[,kk] = create_val(parameters[[kk]])
}
colnames(LUT) = c("N", "Cab", "Car", "Anth", "Cw", "Cm", "brown")

# Couple values of LMA and EWT to enhace LMA retrivals; only under well watered conditions
# Weiss, M., Baret, F., Myneni, R., Pragnere, A. & Knyazikhin, Y. 2000. Investigation of a model inversion technique to estimate
#        canopy biophysical variables from spectral and direc- tional reflectance data.Agronomie 20: 3–22.
# Combal, B., Baret, F., Weiss, M., Trubuil, A., Mace, D., Pragnere, A., Mynenic, R., Knyazikhin, Y. & Wang, L. 2003.
#        Retrieval of canopy biophysical variables from bidirectional reflectance: using prior information to solve the ill-posed inverse problem.
#        Remote Sensing of Environment 84: 1–15.
for (i in 1:nrow(LUT)){
  LUT[i,6] = LUT[i,5]/sample(seq(2.5,5,length.out=1000), 1, replace=TRUE)
}

#: simulate LUT spectra
################################################
LUT_spectra = matrix(NA,ncol=2101, nrow=LUTsize)
for (iiii in 1:N_sim){
  LUT_spectra[iiii,] = prospect_DB(N=LUT[iiii,1],Cab=LUT[iiii,2],Car=LUT[iiii,3], Ant=LUT[iiii,4],Cw=LUT[iiii,5],Cm=LUT[iiii,6], Brown=LUT[iiii,7])
  print(paste("LUT entry", iiii, "of", LUTsize))
  flush.console()
}
LUT_spectra=LUT_spectra[,(minband:maxband)-399]

if(wavelet_yesno==1){
  #calculate wavelets for LUT
  wavelet_LUT3 = matrix(NA, nrow=nrow(LUT_spectra), ncol=ncol(LUT_spectra))
  wavelet_LUT4 = matrix(NA, nrow=nrow(LUT_spectra), ncol=ncol(LUT_spectra))
  wavelet_LUT5 = matrix(NA, nrow=nrow(LUT_spectra), ncol=ncol(LUT_spectra))
  wavelet_LUT6 = matrix(NA, nrow=nrow(LUT_spectra), ncol=ncol(LUT_spectra))

  for (i in 1:nrow(LUT_spectra)){
    wavelet_LUT = as.matrix(wavCWT(x=as.numeric(LUT_spectra[i,]), n.scale=wavelet_no, scale.range=wavelet_range))
    wavelet_LUT3[i,] = wavelet_LUT[,3]
    wavelet_LUT4[i,] = wavelet_LUT[,4]
    wavelet_LUT5[i,] = wavelet_LUT[,5]
    wavelet_LUT6[i,] = wavelet_LUT[,6]
    print(paste("LUT wavelet", i, "of", LUTsize))
    flush.console()
  }
}

#: inversion
##################################
N_sim = LUTsize
LUTcost_vecorder = 1:N_sim


results_all = list()
cl<-makeCluster(no_cores)
registerDoParallel(cl)

# run inversion in parallel processing
results_all <- foreach(i = 1:nrow(spectra)) %dopar% {

  results_run=list()
  #for (i in 1:nrow(spectra)){

  if(wavelet_yesno==1){
    wavelet_spectra = as.matrix(wmtsa::wavCWT(x=as.numeric(spectra[i,]), n.scale=wavelet_no, scale.range=wavelet_range))
    LUTcost_vec3 = sqrt(1/length(wavelet_spectra[,3])*rowSums(sweep(-wavelet_LUT3, 2, -wavelet_spectra[,3])^2))
    LUTcost_vec4 = sqrt(1/length(wavelet_spectra[,4])*rowSums(sweep(-wavelet_LUT4, 2, -wavelet_spectra[,4])^2))
    LUTcost_vec5 = sqrt(1/length(wavelet_spectra[,5])*rowSums(sweep(-wavelet_LUT5, 2, -wavelet_spectra[,5])^2))
    LUTcost_vec6 = sqrt(1/length(wavelet_spectra[,6])*rowSums(sweep(-wavelet_LUT6, 2, -wavelet_spectra[,6])^2))

    #LUTcost_vec = LUTcost_vec3 + LUTcost_vec4 + LUTcost_vec5 + LUTcost_vec6
    LUTcost_vec = LUTcost_vec4 + LUTcost_vec5 + LUTcost_vec6
  }else{
    LUTcost_vec = sqrt(1/nrow(LUT_spectra)*rowSums(sweep(-LUT_spectra, 2, -as.numeric(spectra[i,]))^2))
  }


  LUTcost_matselect = LUTcost_vecorder[order(LUTcost_vec)]
  LUTcost_matselect = LUTcost_matselect[1:(N_sim/100*percent)]

  ordered_RMSE = LUTcost_vec[LUTcost_matselect]
  weight = (1/ordered_RMSE)/sum(1/ordered_RMSE)

  #: return estimated parameters
  ################################################
  #weighted parameter = sum(weight * selected variables); weight = (1/RMSE)/sum(1/RMSE) (see Vohland et al. 2010)
  estimates_traits = colSums(apply(LUT[LUTcost_matselect,], MARGIN=2, function(x) x*weight))
  estimates_sd = apply(LUT[LUTcost_matselect,], 2, sd)
  estimates_spectra = colSums(LUT_spectra[LUTcost_matselect,]*weight)
  estimates_meanRMSE = sqrt((mean(as.numeric(spectra[i,]) - estimates_spectra)^2))#estimates_spectra[i,]

  png(file = paste("inverted_spectra/inversion_PD_",i,".png",sep=""),width = 400,units="px", height = 400, res = 80,bg = "white")
  plot(minband:maxband,as.numeric(spectra[i,]), xlab="wavelength [nm]", ylab="reflectance [%]", col="red", type="l", ylim=c(0,1), main=paste("pot: ",rownames(spectra)[i],"  RMSE: ",round(estimates_meanRMSE[i],3), sep=""))
  lines(minband:maxband, estimates_spectra, col="black")#estimates_spectra[i,]
  grid()
  dev.off()

  results_run[[1]] = estimates_traits
  results_run[[2]] = estimates_sd
  results_run[[3]] = estimates_spectra
  results_run[[4]] = estimates_meanRMSE

  return(results_run)
}

# bind all results
estimates_traits = do.call(rbind, estimates_traits)

colnames(estimates_traits) = c("N", "Cab", "Car", "Anth", "CW", "Cm", "Cbrown")
rownames(estimates_traits) = rownames(spectra)

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
