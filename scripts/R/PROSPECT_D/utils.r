# ********************************************************************************
# main functions for LUT generation and inversion of PROSPECT-D
# ********************************************************************************
# _______________________________________________________________________
#
# Author: Javier Lopatin
# javier.lopatin@uai.cl
# Universidad Adolfo Ibanez
# Santiago
# Chile
# _______________________________________________________________________
# ********************************************************************************
# Main script running with PROSPECT-DB version 6.0 (16 January 2017)
# ********************************************************************************

require(wmtsa)

########################################################
# Function to get simualted LUT
########################################################
getProsailLUT = function(parameters, LUTsize, wavelet_yesno=1, wavelet_no = 6, wavelet_range = c(1,200)){

    # progress bar
    pb <- txtProgressBar(min = 0,
                         max = LUTsize,
                         style = 3,
                         width = 50,
                         char = "=")

    #: create vectors of possible parameter values with uniform distribution and replacement
    ################################################
    create_val = function(variable){
      variablevec = seq(from = variable[1], to = variable[2], length.out = variable[3])
      variablevec = sample(variablevec, size = LUTsize, replace=TRUE)
      return(variablevec)
    }
    #: fill LUT with sampled parameter values
    ################################################
    LUT = matrix(NA, nrow=LUTsize, ncol=length(parameters))
    for (kk in 1:length(parameters)){
      LUT[,kk] = create_val(parameters[[kk]])
    }
    colnames(LUT) = c("N", "Cab", "Car", "Anth", "Cw", "Cm", "brown")
    # Couple values of LMA and EWT to enhace LMA retrivals; only under well watered conditions
    # Weiss, M., Baret, F., Myneni, R., Pragnere, A. & Knyazikhin, Y. 2000. Investigation of a model inversion technique to estimate
    #        canopy biophysical variables from spectral and directional reflectance data. Agronomie 20: 3–22.
    # Combal, B., Baret, F., Weiss, M., Trubuil, A., Mace, D., Pragnere, A., Mynenic, R., Knyazikhin, Y. & Wang, L. 2003.
    #        Retrieval of canopy biophysical variables from bidirectional reflectance: using prior information to solve the ill-posed inverse problem.
    #        Remote Sensing of Environment 84: 1–15.
    for (i in 1:nrow(LUT)){
      LUT[i,6] = LUT[i,5]/sample(seq(2.5,5,length.out=1000), 1, replace=TRUE)
    }

    #: simulate LUT spectra
    ################################################
    LUT_spectra = matrix(NA,ncol=2101, nrow=LUTsize)
    for (iiii in 1:LUTsize){
      LUT_spectra[iiii,] = prospect_DB(N=LUT[iiii,1],Cab=LUT[iiii,2],Car=LUT[iiii,3],
                                       Ant=LUT[iiii,4],Cw=LUT[iiii,5],Cm=LUT[iiii,6], Brown=LUT[iiii,7])
    }
    LUT_spectra=LUT_spectra[,(minband:maxband)-399]

    if(wavelet_yesno==1){
      #calculate wavelets for LUT
      wavelet_LUT3 = matrix(NA, nrow=nrow(LUT_spectra), ncol=ncol(LUT_spectra))
      wavelet_LUT4 = matrix(NA, nrow=nrow(LUT_spectra), ncol=ncol(LUT_spectra))
      wavelet_LUT5 = matrix(NA, nrow=nrow(LUT_spectra), ncol=ncol(LUT_spectra))
      wavelet_LUT6 = matrix(NA, nrow=nrow(LUT_spectra), ncol=ncol(LUT_spectra))

      for (i in 1:nrow(LUT_spectra)){
        wavelet_LUT = as.matrix(wmtsa::wavCWT(x=as.numeric(LUT_spectra[i,]),
                                              n.scale=wavelet_no, scale.range=wavelet_range))
        wavelet_LUT3[i,] = wavelet_LUT[,3]
        wavelet_LUT4[i,] = wavelet_LUT[,4]
        wavelet_LUT5[i,] = wavelet_LUT[,5]
        wavelet_LUT6[i,] = wavelet_LUT[,6]

        setTxtProgressBar(pb, i) # progress bar
      }
    }

    close(pb) # close progress bar conection

    out = list(wavelet_LUT3, wavelet_LUT4, wavelet_LUT5, wavelet_LUT6)
    return(out)
}


########################################################
# Function to invert PROSPECT-D with wavelets
########################################################
invertProsail = function(parameters, LUTsize, LUT, wavelet_yesno=1, percent = 1,
  wavelet_no = 6, wavelet_range = c(1,200), printPlots = FALSE){

    LUTcost_vecorder = 1:LUTsize

    # LUTs generated
    wavelet_LUT3 = LUT[[1]]
    wavelet_LUT4 = LUT[[2]]
    wavelet_LUT5 = LUT[[3]]
    wavelet_LUT6 = LUT[[4]]

    # cost function
    estimates_traits   = matrix(NA, nrow=nrow(spectra), ncol=length(parameters))
    estimates_sd       =  matrix(NA, nrow=nrow(spectra), ncol=length(parameters))
    estimates_meanRMSE = c()
    estimates_spectra  =  matrix(NA, nrow=nrow(spectra), ncol=ncol(spectra))

    # progress bar
    pb <- txtProgressBar(min = 0,
                         max = LUTsize,
                         style = 3,
                         width = 50,
                         char = "=")

    for (i in 1:nrow(spectra)){

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

      # which observation has less error?
      LUTcost_matselect = LUTcost_vecorder[order(LUTcost_vec)]
      LUTcost_matselect = LUTcost_matselect[1:(LUTsize/100*percent)]

      ordered_RMSE = LUTcost_vec[LUTcost_matselect]
      weight = (1/ordered_RMSE)/sum(1/ordered_RMSE)

      #: return estimated parameters
      ################################################
      #weighted parameter = sum(weight * selected variables); weight = (1/RMSE)/sum(1/RMSE) (see Vohland et al. 2010)
      estimates_traits[i,] = colSums(apply(LUT[LUTcost_matselect,], MARGIN=2, function(x) x*weight))
      estimates_sd[i,] = apply(LUT[LUTcost_matselect,], 2, sd)
      estimates_spectra[i,] = colSums(LUT_spectra[LUTcost_matselect,]*weight)
      estimates_meanRMSE[i] = sqrt((mean(as.numeric(spectra[i,]) - estimates_spectra[i,])^2))

      if (printPlots ==TRUE) {

      png(file = paste("inverted_spectra/inversion_PD_",i,".png",sep=""),width = 400,units="px", height = 400, res = 80,bg = "white")
      plot(minband:maxband,as.numeric(spectra[i,]), xlab="wavelength [nm]", ylab="reflectance [%]",
           col="red", type="l", ylim=c(0,1), main=paste("pot: ",rownames(spectra)[i],"  RMSE: ",round(estimates_meanRMSE[i],3), sep=""))
      lines(minband:maxband,estimates_spectra[i,], col="black")
      grid()
      dev.off()
      }

      setTxtProgressBar(pb, i)
    }

    close(pb) # close progress bar conection

    colnames(estimates_traits) = c("N", "Cab", "Car", "Anth", "CW", "Cm", "Cbrown")
    rownames(estimates_traits) = rownames(spectra)

    return(estimates_traits)
}
