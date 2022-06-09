# ********************************************************************************
# main routine to generate leaf optical properties using PROSPECT-D
# ********************************************************************************
# _______________________________________________________________________
# for any question or request, please contact: 
#
# Jean-Baptiste FERET
# UMR-TETIS, IRSTEA Montpellier
# Maison de la Télédétection
# 500 rue Jean-Fracois Breton
# 34093 Montpellier cedex 5
# E-mail: jb.feret@teledetection.fr
#
# Stéphane JACQUEMOUD
# Université Paris Diderot / Institut de Physique du Globe de Paris
# 35 rue Hélène Brion
# 75013 Paris, France
# E-mail: jacquemoud@ipgp.fr
#
# http://teledetection.ipgp.fr/prosail/
# _______________________________________________________________________
# ********************************************************************************
# Main script running with PROSPECT-DB version 6.0 (16 January 2017)
# ********************************************************************************


source("scripts/R/PROSPECT_D/dataSpec_PDB.r")
source("scripts/R/PROSPECT_D/calctav.r")
source("scripts/R/PROSPECT_D/prospect_DB.r")
require(expint)

# load example biochemistry and structure variables
N      = 1.5		# structure coefficient
Cab    = 20		# chlorophyll content (µg.cm-2) 
Car    = 10		# carotenoid content (µg.cm-2)
Ant   = 0.5		# Anthocyanin content (µg.cm-2)
Brown = 0.1		# brown pigment content (arbitrary units)
Cw     = 0.015	# EWT (cm)
Cm     = 0.009	# LMA (g.cm-2)



# load('leaf_parameter.txt')
# N     = leaf_parameter[1]
# Cab   = leaf_parameter[3]
# Car   = leaf_parameter[2]
# Anth  = leaf_parameter[4]
# Cbrown= leaf_parameter[5]
# Cw    = leaf_parameter[6]
# Cm    = leaf_parameter[7]

refl=prospect_DB(N,Cab,Car,Ant,Brown,Cw,Cm)
plot(refl, type="l", xlab="wavenlength [nm]", ylab="reflectance")
