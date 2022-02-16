# Make sure to set the working directory to the location that contains the Code folder.
# setwd(...)
source("./Code/Dependencies.R")
set.seed(1)

# We load a LAT region
data = read.csv("./Experimental data/Linker for Activation of T cells/WT/Central/Data_WTCentral1_1.csv")
data = data[c('x', 'y', 'frame', 'uncertainty_xy_nm')]

# The algorithm estimates rate parameters and performs clustering to data. 
# Expect it to take 2-4 minutes.
corrected = MBC(data, framerate = 25, plot = T)
