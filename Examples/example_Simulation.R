# Make sure to set the working directory to the location that contains the Code folder.
# setwd(...)
source("./Code/Dependencies.R")
set.seed(1)

# The simulation function takes as input (X,Y)-coordinates of real 'protein' locations.
# Here, we simulate some Poisson proteins.
proteins = rpoispp(0.001, win = square(1000))
proteins_df = data.frame(X = proteins$x, Y = proteins$y)

par(mfrow = c(1,3))
plot(proteins_df, asp = 1, pch = 19, main = "Proteins")

# We add blinking following the 4-state blinking model, with specified rates and camera framerate.
# Additionally, we specify an empirical distribution of localisation uncertainties (has a reasonable default).
blinking_params = c(0.005, 5, 0.75, 2.5) # (r_F, r_D, r_R, r_B)
framerate = 25
uncertainties = rgamma(1000, 4, 4/20)

blinking = Add_blinking(proteins_df, blinking_params, framerate, uncertainties)

# Plotting xy blinking data
plot(blinking[,1:2], pch = 19, asp = 1, cex = 0.5, main = 'xy-blinking') 

# Plotting frame~x for the first minute of the recording (nicer visuals)
plot(blinking[,c(1,3)], pch = 19, cex = 0.5, ylim = c(0,25*60), main = 'tx-blinking')

# We perform MBC to go back to the true protein locations.
corrected = MBC(blinking, framerate = framerate)

# Comparing
par(mfrow = c(1,1))
plot(proteins_df[,1:2], col = 2, pch = 18, asp = 1, main = 'Truth and corrected')
points(corrected[,1:2], pch = "+")

