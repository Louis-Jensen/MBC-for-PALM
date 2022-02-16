# Multiple-blinking correction

This is the R code for the paper:
**Correction of multiple-blinking artefacts in photoactivated localisation microscopy**
The MBC-code outputs a set of corrected localisations, obtained after estimating blinking dynamics, and clustering input localisations accordingly.
Additionally, code for simulation of blinking clusters is also included.

The example files 'example_LAT.R' and 'example_Simulation.R' are a gentle introduction to both simulation and MBC on real and simulated data. The interested reader can find additional information on the code below.

## Coding environment
We assume that RStudio and R are installed, as well as the following R packages:
- tidyverse
- magrittr
- wrapr
- spatstat
- multiway
- fastcluster
- dbscan
- markovchain

For the examples to work, first set the RStudio working directory to the top-most folder of the downloaded repo (where e.g. the 'Code' folder is).
After setting the working directory, run the command "source('./Code/Dependencies.R')" to load the necessary functions for MBC and data simulation.

# The MBC() function
The MBC() function takes as input data in a dataframe with the fields (x,y,time,sd) (in that order) representing x, y, frame number and localisation uncertainty (standard deviation) respectively. 

## User-defined parameters 
MBC() expects the framerate used to record the data, e.g. 25.

Other parameters that can be specified (but have reasonable defaults):
- alpha   : Fraction of localisations that are not background noise (see the paper for details on how this can be estimated).
- t.thres : Grid of values to examine when looking for the best time-dilation, S. If the best value is found at either end of the interval (check plot), consider widening the default interval and rerunning the clustering.

## Output
A dataframe containing the fields (x,y,uncertainty) representing corrected x, y coordinates, and updated localisation uncertainty (standard deviation).

# The Add_blinking() function
The Add_blinking() function takes as input a dataframe with the fields (x,y) (in that order) representing the real x, y coordinates of 'proteins', i.e. before any blinking artefacts are added. 

## User-defined parameters
The user can specify the parameters that govern data artefacts:
- parameters: the blinking rates, as a vector of 4 values corresponding to blinking rates in the order (r_F, r_D, r_R, r_B).
- framerate: the framerate of the camera that 'records' the PALM data.
- stds: empirical distribution of localisation uncertainties. The localisation uncertainties in the simulation are drawn with replacement from stds. Has a reasonable default.
- stop_if_unusual: when True, the function checks for parameters that result in very unrealistic (and presumably unintended) artefact behavior. Set to False if you wish to explore extreme circumstances.

## Output
A dataframe containing the fiels (x, y, time, sds) representing the x,y, coordinates of blinks, the camera frame on which each point was observed, and the localisation uncertainty.
