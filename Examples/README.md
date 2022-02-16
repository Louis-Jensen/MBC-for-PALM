# Multiple-blinking correction

This is the R code for the paper:
**Correction of multiple-blinking artefacts in photoactivated localisation microscopy**

The code outputs a set of corrected localisations, obtained after estimating blinking dynamics, and clustering input localisations accordingly.

We provide two examples of fluorophore proteins on a grid, with either light or heavy blinking (Data_light.csv, Data_heavy.csv). 
Additionally, we provide an example of experimental reference data on nuclear pore complexes (Data_NPC.csv).

## Coding environment
We assume that RStudio and R are installed, as well as the following R packages:
- tidyverse
- magrittr
- wrapr
- spatstat
- multiway
- fastcluster
- dbscan

Extract the `mbc.zip` archive and set the RStudio working directory as the extracted folder (where main.R is).

In RStudio, running `source("example.R")` will perform the multiple-blinking correction for the light blinking data, and output to Data_corrected.csv, a set of corrected fluorophore coordinates. Change the file name in example.R from Data_light.csv to Data_heavy.csv to see the corrected output from a heavy blinking input.

## Data input
The code takes as input data in a Comma Separated Value file (.csv) with the fields (x,y,time,sd) representing x, y, frame number and localisation uncertainty (standard deviation) respectively.

## User-defined parameters 
Set the camera framerate, default 25.

`framerate <- 25`

Expected proportion of non-noise points, default 1 (corresponding to no noise)

`alpha <- 1`

## Output
A .csv file containing the fields (x,y,sd,n) representing corrected x, y coordinates, updated localisation uncertainty (standard deviation) and number of input points merged to one, corrected localisation.
