#Setting the working directory
setwd("/Users/ajs0428/OneDrive - UW/University of Washington/Teaching/SEFS433_Lidar/Labs")

#Installing lidR and loading the library
install.packages("lidR")
library(lidR)

#Installing rgl and loading the library
install.packages("rgl")
library(rgl)

#creating a file path string for the .laz lidar file
LASfile <- ("Lab 2/WA_031_rp.laz")
#reading in the the LASfile with the readLAS function
las <- readLAS(LASfile)
#plotting the las file
plot(las)
#printing the las file properties
print(las)
#printing more details
print(las@header)

las_check(las)

plot(las, color = "Intensity")
plot(las, color = "Intensity", breaks = "quantile", nbreaks = 10)
plot(las, color = "RGB")


plot_print <- function(x){
  plot(x)
  print(x)
  las_check(x)
}

plot_print(las)
