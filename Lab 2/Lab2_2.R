setwd("C:/Users/ajs0428/OneDrive - UW/University of Washington/Teaching/SEFS433_Lidar/Labs/Lab 2/")

install.packages("lidR")

library(lidR)
library(rgl)


LASfile <- ("WA_031_rp.laz")
las <- readLAS(LASfile)

plot(las)
rgl::close3d()

plot(las, color = "RGB")
?lidR::clip
