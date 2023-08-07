#First you need to set your working directory and load libraries
setwd("//fs-persona.sefs.uw.edu/student_redirect$/JonBatch/Desktop/ESRM433/")
library(lidR) #loads the lidR package into RStudio

#load las file into R and set the CRS
LASfile <- ("NoClassified.laz") #tells R where your las file is
las <- readLAS(LASfile) #reads in your lidar file into Rstudio
NEWcrs <- sp::CRS("init=epsg:2927") #defines NEWcrs as EPSG 2927
projection(las) <- NEWcrs #defines the CRS for the lidar data

#you must define either pmf or csf to classify the ground points
mycsf <- csf(FALSE, 1, 1, time_step = 1) #defines arguments for csf
lascsf <- classify_ground(las, mycsf) #uses csf to classify ground

#you must define knnidw, tin, or kriging to create the ground model
dtmcsf <- grid_terrain(lascsf, res = 3.28, algorithm = tin())

#checking the output is nessesary
plot_dtm3d(dtmcsf) #quality controll to make sure the dtm is good

#last, write out the dtm so you can use it in GIS programs
writeRaster(dtmcsf,filename =file.path("NEWdtm.tif"),format="GTiff")


