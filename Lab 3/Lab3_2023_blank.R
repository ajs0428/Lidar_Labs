
#### Before starting the lab here are some tips####
# character or strings are contained in quotes such as 
"string" 
# or 
'character'

# variables don't use any formatting and have no spaces
# They can be set with the <- or =
    # R syntax commonly use <-
a_string <- "string"
print(a_string) # 
# You can use the 'print' function to print strings/characters or the output
# of other variables

# Use ? to look up functions such as:
?readLAS #look at the help window that should open

# Also functions take arguments within parentheses ()
  # such as the print function that takes 1 argument to print but can take others
?print
## S3 method for class 'factor'
# print(x, quote = FALSE, max.levels = NULL,
#       width = getOption("width"), ...)



#### Start of the lab 3 code ####

#Set the working directory which is the location of all working files
setwd("")

#install packages if needed
#install.packages("lidR")

#loading libraries needed after installing
library(lidR)
library(rgl)
library(terra)

#set the variable 'LASfile' to the path of the file 
LASfile <- ("")
#read the 'LASfile' with the lidR readLAS function that opens the point clouds
las <- readLAS(________)

#check the las point cloud with a plot
plot(______)

#print the las file and examine output
print(______)
#the projection function prints a long string that defines the coordinate projection system
    # look at the familiar numbers
projection(___________)

#Set the coordinate reference system (crs) variable to the EPSG code 2927
crs <- st_crs(___________)
#Then set the projection(las) to the new coordinate reference system variable
projection(las) <- _________
  
#check out the projection now  
  # is there a Washington South and (ftUS)?
projection(las)

#Check out other plot variations herer
# plot(las, color = "Intensity")
# plot(las, color = 'Intensity', breaks = "quantile", nbreaks = 10)
# plot(las, color = "RGB")

#Check out the rasterize_terrain function
?rasterize_terrain

#create a dtm based on the triangulation 
dtm_tin <- rasterize_terrain(_________________________________)
#check out the result of your dtm in 3D
plot_dtm3d(dtm_tin)
#now check it out in 2D
plot(dtm_tin)

#create a dtm based on k-nearest neighbors inverse distance weighting
dtm_knnidw <- rasterize_terrain(_________________________________)
#similarly check out in 3D or 2D
plot_dtm3d(dtm_knnidw) 
plot(dtm_knnidw)

#print out las again to find the extent of the lidar data
print(las)

#now with that extent in mind use the clip_circle function to clip out a small
  #section of the las file
las_clip <- clip_circle(_________________________________)
#check out your section
plot(las_clip)

# These are the arguments that go into the classify ground function
ws <- seq(3, 12, 3) #This makes a sequence of numbers from 3 to 12 by 3
# this also makes a sequence of numbers from 0.1 to 1.5 but evenly 
th <- seq(0.1, 1.5, length.out = length(ws)) # it also sets the length of the sequence to the same as ws
# Now run for the pmf algorithm
laspmf <- classify_ground(_________, pmf(_________)) #replace the _______ here 
plot(laspmf, color = "Classification")

# And run for the csf algorithm
mycsf <- csf(TRUE, 1, 1, time_step = 1)
lascsf <- classify_ground (_________, _________) #replace the ____________ here

#Create a DTM from your clips by running:
  #use your laspmf
dtmPMF = rasterize_terrain(_______, res=3.28, algorithm = tin())
plot_dtm3d(dtmPMF)
  #use your lascsf
dtmCSF = rasterize_terrain(________, res=3.28, algorithm = tin())
plot_dtm3d(dtmCSF)

#We can write the laz file out to our computer if we want to bring it d
  # to cloudcompare
writeLAS(laspmf, file = "Labs/Lab 3/laspmf_test.laz")


#Bring in the sf library for using with readLAScatalog
library(sf)
#change the location of files for the readLAScatalog function
ctg<- readLAScatalog("")
las_check(ctg)
#plot the catalog
plot(ctg)

#use the same rasterize_terrain as earlier for the ctg catalog variable, 
  #but it can recognize 
  #that this is a collection of tiles and we can make it a variable called 
    #TileDTM
TileDTM <- rasterize_terrain()
#plot 2D
plot(_________, main = "Tile")
#plot 3D
plot_dtm3d(________)

#check the projection of TileDTM
projection(TileDTM)

#lets fix the projection with the new EPSG code using the terra package

#install.packages("terra") #just in case it's not installed
library(terra) #just making sure it's loaded
terra::crs(TileDTM) <- ""
#check the projection again
projection(TileDTM)

#write it to a new raster file on the computer using terra
terra::writeRaster(_________, filename = "", overwrite=TRUE)


