setwd("//fs-persona.sefs.uw.edu/student_redirect$/JonBatch/Desktop/ESRM433/")

install.packages("lidR") #downloads and installs lidR and dependent packages
install.packages("rgl") #downloads and installs
install.packages("gstat") #optional, only needed for lasground
install.packages("sf") #needed for readLAScatalog. DO NOT Sellect "yes" for binary
install.packages("mapview") #needed to view tiles with world imagry
install.packages("BiocManager")#needed to install EBImage
BiocManager::install("EBImage") #needed for 
install.packages("rgdal") #required for mapping, should already be installed.
install.packages("concaveman")#optional, needed to make detailed, non-overlaping hulls 

library(lidR) #loads the lidR package into RStuio 
library(rgl) #loads package into RStudio 
library(gstat) #optional
library(sf) #needed for readLAScatalog
library(mapview) 
library(BiocManager) #loading BiocManager
library(EBImage)#loading EBImage
library(rgdal) #run if you get an error "rgdal required"
library(concaveman) #optional

LASfile <- ("LAB3/q47123g8201.laz") #tells R that 'LASfile' is the file q47123g8201.laz
las <- readLAS(LASfile) #uses function 'readLAS' on 'LASfile' 
plot(las)  #plot the point cloud 'las' with a color height ramp applied

print(las)  
print(las@header)  # Basic information such as file size, extent, projection, and more

lascheck(las)

plot(las, color="Intensity") #different methods of visualizing the point cloud
plot(las, color="Intensity", trim=100) 
plot(las, color="RGB")

############################ LAB3 ############################

help(package = lidR) #commands to access the help documentation for each package 
help(package = raster)
help(package = sp)

crs <- sp::CRS("+init=epsg:2927") #associates 'crs' as representing projection with the EPSG code 2927

?projection #help documentation for the 'projection' function
projection(las) #checks the projection of the 'las' object
projection(las) <- crs #changes the projection of the 'las' object to EPSG:2927 

?grid_terrain #help documentation
dtm1 = grid_terrain(las, res = 3.28, algorithm = knnidw(k = 6L, p = 2)) #uses nearest neighbor to create DTM 
dtm2 = grid_terrain(las, res = 3.28, algorithm = tin()) #uses triangulated irregular network to creat DTM
dtm3 = grid_terrain(las, res = 3.28, algorithm = kriging(k = 10L)) #optional, long processing time

plot(dtm1, main='knnidw')
plot(dtm2, main='tin')

?plot_dtm3d
plot_dtm3d(dtm1)
plot_dtm3d(dtm2)

?lasclip
lasclip <- lasclipCircle(las,797000,944000,100) #creates a circle clip of 100unit (ft) radius
plot(lasclip) #plots the 3D point cloud in the clip with a height color ramp

?lasground
?pmf
ws  <- seq(3,12, 3)
th  <- seq(0.1, 1.5, length.out = length(ws))
laspmf <- lasground(lasclip, pmf(ws, th))
plot(laspmf, color = "Classification")

?csf
mycsf <- csf(FALSE, 1, 1, time_step = 1)
lascsf <- lasground(lasclip, mycsf)
plot(lascsf, color = "Classification")

dtmPMF = grid_terrain(laspmf, res = 3.28, algorithm = tin())
plot_dtm3d(dtmPMF)

dtmcsf = grid_terrain(lascsf, res = 3.28, algorithm = tin())
plot_dtm3d(dtmcsf)

?writeLAS #help info on 'writeLAS'
writeLAS(laspmf, file = "LAB3/laspmf.laz") #the name and location to write the las (or laz) file

?readLAScatalog #help info on 'readLAScatalog'
ctg <- readLAScatalog("LAB3/TILES/") # looks at the folder TILES and loads all las files

lascheck(ctg) #runs basic check on all las files in TILES folder
plot(ctg) #plots out the las files geographically

TILEdtm = grid_terrain(ctg, res = 3.28, algorithm = tin()) #creates a singel, 1m (3.28ft) resolution DTM from all tiles

plot(TILEdtm, main='TILE') #produces height map of the TILEdtm
plot_dtm3d(TILEdtm) #3D model of our output DTM

projection(TILEdtm) #checks the projection information on our output DTM
projection(TILEdtm) <-crs #corrects the projection using predefined 'crs' 

?writeRaster #help info on 'writeRaster' function
writeRaster(TILEdtm,filename =file.path("LAB3/TILEdtm.tif"),format="GTiff", overwrite=TRUE) #location and file type for output

############################ LAB4 ############################


ctg <- readLAScatalog("LAB3/TILES/") # looks at the folder TILES and loads all las files

?mapview
plot(ctg, mapview = TRUE, map.types = "Esri.WorldImagery") #places tiles onto satalite imagry
#This is awesome to check to make sure the projection is correct, and it is not!

NEWcrs <- sp::CRS("+init=epsg:2927") #setting the objet "NEWcrs" as the crs EPSG:2927

crs(ctg) #checking what the crs is for the las files in the catalog "ctg"

crs (ctg)<- NEWcrs #changing the defined projection informaton to EPSG 2927 

plot(ctg) #quick way to plot the position of the tiles, we want a clip

lasclip <- lasclipCircle(ctg,798000,944000,200) #creates a circle clip of 200unit (ft) radius
plot(lasclip) 

chmHILL = grid_canopy(lasclip, 3.28, p2r()) 
#Creates a DSM of the lasclip with a pixel resolution of 1m
plot(chmHILL, col = height.colors(50)) #plots the DSM using a defined height color ramp
plot_dtm3d(chmHILL)

?lasnormalize
NORMclip <- lasnormalize(lasclip,tin()) #runs lasnormalize on "lasclip" using the "tin" method. 
plot(NORMclip)
chmNORM = grid_canopy(NORMclip, 3.28, p2r()) 
#Creates a DSM of the lasclip with a pixel resolution of 1m
plot(chmNORM, col = height.colors(50)) #plots the DSM using a defined height color ramp
plot_dtm3d(chmNORM)

?grid_canopy
chm1 = grid_canopy(NORMclip, 1, p2r())
plot(chm1, col = height.colors(50), main = "Res1, p2r0")
plot_dtm3d(chm1)

chm2 = grid_canopy(NORMclip, 3.28, p2r())
plot(chm2, col = height.colors(50), main = "Res3.28, p2r0")
plot_dtm3d(chm2)

chm3 = grid_canopy(NORMclip, 1, p2r(3))
plot(chm3, col = height.colors(50), main = "Res1, p2r3")
plot_dtm3d(chm3)

chm4 = grid_canopy(NORMclip, 1, dsmtin())
plot(chm4, col = height.colors(50), main = "dsmtin")
plot_dtm3d(chm4)

chm5 = grid_canopy(NORMclip, 1, pitfree(c(0,2,5,10,15), c(0,1), subcircle = 2))
plot(chm5, col = height.colors(50), main = "pitfree")
plot_dtm3d(chm5)


# smoothing post-process 
?focal
chmSmooth = focal(chm4, w = matrix(1,5,5), fun = mean)
plot(chmSmooth, col = height.colors(50), main = "Smoothed")
plot_dtm3d(chmSmooth)

help(package = EBImage)

?lidR::watershed
?lastrees
algo = lidR::watershed(chm5, th = 20) #defines the algorythm to be used with lastrees
algoS = lidR::watershed(chmSmooth, th = 20) #overwriteing to use chmSmooth
lasT  = lastrees(NORMclip, algo) #adds a new scalar field to the las data "TreeID"
lasTS  = lastrees(NORMclip, algoS)
#this requires EBImage to be installed

trees = lasfilter(lasT, !is.na(treeID))#remove points that are not assigned to a tree
treesS = lasfilter(lasTS, !is.na(treeID))

plot(trees, color = "treeID", colorPalette = pastel.colors(100))
plot(treesS, color = "treeID", colorPalette = pastel.colors(100))

writeLAS(trees, file = "LAB4/trees.las") #you can write your trees las file.

?tree_hulls
hulls  = tree_hulls(lasT, func = ~list(Zmax = max(Z)))
hullsS  = tree_hulls(lasTS, func = ~list(Zmax = max(Z)))
#asigns each hull a value equal to the maximum Z value within it (i.e. tree height)
spplot(hulls,"Zmax") #super cool plot of the hulls colored by height
spplot(hullsS,"Zmax")
print(hulls) #information about extent, projeciton, and number of polygons
print(hullsS)

crowns = lidR::watershed(chm5, th = 20)()
crownsS = lidR::watershed(chmSmooth, th = 20)()
#creates a raster of the treesegmentation, th is th minimum size a "tree" must be
crownsS = lidR::watershed(chmSmooth, th = 20)()
plot(crowns, col = pastel.colors(100)) #plots the segmentation with unique tree IDs
plot(crownsS, col = pastel.colors(100))

contour = rasterToPolygons(crowns, dissolve = TRUE)

plot(chm5, col = height.colors(50))
plot(contour, add = T)

?mapview
mapview(chm5)
mapview(contour)
mapview(hulls,"Zmax")

### Arboretum Code ###
ArborLASfile <- ("LAB4/q47122f3417_rp.laz")
Arborlas <- readLAS(ArborLASfile)
print(Arborlas)
lascheck(Arborlas)

ArborClip <- lasclipRectangle(Arborlas, 1197000, 845300, 1198200, 846700)
plot(ArborClip)

ArborNORM <- lasnormalize(ArborClip,tin()) 
crs(ArborNORM)
NEWcrs <- sp::CRS("+init=epsg:2927")
crs(ArborNORM) <-NEWcrs

ArborCHM = grid_canopy(ArborNORM, 1.64, p2r(2))
plot(ArborCHM)
mapview(ArborCHM)
plot_dtm3d(ArborCHM)

ArborSmooth = focal(ArborCHM, w = matrix(1,3,3), fun = max)
plot(ArborSmooth, col = height.colors(50), main = "ArborSmoothed")
plot_dtm3d(ArborSmooth)

algo = lidR::watershed(ArborSmooth, th = 10) 
ArborT  = lastrees(ArborNORM, algo)

ArborTrees = lasfilter(ArborT, !is.na(treeID))

plot(ArborTrees, color = "treeID", colorPalette = pastel.colors(100))
writeLAS(ArborTrees, file = "LAB4/Arbortrees.las")
ArborHulls  = tree_hulls(ArborTrees, func = ~list(Zmax = max(Z)))
#ArborHC  = tree_hulls(ArborTrees, type = "concave", concavity = 2, func = ~list(Zmax = max(Z)))
#Above will output more detailed polygons for ArborHulls, but requires the package "concaveman"
#concaveman takes much longer to run on all steps due to the much higher detail in the polygons
#but you are incouraged to try it, not required.

spplot(ArborHulls,"Zmax") 
print(ArborHulls)
mapview(ArborHulls,"Zmax")

#create a shp file to import into a GIS program (ArcGIS Pro)
shapefile(ArborHulls,filename='LAB4/ArborHulls.shp',overwrite=TRUE)

##############################################################################
# For the code above, all that is actually needed to go from a downloaded laz tile to a
# map of tree crowns and heights is the code below which you could sellect all the lines
# and run them all at once.
ArborLASfile <- ("LAB4/q47122f3417_rp.laz")
Arborlas <- readLAS(ArborLASfile)
ArborClip <- lasclipRectangle(Arborlas, 1197000, 845300, 1198200, 846700)
ArborNORM <- lasnormalize(ArborClip,tin()) 
NEWcrs <- sp::CRS("+init=epsg:2927")
crs(ArborNORM) <-NEWcrs
ArborCHM = grid_canopy(ArborNORM, 1.64, p2r(2))
ArborSmooth = focal(ArborCHM, w = matrix(1,3,3), fun = max)
algo = lidR::watershed(ArborSmooth, th = 10) 
ArborT  = lastrees(ArborNORM, algo)
ArborTrees = lasfilter(ArborT, !is.na(treeID))
ArborHulls  = tree_hulls(ArborTrees, func = ~list(Zmax = max(Z)))
shapefile(ArborHulls,filename='LAB4/ArborHulls.shp',overwrite=TRUE) 
##############################################################################

############################ LAB5 ############################

LASfile <- ("LAB5/Sequoia.las")
las <- readLAS(LASfile)
lasNORM <- lasnormalize(las,tin())
writeLAS(lasNORM, file = "LAB5/SequoiaNorm.las")

?cloud_metrics
cloud_metrics(lasNORM, ~max(Z))

CMSequoia <- cloud_metrics(lasNORM, .stdmetrics)
write.csv(CMSequoia, file = "LAB5/CMSequoia.csv")

over2 = lasfilter(lasNORM, Z > 6.56)
CMSequoia2 <- cloud_metrics(over2, .stdmetrics)
write.csv(CMSequoia2, file = "LAB5/CMSequoia2.csv")

field <-(read.csv("LAB5/plotdata.csv"))
lidar <- (read.csv("LAB5/cloudmetrics.csv"))

?lm
Y <- field$QMD.GT5 #Dependent Variable.Identifies QMD.GT5 column from the field dataset
X <- lidar$Elev.P95 #Independent Variable. Identifes Elev.P95 column from the lidar dataset 

reg <- lm (Y ~ X) # reg is the regression output from Y & X
sum <- summary(reg) # sum is the summary statistics of the linear regression of Y & X 
sum # prints the summary statistics. p-value and R-squared value are of the most interest to us

coeff <- coefficients(reg) #coefficients are the slope of the regression line
#eq will be the title of the plot with slope and R-Squared values 
eq <- paste0("Slope=",round(coeff[2],1),"*x+",round(coeff[1],1)," R^2=",round(sum$r.squared, 2))
plot(X, Y, xlab="P95", ylab="QMD", main=eq) #plots the points and title
abline(reg, col="blue") #adds the regression line to the plot

Mreg <-lm(field$Carbon.AB ~ lidar$Elev.stddev + lidar$Elev.mean + 
          lidar$Percentage.all.returns.above.6.56 + lidar$Elev.P95 + lidar$Elev.P25)
summary (Mreg)
plot(Mreg)

LASfile <- ("LAB5/CloudSection.las") 
las <- readLAS(LASfile) 
plot(las)

#create 5 circle clips and then normalize the clips
P1 <- lasclipCircle(las,1636,-2636,60)
P1n <- lasnormalize(P1,tin())
plot(P1n)
P2 <- lasclipCircle(las,1430,-2230,60)
P2n <- lasnormalize(P2,tin())
plot(P2n)
P3<- lasclipCircle(las,1216,-2425,60)
P3n <- lasnormalize(P3,tin())
plot(P3n)
P4 <- lasclipCircle(las,1279,-2725,60)
P4n <- lasnormalize(P4,tin())
plot(P4n)
P5 <- lasclipCircle(las,1139,-2174,60)
P5n <- lasnormalize(P5,tin())
plot(P5n)

#Combines all the non-normalized clips
#and creates a single las file from the clips
PLOTSlas <- rbind(P1,P2,P3,P4,P5)
plot(PLOTSlas)
writeLAS(PLOTSlas, file = "LAB5/PLOTSlas.laz")

#Creates standard cloud metrics from the normalized clips
#only uses the points above 2m (6.56ft)
?cloud_metrics
P1m <- cloud_metrics((lasfilter(P1n,Z > 6.56)), .stdmetrics)
P2m <- cloud_metrics((lasfilter(P2n,Z > 6.56)), .stdmetrics)
P3m <- cloud_metrics((lasfilter(P3n,Z > 6.56)), .stdmetrics)
P4m <- cloud_metrics((lasfilter(P4n,Z > 6.56)), .stdmetrics)
P5m <- cloud_metrics((lasfilter(P5n,Z > 6.56)), .stdmetrics)

#combines all of the clip cloud metrics into a single data frame
#rbind is merging (bind) data by stacking rows (r)
CM <- rbind.data.frame(P1m,P2m,P3m,P4m,P5m)

#creates a metric for canopy cover by counting the number of points 
#above 2m (sum(Z>6.56)) and deviding by the number of all points (sum(Z>-1))
#number of points is >-1 as a few poitns may have a slightly negitive value
#due to the tin used to create the normalized point cloud.
P1c <- cloud_metrics(P1n, ~sum(Z>6.56)/sum(Z>-1))
P2c <- cloud_metrics(P2n, ~sum(Z>6.56)/sum(Z>-1))
P3c <- cloud_metrics(P3n, ~sum(Z>6.56)/sum(Z>-1))
P4c <- cloud_metrics(P4n, ~sum(Z>6.56)/sum(Z>-1))
P5c <- cloud_metrics(P5n, ~sum(Z>6.56)/sum(Z>-1))

#combineds our canopy clousure values into one data frame
CC <- rbind(P1c,P2c,P3c,P4c,P5c)

#reads in the csv file of our plot data and creates a data frame
plot2 <- read.csv("LAB5/plot2.csv")

#combines our outputs togeter for one datafrom. 
#bindes the columns instead of the rows
D <-cbind.data.frame(plot2,CC,CM)

Y <- D$BASAL #Dependent Variable. Basal area from our field data
X <- D$zq95 #Independent Variable. Height where 95% of all points are below

reg <- lm (Y ~ X) # reg is the regression output from Y & X
sum <- summary(reg) # sum is the summary statistics of the linear regression of Y & X 
sum # prints the summary statistics. p-value and R-squared value are of the most interest to us

coeff <- coefficients(reg) #coefficients are the slope of the regression line
#eq will be the title of the plot with slope and R-Squared values 
eq <- paste0("Slope=",round(coeff[2],1),"*x+",round(coeff[1],1)," R^2=",round(sum$r.squared, 2))
plot(X, Y, xlab="zq95", ylab="Basal", main=eq) #plots the points and title
abline(reg, col="blue") #adds the regression line to the plot

Freg <-lm( D$QMD ~ D$zmean + D$zq95 + D$CC + D$zq25 + D$zsd)#don't run all these at once
summary (Freg)

############################ LAB6 ############################



