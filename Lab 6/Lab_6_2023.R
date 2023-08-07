install.packages("lidR") #downloads and installs lidR and dependent packages
install.packages("rgl") #downloads and installs graphis package
install.packages("gstat") #optional, only needed for lasground
install.packages("sf") #needed for readLAScatalog. DO NOT Sellect "yes" for binary
install.packages("mapview") #needed to view tiles with world imagry
install.packages("BiocManager")#needed to install EBImage
BiocManager::install("EBImage") #Installs EBImage
install.packages("rgdal") #required for mapping, should already be installed.
install.packages("concaveman")#optional, needed to make detailed, non-overlaping hulls


library(lidR) #loads the lidR package into RStuio 
library(rgl) #loads package into RStudio 
library(gstat) #optional
library(sf) #needed for readLAScatalog
library(mapview) #needed for ploting files on basemaps
library(BiocManager) #loading BiocManager
library(EBImage)#loading EBImage
library(rgdal) #run if you get an error "rgdal required"
library(concaveman) #optional
library(leaflet)

setwd <- "/Users/Anthony/OneDrive - UW/University of Washington/Teaching/SEFS433_Lidar/Labs/"

ALS <- readLAS("Lab 5/LAB5Data/CloudSection.las")
DAP <- readLAS("Lab 6/LAB6Data/DAPClip.las")

DTM <- rast("Lab 6/LAB6Data/VendorDTM_Clip.tif")


?normalize_height


ALSN <- normalize_height(ALS, DTM)
DAPN <- normalize_height(DAP, DTM)

P1 <- clip_circle(ALSN,1636,-2636,60)
P2 <- clip_circle(ALSN,1430,-2230,60)
P3 <- clip_circle(ALSN,1216,-2425,60)
P4 <- clip_circle(ALSN,1279,-2725,60)
P5 <- clip_circle(ALSN,1139,-2174,60)

chm1 <- rasterize_canopy(P1, 1.64, p2r(3.28))
chm1R <- rumple_index(chm1)
chm2 <- rasterize_canopy(P2, 1.64, p2r(3.28))
chm2R <- rumple_index(chm2)
chm3 <- rasterize_canopy(P3, 1.64, p2r(3.28))
chm3R <- rumple_index(chm3)
chm4 <- rasterize_canopy(P4, 1.64, p2r(3.28))
chm4R <- rumple_index(chm4)
chm4 <- rasterize_canopy(P4, 1.64, p2r(3.28))
chm4R <- rumple_index(chm4)
chm5 <- rasterize_canopy(P5, 1.64, p2r(3.28))
chm5R <- rumple_index(chm5)
ALSRumple <- rbind(chm1R,chm2R,chm3R,chm4R,chm5R)
plot(ALSRumple)

#create 5 circle clips with the normalized DAP cloud
D1 <- clip_circle(DAPN,1636,-2636,60)
D2 <- clip_circle(DAPN,1430,-2230,60)
D3<- clip_circle(DAPN,1216,-2425,60)
D4 <- clip_circle(DAPN,1279,-2725,60)
D5 <- clip_circle(DAPN,1139,-2174,60)
DAP_PLOTS <- rbind(D1,D2,D3,D4,D5)


chmD1 <- rasterize_canopy(D1, 1.64, p2r(3.28))
chmD1R <- rumple_index(chmD1)
chmD2 <- rasterize_canopy(D2, 1.64, p2r(3.28))
chmD2R <- rumple_index(chmD2)
chmD3 <- rasterize_canopy(D3, 1.64, p2r(3.28))
chmD3R <- rumple_index(chmD3)
chmD4 <- rasterize_canopy(D4, 1.64, p2r(3.28))
chmD4R <- rumple_index(chmD4)
chmD5 <- rasterize_canopy(D5, 1.64, p2r(3.28))
chmD5R <- rumple_index(chmD5)
DAPRumple <- rbind(chmD1R,chmD2R,chmD3R,chmD4R,chmD5R)
plot(DAPRumple)


reg <- lm(ALSRumple ~ DAPRumple)
summary(reg)


plot(DAPRumple, ALSRumple, xlab = "DAP", ylab = "ALS")
abline(reg, col="blue") #regression line
abline(0,1, col="red", lty = 3) #1 to 1 line


ALSc <- clip_rectangle(ALSN, 700, -3300, 2400, -1700)
DAPc <- clip_rectangle(DAPN, 700, -3300, 2400, -1700)


?pixel_metrics


ALSrasters <- pixel_metrics(ALSc, ~stdmetrics_z(Z), res = 32.8)
plot(ALSrasters)
DAPrasters <- pixel_metrics(DAPNc, ~stdmetrics_z(Z), res = 32.8)  
plot(DAPrasters)


Zdiff <- (ALSrasters$zmax - DAPrasters$zmax)
plot(Zdiff)

plot(ALSrasters$zmax, DAPrasters$zmax)
abline(0,1, col = "red")






?readLAScatalog


ALS <- readLAScatalog("Lab 6/LAB6Data/pack_forest_2013/laz/")

plot(ALS , mapview = T)

DTM <- rast("Lab 6/LAB6Data/pack_forest_2013/dtm/pack_forest_2013_dtm_1.tif")

NEWCRS <- "EPSG:2927"
crs(ALS) <- NEWCRS
crs(DTM) <- NEWCRS


mapview(raster::raster(DTM), map.types = "Esri.WorldImagery")

#the chunk size. 0 means the entire tile. Since the tiles are small
opt_chunk_size(ALS) <- 0 
opt_output_files(ALS) <- "Lab 6/LAB6Data/{ORIGINALFILENAME}_Norm" 
opt_laz_compression(ALS) <- TRUE

ALSN <- normalize_height(ALS, DTM)

plot(ALSN)

opt_chunk_size(ALSN) <- 0
opt_output_files(ALSN) <- "Lab 6/LAB6Data/{ORIGINALFILENAME}_pixels"

opt_filter(ALSN) <- "-drop_z_below 0"

PackCC <- pixel_metrics(ALSN, ~sum(Z>6.56)/sum(Z>=0), 32.8) #canopy cover

mapview(raster::raster(PackZ$zmean_lidar), map.types = "Esri.WorldImagery")



fieldplots <- read.csv("Lab 5/LAB5Data/plotdata.csv")
lidarplots <- read.csv("Lab 5/LAB5Data/cloudmetrics.csv")


BA_field <- fieldplots$BA
zmean_lidar <- (lidarplots$Elev.mean)*3.28
above2percent_lidar <- lidarplots$Percentage.all.returns.above.6.56


BA_mod <- lm(BA_field ~ above2percent_lidar + zmean_lidar)
summary(BA_mod)

names(PackZ) <- c("above2percent_lidar", "zmean_lidar")

BA_predict <- predict(PackZ, BA_mod)
plot(BA_predict, main = "Basal Area")


BA_predict_optimal <- (BA_predict > 160 & PackCC > 0.7)
BA_predict_suboptimal <- (BA_predict < 160 & PackCC < 0.7)


plot(BA_predict_optimal)
mapview(raster::raster(BA_predict_optimal), map.types = "Esri.WorldImagery", col.regions = c("blue", "red"), method = "ngb", layer.name = "Optimal Habitat")

