setwd(".../LAB10_Data")
library(lidR)
library(mapview)
library(EBImage)
library(rGEDI)

########################################################
################## ALS & TLS ###########################
########################################################

NEWCRS <- sp::CRS("+init=epsg:2927")
PFTLS <- readLAS('.../2019-08-18 Pack Forest TMLS.laz')
HMTLS <- readLAS('.../2019-08-22 Hall of Mosses TMLS.laz')
PFALS <- readLAS('.../WA_031_rp.laz')
HMALS <- readLAS('.../q47123g8201.laz')
crs(PFTLS) <- NEWCRS
crs(HMTLS) <- NEWCRS
crs(PFALS) <- NEWCRS
crs(HMALS) <- NEWCRS

?clip_circle
PFTLSc <- clip_circle(PFTLS,1185950,555180,50)
HMTLSc <- clip_circle(HMTLS,797785,941189,50)
PFALSc <- clip_circle(PFALS,1185950,555180,50)
HMALSc <- clip_circle(HMALS,797785,941189,50)


plot(PFTLSc)
plot(PFALSc)
plot(HMTLSc)
plot(HMALSc)

writeLAS(PFALSc,file = "PFALSc.laz")
writeLAS(PFTLSc,file = "PFTLSc.laz")
writeLAS(HMALSc,file = "HMALSc.laz")
writeLAS(HMTLSc,file = "HMTLSc.laz")

mycsf <- csf(FALSE, 0.2, 0.2, rigidness = 3L)
PFTLSc <- classify_ground(PFTLSc, mycsf)
PFDTM1 <-grid_terrain(PFTLSc, res = 1, algorithm = knnidw())
plot_dtm3d(PFDTM1)

HMTLSc <- classify_ground(HMTLSc, mycsf)
HMDTM1 <-grid_terrain(HMTLSc, res = 1, algorithm = knnidw())
plot_dtm3d(HMDTM1)

PFTLScN <- normalize_height(PFTLSc,PFDTM1)
HMTLScN <- normalize_height(HMTLSc,HMDTM1)
PFALScN <- normalize_height(PFALSc,PFDTM1)
HMALScN <- normalize_height(HMALSc,HMDTM1)
plot(PFTLScN)
plot(PFALScN)
plot(HMTLScN)
plot(HMALScN)

?grid_canopy
PFTLS_DSM <- grid_canopy(PFTLScN,res = 1.5, p2r(2))
PFALS_DSM <- grid_canopy(PFALScN,res = 1.5, p2r(2))  
HMTLS_DSM <- grid_canopy(HMTLScN,res = 1.5, p2r(2))
HMALS_DSM <- grid_canopy(HMALScN,res = 1.5, p2r(2)) 

plot_dtm3d(PFTLS_DSM)
plot_dtm3d(PFALS_DSM)
plot_dtm3d(HMTLS_DSM)
plot_dtm3d(HMALS_DSM)

?rumple_index
rumple_index(PFTLS_DSM)
rumple_index(PFALS_DSM)
rumple_index(HMTLS_DSM)
rumple_index(HMALS_DSM)

?focal
PFTLS_DSMs <- focal(PFTLS_DSM, w = matrix(1,3,3), fun = mean)
PFALS_DSMs <- focal(PFALS_DSM, w = matrix(1,3,3), fun = mean)               
HMTLS_DSMs <- focal(HMTLS_DSM, w = matrix(1,3,3), fun = mean)
HMALS_DSMs <- focal(HMALS_DSM, w = matrix(1,3,3), fun = mean) 
                      
plot_dtm3d(PFTLS_DSMs)
plot_dtm3d(PFALS_DSMs)
plot_dtm3d(HMTLS_DSMs)
plot_dtm3d(HMALS_DSMs)


?segment_trees
PFTLStrees <- segment_trees(PFTLScN, lidR::watershed(PFTLS_DSMs, th = 10))
plot(PFTLStrees,color = "treeID", colorPalette = pastel.colors(100))
PFALStrees <- segment_trees(PFALScN, lidR::watershed(PFALS_DSMs, th = 10))
plot(PFALStrees,color = "treeID", colorPalette = pastel.colors(100))
HMTLStrees <- segment_trees(HMTLScN, lidR::watershed(HMTLS_DSMs, th = 10))
plot(HMTLStrees,color = "treeID", colorPalette = pastel.colors(100))
HMALStrees <- segment_trees(HMALScN, lidR::watershed(HMALS_DSMs, th = 10))
plot(HMALStrees,color = "treeID", colorPalette = pastel.colors(100))

PFTLSHulls  = delineate_crowns(PFTLStrees, func = ~list(Zmax = max(Z)))
print(PFTLSHulls)
mapview(PFTLSHulls,"Zmax")
PFALSHulls  = delineate_crowns(PFALStrees, func = ~list(Zmax = max(Z)))
print(PFALSHulls)
mapview(PFALSHulls,"Zmax")
HMTLSHulls  = delineate_crowns(HMTLStrees, func = ~list(Zmax = max(Z)))
print(HMTLSHulls)
mapview(HMTLSHulls,"Zmax")
HMALSHulls  = delineate_crowns(HMALStrees, func = ~list(Zmax = max(Z)))
print(HMALSHulls)
mapview(HMALSHulls,"Zmax")

?cloud_metrics
PFTLS_cm <- cloud_metrics((filter_poi(PFTLScN,Z > 6.56)), ~stdmetrics_z(Z))
PFTLS_cm
PFALS_cm <- cloud_metrics((filter_poi(PFALScN,Z > 6.56)), ~stdmetrics_z(Z))
PFALS_cm
HMTLS_cm <- cloud_metrics((filter_poi(HMTLScN,Z > 6.56)), ~stdmetrics_z(Z))
HMTLS_cm
HMALS_cm <- cloud_metrics((filter_poi(HMALScN,Z > 6.56)), ~stdmetrics_z(Z))
HMALS_cm

########################################################
################## ALS & GEDI ##########################
########################################################

library(rGEDI)
library(raster)
library(mapview)
library(dplyr)

# Reading GEDI data
?readLevel1B
?readLevel2A
?readLevel2B

#Pack Forest clips for comparison
NEWCRS <- sp::CRS("+init=epsg:2927")
PFALS <- readLAS('.../WA_031_rp.laz')
PFALS2 <- readLAS('.../WA_038_rp.laz')
ONPALS <- readLAS('.../q47123g8206.laz')
PF_DTM <- raster('.../pack_forest_2013_dtm_1.tif')
ONP_DTM <- raster('.../hoh_2013_dtm_4.tif')
ONP_DSM <- raster ('.../hoh_2013_dsm_4.tif')
crs(PFALS2) <- NEWCRS
crs(ONPALS2) <- NEWCRS

?clip_circle
PF_Forested <- clip_circle(PFALS2,1192677,553791,50) #Shot 72230300300100204
PF_Light <- clip_circle(PFALS,1186072,557265,50) #Shot 72230300300100165
PF_ForestedN <- normalize_height(PF_Forested, PF_DTM)
PF_LightN <- normalize_height(PF_Light, PF_DTM)
plot(PF_ForestedN)
plot(PF_LightN)

#### You will compare the following cloud metrics with figures from lab 9#############

?cloud_metrics
PF_Forested_CM <- cloud_metrics(PF_ForestedN,~stdmetrics_z(Z))
PF_Forested_CM
PF_Light_CM <- cloud_metrics(PF_LightN,~stdmetrics_z(Z))
PF_Light_CM

#####################################################################################
  
ONP_Forested <- clip_circle(ONPALS,797610,940475,50) #Shot 64580600200464691
ONP_Light <- clip_circle(ONPALS,797405,937782,50) #Shot 64580800200137203
ONP_ForestedN <- normalize_height(ONP_Forested, ONP_DTM)
ONP_LightN <- normalize_height(ONP_Light, ONP_DTM)
plot(ONP_ForestedN)
plot(ONP_LightN)

#### You will compare the following cloud metrics with figures Made with code below ##

?cloud_metrics
ONP_Forested_CM <- cloud_metrics(ONP_ForestedN,~stdmetrics_z(Z))
ONP_Forested_CM
ONP_Light_CM <- cloud_metrics(ONP_LightN,~stdmetrics_z(Z))
ONP_Light_CM

#####################################################################################
gedilevel1b<-readLevel1B(level1Bpath = file.path ("ONP_GEDI/HM_level1b_clip.h5"))
gedilevel2a<-readLevel2A(level2Apath = file.path ("ONP_GEDI/HM_level2a_clip.h5"))
gedilevel2b<-readLevel2B(level2Bpath = file.path ("ONP_GEDI/HM_level2b_clip.h5"))

#Get GEDI Pulse Full-Waveform Geolocation (GEDI Level1B)
?getLevel1BGeo
level1bGeo<-getLevel1BGeo(level1b=gedilevel1b,select=c("elevation_bin0"))
head(level1bGeo)

# Converting shot_number as "integer64" to "character"
level1bGeo$shot_number<-paste0(level1bGeo$shot_number)

# Converting level1bGeo as data.table to SpatialPointsDataFrame
level1bGeo_spdf<-SpatialPointsDataFrame(cbind(level1bGeo$longitude_bin0, 
                                              level1bGeo$latitude_bin0), data=level1bGeo)

#GEDI shots and ROIs are treated as Geographic (EPSG:4326) coordinate reference system
GEDICRS <- sp::CRS("+init=epsg:4326") #this is WGS84 lat/lon
crs(level1bGeo_spdf)
crs(level1bGeo_spdf)<-GEDICRS

#Get GEDI Elevation and Height Metrics (GEDI Level2A)
# Get GEDI Elevation and Height Metrics
?getLevel2AM
level2AM<-getLevel2AM(gedilevel2a)
head(level2AM[,c("beam","shot_number","elev_highestreturn","elev_lowestmode","rh100")])

# Converting shot_number as "integer64" to "character"
level2AM$shot_number<-paste0(level2AM$shot_number)
# Converting Elevation and Height Metrics as data.table to SpatialPointsDataFrame
level2AM_spdf<-SpatialPointsDataFrame(cbind(level2AM$lon_lowestmode,level2AM$lat_lowestmode),
                                      data=level2AM)

# Exporting Elevation and Height Metrics as ESRI Shapefile
crs(level2AM_spdf)<-GEDICRS
shapefile(level2AM_spdf,paste0("ONP_GEDI/GEDI_2A_EHM"), overwrite=TRUE)

#Plot waveform with RH metrics
?plotWFMetrics
plotWFMetrics(gedilevel1b, gedilevel2a, "64580600200464691", rh=c(25, 50, 75, 90),main="Forest")
plotWFMetrics(gedilevel1b, gedilevel2a, "64580800200137203", rh=c(25, 50, 75, 90),main="Light")

#Get GEDI Vegetation Biophysical Variables (GEDI Level2B)
?getLevel2BVPM
level2BVPM<-getLevel2BVPM(gedilevel2b)
head(level2BVPM[,c("beam","shot_number","pai","fhd_normal","omega","pgap_theta","cover")])

# Converting shot_number as "integer64" to "character"
level2BVPM$shot_number<-paste0(level2BVPM$shot_number)
# Converting GEDI Vegetation Profile Biophysical Variables as data.table to SpatialPointsDataFrame
level2BVPM_spdf<-SpatialPointsDataFrame(cbind(level2BVPM$longitude_lastbin,
                                              level2BVPM$latitude_lastbin),data=level2BVPM)

# Exporting GEDI Vegetation Profile Biophysical Variables as ESRI Shapefile
crs(level2BVPM_spdf)<-GEDICRS
shapefile(level2BVPM_spdf,paste0("ONP_GEDI/GEDI_2B_VPM"),overwrite=TRUE)

Forest <- level2BVPM %>% filter(shot_number =="64580600200464691")
Forest$pai
Lightveg <- level2BVPM %>% filter(shot_number =="64580800200137203")
Lightveg$pai

#Get Plant Area Index (PAI) Profiles (GEDI Level2B)
?getLevel2BPAIProfile
level2BPAIProfile<-getLevel2BPAIProfile(gedilevel2b)
head(level2BPAIProfile[,c("beam","shot_number","pai_z0_5m","pai_z5_10m")])

#Plot Plant Area Index (PAI) Profiles
# Plot Level2B PAI Profile
?plotPAIProfile
#Full power beams are: 0101, 0110, 1000, 1011
gPAIprofile<-plotPAIProfile(level2BPAIProfile, beam="BEAM0101", elev=FALSE)
gPAIprofile<-plotPAIProfile(level2BPAIProfile, beam="BEAM0101", elev=TRUE)
#Coverage beams are: 0000, 0001, 0010, 0011
gPAIprofile<-plotPAIProfile(level2BPAIProfile, beam="BEAM0001", elev=FALSE)
gPAIprofile<-plotPAIProfile(level2BPAIProfile, beam="BEAM0001", elev=TRUE)

# Define your own function
mySetOfMetrics = function(x)
{
  metrics = list(
    min =min(x), # Min of x
    max = max(x), # Max of x
    mean = mean(x), # Mean of x
    sd = sd(x)# Sd of x
  )
  return(metrics)
}

# Computing a serie of statistics of GEDI RH100 metric
# res is in degrees of lat / long. 
?gridStatsLevel2AM
rh100metrics<-gridStatsLevel2AM(level2AM = level2AM, func=mySetOfMetrics(rh100), res=0.005)

?levelplot
plot(rh100metrics$min, main="rh100 Min", xlab = "Longitude (degree)", ylab = "Latitude (degree)")
plot(rh100metrics$max, main="rh100 Max", xlab = "Longitude (degree)", ylab = "Latitude (degree)")
plot(rh100metrics$mean, main="rh100 Mean", xlab = "Longitude (degree)", ylab = "Latitude (degree)")
plot(rh100metrics$sd, main="rh100 SD", xlab = "Longitude (degree)", ylab = "Latitude (degree)")

crs(rh100metrics)
crs(rh100metrics) <- GEDICRS
?writeRaster #help info on 'writeRaster' function
#location and file type for output
writeRaster(rh100metrics$max,filename =file.path("LAB9/rh100_max.tif"),format="GTiff") 

# Computing a serie of statistics of Total Plant Area Index
level2BVPM$pai[level2BVPM$pai==-9999]<-NA # assing NA to -9999
?gridStatsLevel2BVPM
pai_metrics<-gridStatsLevel2BVPM(level2BVPM = level2BVPM, func=mySetOfMetrics(pai), res=0.005)

# View maps
plot(pai_metrics$min, main="pai Min", xlab = "Longitude (degree)", ylab = "Latitude (degree)")
plot(pai_metrics$max, main="pai Max", xlab = "Longitude (degree)", ylab = "Latitude (degree)")
plot(pai_metrics$mean, main="pai Mean", xlab = "Longitude (degree)", ylab = "Latitude (degree)")
plot(pai_metrics$sd, main="pai SD", xlab = "Longitude (degree)", ylab = "Latitude (degree)")

crs(pai_metrics)
crs(pai_metrics) <- GEDICRS
#location and file type for output
writeRaster(pai_metrics$max,filename =file.path("LAB9/PAI_max.tif"),format="GTiff") 


#check if there is a relationship between pai and rh100
plot(pai_metrics$max,rh100metrics$max)


#use the GEDI point data to get values from the ALS DTM
?raster::extract
ONP_DTM_GEDI <- raster::extract(ONP_DTM,level1bGeo_spdf)

#check if there is a relationship between the GEDI and ALS elevation values
reg <- lm(ONP_DTM_GEDI ~ level1bGeo_spdf$elevation_bin0)
summary(reg)
plot(level1bGeo_spdf$elevation_bin0, ONP_DTM_GEDI, xlab="GEDI (m)", ylab="ALS (ft)", 
     main="Elevation") 
abline(reg, col="blue") #regression line

#lets remove the weird GEDI elvevation values that are > 1200
GEDIrh100 <- level1bGeo_spdf$elevation_bin0
GEDIrh100[GEDIrh100 > 1200] <- NA
GEDIrh100
#check if there is a relationship between the corrected GEDI and ALS elevation values
reg <- lm(ONP_DTM_GEDI ~ GEDIrh100)
summary(reg)
plot(GEDIrh100, ONP_DTM_GEDI, xlab="GEDI (m)", ylab="ALS (ft)", 
     main="Elevation") 
abline(reg, col="blue") #regression line

#you can subtract the DTM from the DSM for an easy Canopy Height Model
ONPCHM <- ONP_DSM-ONP_DTM
plot(ONPCHM)

GEDIrh100 <- level2AM_spdf$rh100
GEDIrh100[GEDIrh100 == 0] <- NA

reg <- lm(ONP_DTM_GEDI ~ GEDIrh100)
summary(reg)
plot(GEDIrh100, ONP_DTM_GEDI, xlab="GEDI rh100 (m)", ylab="ALS CHM (ft)", 
     main="Veg Height") 
abline(reg, col="blue")

#You can then use the CHM to test against the RH100 GEDI variable


