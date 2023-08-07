setwd("")
#install.packages("lidR")
library(lidR)
library(rgl)


LASfile <- ("Lab 3/ONPlidar/datasetsA/hoh_2013/laz/q47123g8201.laz")
las <- readLAS(LASfile)
plot(las)

print(las)
projection(las)

crs <- st_crs(2927)
projection(las) <- crs
  
  
print(las@header)
las_check(las)

# plot(las, color = "Intensity")
# plot(las, color = 'Intensity', breaks = "quantile", nbreaks = 10)
# plot(las, color = "RGB")


dtm_tin <- rasterize_terrain(las, res = 3.28, algorithm = tin())
plot_dtm3d(dtm_tin)
plot(dtm_tin)
dtm_knnidw <- rasterize_terrain(las, res = 3.28, algorithm = knnidw(k=6L, p = 2))
plot_dtm3d(dtm_knnidw) 


print(las)
las_clip <- clip_circle(las, 797000, 944000, 100)
plot(las_clip)

ws <- seq(3, 12, 3)
th <- seq(0.1, 1.5, length.out = length(ws))
laspmf <- classify_ground(las_clip, pmf(ws, th))
plot(laspmf, color = "Classification")


library(sf)
ctg<- readLAScatalog("Labs/Lab 3/ONPlidar/datasetsA/hoh_2013/laz/")
las_check(ctg)

TileDTM <- rasterize_terrain(ctg, res = 3.28, algorithm = tin())
plot(TileDTM, main = "Tile")
plot_dtm3d(TileDTM)

projection(TileDTM)
terra::crs(TileDTM) <- "EPSG:2927"
projection(TileDTM)

terra::writeRaster(TileDTM, filename = "Lab 3/TileDTM.tif", overwrite=TRUE)


