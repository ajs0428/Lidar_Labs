setwd("OneDrive - UW/University of Washington/Teaching/SEFS433_Lidar/Labs/")

library(lidR)
library(rgl)
library(sf)
library(terra)


las <- readLAS("Lab 5/LAB5Data/Sequoia.las")
lasNORM <- normalize_height(las, tin())
writeLAS(lasNORM, "Lab 5/LAB5Data/SequoiaNORM.las")

plot(lasNORM)
?cloud_metrics

CMSequoia <- as.data.frame(cloud_metrics(lasNORM, .stdmetrics))
?stdmetrics

CMSequoia$zmax

hist(lasNORM$Z, breaks = 100)

over2m <- filter_poi(lasNORM, Z>=6.56)
CMSequoia_over2m <- as.data.frame(cloud_metrics(over2m,.stdmetrics))
hist(over2m$Z, breaks = 100)

plot(over2m)


las_proj <- readLAS("Lab 5/LAB5Data/CloudSection.las")
plot(las_proj)

P1 <- clip_circle(las_proj, 1636, -2636, 60)
P1n <- normalize_height(P1, tin())

P2 <- clip_circle(las_proj, 1430, -2230, 60)
P2n <- normalize_height(P2, tin())

P3 <- clip_circle(las_proj, 1216, -2425, 60)
P3n <- normalize_height(P3, tin())

P4 <- clip_circle(las_proj, 1279, -2725, 60)
P4n <- normalize_height(P4, tin())

P5 <- clip_circle(las_proj, 1139, -2174, 60)
P5n <- normalize_height(P5, tin())

lasPLOTS <- rbind(P1n, P2n, P3n, P4n, P5n)
plot(lasPLOTS)



P1m <- cloud_metrics(filter_poi(P1n, Z > 6.56), .stdmetrics)
P2m <- cloud_metrics(filter_poi(P2n, Z>6.56), .stdmetrics)
P3m <- cloud_metrics(filter_poi(P3n, Z>6.56), .stdmetrics)
P4m <- cloud_metrics(filter_poi(P4n, Z>6.56), .stdmetrics)
P5m <- cloud_metrics(filter_poi(P5n, Z>6.56), .stdmetrics)


P_CM <- rbind.data.frame(P1m, P2m, P3m, P4m, P5m)

# Canopy Cover 
P1c <- cloud_metrics(P1n, ~sum(Z>6.56)/sum(Z>-1))
P2c <- cloud_metrics(P2n, ~sum(Z>6.56)/sum(Z>-1))
P3c <- cloud_metrics(P3n, ~sum(Z>6.56)/sum(Z>-1))
P4c <- cloud_metrics(P4n, ~sum(Z>6.56)/sum(Z>-1))
P5c <- cloud_metrics(P5n, ~sum(Z>6.56)/sum(Z>-1))

P_CC <- rbind(P1c, P2c, P3c, P4c, P5c)

plot2 <- read.csv("Lab 5/LAB5Data/plot2.csv")

df <- cbind.data.frame(plot2, P_CC, P_CM)


plot(df$zq95, df$BASAL)
abline(lm(df$BASAL ~ df$zq95))

summary(lm(df$BASAL ~ df$zq95))


plotdata <- read.csv("Lab 5/LAB5Data/plotdata.csv")
lidardata <- read.csv("Lab 5/LAB5Data/cloudmetrics.csv")

str(plotdata)


lmod <- lm(plotdata$BA.GT5 ~ lidardata$Percentage.all.returns.above.6.56)

plot(lidardata$Percentage.all.returns.above.6.5, plotdata$BA.GT5)
abline(lm(plotdata$BA.GT5 ~ lidardata$Percentage.all.returns.above.6.56))

summary(lmod)


lmod_all <- lm(plotdata$BA ~  lidardata$Percentage.all.returns.above.6.56 +
                   lidardata$Elev.mean 
                   )

summary(lmod_all)
