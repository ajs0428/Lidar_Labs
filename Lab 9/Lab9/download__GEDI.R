library(rGEDI)
library(sf)

stdy <- st_read(r"(D:\Courses\CEWA_567\project\data\study_area\study_area.shp)")
stdy <- st_transform(stdy, 4326)
ext <- st_bbox(stdy)

ul_lat<- ext["ymax"]
lr_lat<- ext["ymin"]
ul_lon<- ext["xmin"]
lr_lon<- ext["xmax"]

daterange1=c("2020-06-01","2020-08-01")
daterange2=c("2021-06-01","2021-08-01")

gLevel2A_pre<-gedifinder(product="GEDI02_A",ul_lat, ul_lon, lr_lat, lr_lon,version="002",daterange=daterange1)
gLevel2B_pre<-gedifinder(product="GEDI02_B",ul_lat, ul_lon, lr_lat, lr_lon,version="002",daterange=daterange1)

gLevel2A_post<-gedifinder(product="GEDI02_A",ul_lat, ul_lon, lr_lat, lr_lon,version="002",daterange=daterange2)
gLevel2B_post<-gedifinder(product="GEDI02_B",ul_lat, ul_lon, lr_lat, lr_lon,version="002",daterange=daterange2)

outdir <- r"(D:\Courses\CEWA_567\project\data\gedi\new_download\)"

gediDownload(filepath=gLevel2A_pre,outdir=outdir)
gediDownload(filepath=gLevel2B_pre,outdir=outdir)
gediDownload(filepath=gLevel2A_post,outdir=outdir)
gediDownload(filepath=gLevel2B_post,outdir=outdir)

# Converting

fils <- list.files(r"(D:\Courses\CEWA_567\project\data\gedi\new_download\)" ,
                   pattern = ".h5$", full.names = T)

##

fils2a <- fils[grepl("_A_", fils)]

ints <- list()

this_gedi <- readLevel2A(level2Apath = fils2a[[1]])
this_gedi <- getLevel2AM(this_gedi)
this_gedi$shot_number<-paste0(this_gedi$shot_number)
this_gedi <- this_gedi[complete.cases(this_gedi),]
this_gedi <- st_as_sf(this_gedi, coords = c("lon_lowestmode", "lat_lowestmode"), crs = 4326)
this_gedi <- dplyr::filter(this_gedi, degrade_flag != 1)
this_gedi <- dplyr::filter(this_gedi, quality_flag == 1)

all_pnts <- this_gedi

for (i in 2:length(fils2a)) {
  this_gedi <- readLevel2A(level2Apath = fils2a[[i]])
  this_gedi <- getLevel2AM(this_gedi)
  this_gedi$shot_number<-paste0(this_gedi$shot_number)
  this_gedi <- this_gedi[complete.cases(this_gedi),]
  this_gedi <- st_as_sf(this_gedi, coords = c("lon_lowestmode", "lat_lowestmode"), crs = 4326)
  this_gedi <- dplyr::filter(this_gedi, degrade_flag != 1)
  this_gedi <- dplyr::filter(this_gedi, quality_flag == 1)
  
  if(nrow(this_gedi) > 0) {
    all_pnts <- rbind(all_pnts, this_gedi)
    ints <- c(ints, i)
  }
  print(i/length(fils2a))
}

final_pnts_2a <- all_pnts

## level 2

fils2b <- fils[grepl("_B_", fils)]

ints <- list()

this_gedi <- readLevel2B(level2Bpath = fils2b[[1]])
this_gedi <- getLevel2BVPM(this_gedi)
this_gedi$shot_number<-paste0(this_gedi$shot_number)
this_gedi <- this_gedi[complete.cases(this_gedi),]
this_gedi <- st_as_sf(this_gedi, coords = c("longitude_lastbin", "latitude_lastbin"), crs = 4326)
this_gedi <- dplyr::filter(this_gedi, l2b_quality_flag == 1)
this_gedi <- dplyr::filter(this_gedi, algorithmrun_flag == 1)

all_pnts <- this_gedi

for (i in 2:length(fils2b)) {
  this_gedi <- readLevel2B(level2Bpath = fils2b[[i]])
  this_gedi <- getLevel2BVPM(this_gedi)
  this_gedi$shot_number<-paste0(this_gedi$shot_number)
  this_gedi <- this_gedi[complete.cases(this_gedi),]
  this_gedi <- st_as_sf(this_gedi, coords = c("longitude_lastbin", "latitude_lastbin"), crs = 4326)
  this_gedi <- dplyr::filter(this_gedi, l2b_quality_flag == 1)
  this_gedi <- dplyr::filter(this_gedi, algorithmrun_flag == 1)
  
  if(nrow(this_gedi) > 0) {
    all_pnts <- rbind(all_pnts, this_gedi)
    ints <- c(ints, i)
  }
  print(i/length(fils2b))
}

final_pnts_2b <- all_pnts

## joining

st_geometry(final_pnts_2b) <- NULL

final_pnts <- dplyr::inner_join(final_pnts_2a, final_pnts_2b, by = "shot_number")
final_pnts <- dplyr::select(final_pnts, !contains(".y"))
nw_names <- gsub(".x", "", colnames(final_pnts))
names(final_pnts) <- nw_names

final_pnts <- final_pnts[st_intersects(final_pnts, stdy, sparse = F),]

st_write(final_pnts, r"(D:\Courses\CEWA_567\project\data\gedi\gedi_points_study_area\cewa567_gedi.shp)")

