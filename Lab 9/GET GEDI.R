setwd("")
library(rGEDI)
#Study area boundary box coordinates
#Hall of moss
ul_lat<- 47.87
lr_lat<- 47.84
ul_lon<- -123.95
lr_lon<- -123.9
daterange=c("2020-05-01","2020-05-22")
# Get path to GEDI data
?gedifinder
gLevel1B<-gedifinder(product="GEDI01_B",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
gLevel2A<-gedifinder(product="GEDI02_A",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
gLevel2B<-gedifinder(product="GEDI02_B",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
gLevel1B
gLevel2A       
gLevel2B
# Set output dir for downloading the files
outdir=("ONP/NEWset")
#Downloading GEDI data
#you need an https://earthexplorer.usgs.gov/ account 
#you can also just go to the link defined above in gLevel1B etc..
gediDownload(filepath=gLevel1B,outdir=outdir)
gediDownload(filepath=gLevel2A,outdir=outdir)
gediDownload(filepath=gLevel2B,outdir=outdir)
## Clip GEDI data by coordinates
# Study area boundary box 
#Hall of moss
ymin<- 47.80
ymax<- 47.9
xmin<- -124
xmax<- -123.8
## clipping GEDI data within boundary box
gedilevel1b<-readLevel1B(level1Bpath = file.path ("ONP/GEDI01_B_2020032225903_O06458_02_T04012_02_005_01_V002.h5"))
gedilevel2a<-readLevel2A(level2Apath = file.path ("ONP/GEDI02_A_2020032225903_O06458_02_T04012_02_003_01_V002.h5"))
gedilevel2b<-readLevel2B(level2Bpath = file.path ("ONP/GEDI02_B_2020032225903_O06458_02_T04012_02_003_01_V002.h5"))
level1b_clip_bb <- clipLevel1B(gedilevel1b, xmin, xmax, ymin, ymax,output=file.path(outdir, "HM_level1b_clip.h5"))
level2a_clip_bb <- clipLevel2A(gedilevel2a, xmin, xmax, ymin, ymax, output=file.path(outdir, "HM_level2a_clip.h5"))
level2b_clip_bb <- clipLevel2B(gedilevel2b, xmin, xmax, ymin, ymax,output=file.path(outdir, "HM_level2b_clip.h5"))



#Pack Forest
ul_lat<- 46.84
lr_lat<- 46.83
ul_lon<- -122.32
lr_lon<- -122.3
daterange=c("2020-04-20","2020-05-22")
# Get path to GEDI data
gLevel1B<-gedifinder(product="GEDI01_B",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
gLevel2A<-gedifinder(product="GEDI02_A",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
gLevel2B<-gedifinder(product="GEDI02_B",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
gLevel1B
gLevel2A
gLevel2B
# Set output dir for downloading the files
outdir=("PackForestGEDI")
#Downloading GEDI data
#you need an https://earthexplorer.usgs.gov/ account 
#you can also just go to the link defined above in gLevel1B etc..
gediDownload(filepath=gLevel1B,outdir=outdir)
gediDownload(filepath=gLevel2A,outdir=outdir)
gediDownload(filepath=gLevel2B,outdir=outdir)
## Clip GEDI data by coordinates
# Study area boundary box 
#Pack Forest
ymin<- 46.8
ymax<- 46.9
xmin<- -122.4
xmax<- -122.2
## clipping GEDI data within boundary box
gedilevel1b<-readLevel1B(level1Bpath = file.path ("PackForestGEDI/GEDI01_B_2020082072049_O07223_03_T02974_02_005_01_V002.h5"))
gedilevel2a<-readLevel2A(level2Apath = file.path ("PackForestGEDI/GEDI02_A_2020082072049_O07223_03_T02974_02_003_01_V002.h5"))
gedilevel2b<-readLevel2B(level2Bpath = file.path ("PackForestGEDI/GEDI02_B_2020082072049_O07223_03_T02974_02_003_01_V002.h5"))
level1b_clip_bb <- clipLevel1B(gedilevel1b, xmin, xmax, ymin, ymax,output=file.path(outdir, "PF_level1b_clip.h5"))
level2a_clip_bb <- clipLevel2A(gedilevel2a, xmin, xmax, ymin, ymax, output=file.path(outdir, "PF_level2a_clip.h5"))
level2b_clip_bb <- clipLevel2B(gedilevel2b, xmin, xmax, ymin, ymax,output=file.path(outdir, "PF_level2b_clip.h5"))


