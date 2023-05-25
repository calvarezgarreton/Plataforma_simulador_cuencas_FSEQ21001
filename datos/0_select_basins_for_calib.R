# select basins afor mHM calibration

rm(list=ls()); a=getSrcDirectory(function(x){x}) ; setwd(a)
library(readr)
library(rgdal)
library(rgeos)


###### Filter based on DAMS & AREA ###### 
load("/Users/cag/Dropbox/Research_projects/AA_CAMELS_CL/Updating_CAMELScl/v2021/processed_data/v2021_12/camels_v2021_attributes.RData")
ids=attrib
id_processed = list.files('/Volumes/lacie_cami/mHM/gcp_data/mhm_data_basin',recursive = FALSE,full.names=FALSE)
ids=ids[ids$gauge_id %in% id_processed,]
ids=ids[ids$dam_index==0 ,]
# ids=ids[ids$lc_crop<10,]
ids=ids[ids$area_km2>100,]
id_filt_attrib = ids$gauge_id

#######  Filter based on BUDYKO analysis ###### 
# basins_budyko <- read_csv("~/Dropbox/mHM/data_processing/select_basins_calib/budyko/camels-cl_gauges_budyko_analysis_tol_0_from_1979_to_2019.csv")
basins_budyko <- read_csv("/Users/cag/Dropbox/mHM/data_processing/2023_04_consolidated_flow/basin_calib_selection/budyko_analysis_Rdirectory/out/camels-cl_gauges_budyko_analysis_tol_0.8_from_1985_to_2018.csv")
id_filt_bud   <- basins_budyko$gauge_id[basins_budyko$meets_budyko==1]
id_calib      <- id_filt_attrib[id_filt_attrib %in% id_filt_bud]


#######  Filter based on data availability ###### 

load('/Users/cag/Dropbox/Research_projects/AA_CAMELS_CL/Updating_CAMELScl/v2021/processed_data/v2021_12/camels_v2021_qflx_m3s_day.RData')
qdate  = as.Date(paste0(as.character(q_m3s_day$year),'-',as.character(q_m3s_day$month),'-',as.character(q_m3s_day$day)))
qdata  = q_m3s_day[,4:519]
qspin = qdata[qdate>'1981-01-01' & qdate<('1986-01-01'),]
nobspin = colSums(!is.na(qspin))/365
qcalib = qdata[qdate>'1986-01-01' & qdate<('2020-01-01'),] #lai data is only up to 31/12/2019
nobscalib = colSums(!is.na(qcalib))/365
filt_reg = data.frame('gauge_id'=colnames(qdata),'yr_obs_calib'=nobscalib,'yr_obs_spin'=nobspin)

basins_filtered_mhm = basins_budyko
basins_filtered_mhm$meets_dam_filter=1
basins_filtered_mhm$meets_area_filter=1
basins_filtered_mhm$meets_record_filter=1
basins_filtered_mhm$combined_filter=1
# basins_filtered_mhm$calibrated=1

basins_filtered_mhm = merge(basins_filtered_mhm,attrib[,c('gauge_id','dam_index','area_km2','n_obs','record_period_start')],by='gauge_id')
basins_filtered_mhm = merge(basins_filtered_mhm,filt_reg,by='gauge_id')
basins_filtered_mhm$meets_area_filter[basins_filtered_mhm$area_km2<=100]<-0
basins_filtered_mhm$meets_dam_filter[basins_filtered_mhm$dam_index==1]<-0

# basins_filtered_mhm$meets_record_filter[basins_filtered_mhm$yr_obs_calib<(10)|
#                                           basins_filtered_mhm$yr_obs_spin<(3)]<-0
# 
basins_filtered_mhm$meets_record_filter[basins_filtered_mhm$yr_obs_calib<(15)]<-0

basins_filtered_mhm$combined_filter[basins_filtered_mhm$meets_area_filter==0 |
                                      basins_filtered_mhm$meets_dam_filter==0 |
                                      basins_filtered_mhm$meets_record_filter==0 |
                                      basins_filtered_mhm$meets_budyko==0] <- 0

sum(basins_filtered_mhm$combined_filter==1)

####### Manually add chi2 basins and other that should be calibrated #####
basins_filtered_mhm$combined_filter[basins_filtered_mhm$gauge_id%in%c(5710001,5722002,5100001,5101001,7336001,7339001,9414001)]<-1
sum(basins_filtered_mhm$combined_filter==1)


### save selection ###
basins_filtered_mhm_order <- basins_filtered_mhm[order(basins_filtered_mhm$area_km2),]
id_calib=data.frame(basins_filtered_mhm_order$gauge_id[basins_filtered_mhm_order$combined_filter==1 & basins_filtered_mhm_order$gauge_id %in% id_processed])
id_calib_area=data.frame(basins_filtered_mhm_order[basins_filtered_mhm_order$combined_filter==1 & basins_filtered_mhm_order$gauge_id %in% id_processed,])
dev.new();hist(id_calib_area$area_km2)

id_calib_below1000=id_calib_area$gauge_id[id_calib_area$area_km2<1000]
id_calib_above1000=id_calib_area$gauge_id[id_calib_area$area_km2>=1000]

id_calib_below1000_sample = c(id_calib_below1000[1:10],id_calib_below1000[80:100])
id_calib_below1000_sample = id_calib_above1000[1:20]
ids_chi2=c(5100001,5722002,7339001,9414001)

dim(id_calib)
length(id_calib_below1000)
length(id_calib_above1000)

write.table(id_calib,file='id_calib_20230414.txt',row.names = FALSE,col.names = FALSE)
write.table(id_calib_below1000,file='id_calib_20230414_below1000.txt',row.names = FALSE,col.names = FALSE)
write.table(id_calib_above1000,file='id_calib_20230414_above1000.txt',row.names = FALSE,col.names = FALSE)


#######  Write to shapefile ###### 
# basin       <- readOGR('/Users/calvarez/Dropbox/Research_projects/2018_CAMELS_CL/Updating_CAMELScl/processed_data/v2021/boundaries',layer='catchments_camels_cl_v2021')
# basin_calib <- basin[basin$gauge_id %in% id_calib,]
# path        <- "shp"

# writeOGR(obj=basin_calib, layer = 'basin_calib', dsn=path , driver="ESRI Shapefile",overwrite_layer = TRUE)
# write.table(id_calib, file='ids_basin_calib_20211111.txt', append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
