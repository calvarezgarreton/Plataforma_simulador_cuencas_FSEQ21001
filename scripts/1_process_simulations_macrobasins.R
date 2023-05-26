rm(list=ls()); a=getSrcDirectory(function(x){x}) ; setwd(a)
options(scipen = 999) # to avoid scientific nomenclature
options(device="quartz")
quartz.options(dpi=100)
graphics.off()
library(hydroGOF)
library(hydroTSM)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2)
library(reshape2)

source('~/Dropbox/mHM/mHM_results/mhm_read_functions.R')
source('aux_init_variables.R')

dir.create('../figs')
dir.create('../results')

############ SIMULATION HISTORICA #######
main_path_results = '~/Dropbox/mHM/mhm_experiments/2023_05_platform_runs'
main_path_input   = '~/Dropbox/Data/mhm_data_platform'

simulations = data.frame(rbind(c(exp_folder_name='historical_run_macro_basins',
                                 sim_short_name='hist_cr2met_lai_ltm_sat',
                                 sim_pr='meteo_cr2met_v2.5_pr024_pet_hsm_10spin_1960_2021',
                                 sim_pr_name='cr2met',
                                 sim_pet='meteo_cr2met_v2.5_pr024_pet_hsm_10spin_1960_2021',
                                 sim_pet_name='cr2met_hsm',
                                 sim_txn='meteo_cr2met_v2.5_pr024_pet_hsm_10spin_1960_2021',
                                 sim_txn_name='cr2met',
                                 sim_lai_name='AVHRR_VGT_LTM')
                               ))

############ SELECT EXPERIMENT #######
# j=1
nexp=nrow(simulations)

for (j in 1:nexp){
  
  exp_folder_name  = simulations$exp_folder_name[j]
  sim_short_name   = simulations$sim_short_name[j]
  sim_pr_name      = simulations$sim_pr_name[j]
  sim_pet_name     = simulations$sim_pet_name[j]
  sim_txn_name     = simulations$sim_txn_name[j]
  
  sim_pr_path      = simulations$sim_pr[j]
  sim_pet_path     = simulations$sim_pet[j]
  sim_txn_path     = simulations$sim_txn[j]
  sim_lai_name     = simulations$sim_lai[j]
  
  path_results     = paste0(main_path_results,'/',exp_folder_name,'/results')
  ids              = list.files(path = paste0(path_results), all.files = FALSE,full.names = FALSE, recursive = FALSE,include.dirs = FALSE)
  
  # i=ids[1]
  for (i in ids){
    path = paste0(path_results)
    
    if (file.exists(paste0(path,'/',i,'/mHM_Fluxes_States.nc'))){

    attrib_i = attrib[attrib$gauge_id==i,]
    
    #### load forcing data ####
    path_input_i = paste0(main_path_input,'/',i,'/',sim_pr_path,'/')
    varname = 'pre'
    pr_day  = mhm.day.nc.input.to.day.ts(path_input_i,varname)
    pr      = daily2monthly(pr_day,FUN=sum)
    
    path_input_i = paste0(main_path_input,'/',i,'/',sim_pet_path,'/')
    varname = 'pet'
    pet_in_day = mhm.day.nc.input.to.day.ts(path_input_i,varname)
    pet_in     = daily2monthly(pet_in_day,FUN=sum)
    
    path_input_i = paste0(main_path_input,'/',i,'/',sim_txn_path,'/')
    varname = 'tmin'
    tmin_day = mhm.day.nc.input.to.day.ts(path_input_i,varname)
    tmin     = daily2monthly(tmin_day,FUN=mean)
    varname = 'tmax'
    tmax_day = mhm.day.nc.input.to.day.ts(path_input_i,varname)
    tmax     = daily2monthly(tmax_day,FUN=mean)
    varname = 'tavg'
    tavg_day = mhm.day.nc.input.to.day.ts(path_input_i,varname) 
    tavg     = daily2monthly(tavg_day,FUN=mean)
    
    #### load simulated fluxes ####
    
    nc.path       = paste0(path,'/',i,'/mHM_Fluxes_States.nc') 
    et_mon        = mhm.nc.out.to.ts.accum.mon(nc.path,'aET')
    index(et_mon) <- as.Date(time(et_mon))
    
    date_ini    = time(et_mon)[1]
    date_end    = time(et_mon)[length(et_mon)]
    
    qsim        = mhm.qsim.to.ts(paste0(path,'/',i),date_ini,date_end)
    qobs        = mhm.qobs.to.ts(paste0(path,'/',i),date_ini,date_end)
    qobs[qobs<0]<- NA
    qobs_mon    = daily2monthly(qobs,FUN=mean)
    qsim_mon    = daily2monthly(qsim,FUN=mean)
    
    # truncate time series
    pr_mon      = window(pr, start = date_ini, end = date_end)
    pet_mon     = window(pet_in, start = date_ini, end = date_end)
    tmin_mon    = window(tmin, start = date_ini, end = date_end)
    tmax_mon    = window(tmax, start = date_ini, end = date_end)
    tavg_mon    = window(tavg, start = date_ini, end = date_end)
    
    pr_yr    = monthly2annual(pr_mon,FUN=sum)
    pet_yr   = monthly2annual(pet_mon,FUN=sum)
    tmin_yr  = monthly2annual(tmin_mon,FUN=mean)
    tmax_yr  = monthly2annual(tmax_mon,FUN=mean)
    tavg_yr  = monthly2annual(tavg_mon,FUN=mean)
    qobs_yr  = monthly2annual(qobs_mon,FUN=mean)
    qsim_yr  = monthly2annual(qsim_mon,FUN=mean)
    et_yr    = monthly2annual(et_mon,FUN=sum)
    
    nazoo <- zoo(NA*(1:12))

    # processing historical simulations
    
    if (date_ini<'1970-01-01' & date_end<'2050-01-01'){
      
    dini='1960-01-01'
    dend='2020-12-31'
    pr_mon_LTM_1960_2020    = monthlyfunction(window(pr_mon, start = dini, end = dend),FUN=mean)
    pet_mon_LTM_1960_2020   = monthlyfunction(window(pet_mon, start = dini, end = dend),FUN=mean)
    tmin_mon_LTM_1960_2020  = monthlyfunction(window(tmin_mon, start = dini, end = dend),FUN=mean)
    tmax_mon_LTM_1960_2020  = monthlyfunction(window(tmax_mon, start = dini, end = dend),FUN=mean)
    tavg_mon_LTM_1960_2020  = monthlyfunction(window(tavg_mon, start = dini, end = dend),FUN=mean)
    qsim_mon_LTM_1960_2020  = monthlyfunction(window(qsim_mon, start = dini, end = dend),FUN=mean)
    et_mon_LTM_1960_2020    = monthlyfunction(window(et_mon, start = dini, end = dend),FUN=mean)

    if (sum(!is.na(window(qobs_mon, start = dini, end = dend)))/length(is.na(window(qobs_mon, start = dini, end = dend)))>0.8){
    qobs_mon_LTM_1960_2020  = monthlyfunction(window(qobs_mon, start = dini, end = dend),FUN=mean,na.rm=T)
    }else{
    qobs_mon_LTM_1960_2020  = monthlyfunction(window(qobs_mon, start = dini, end = dend),FUN=mean,na.rm=F)
    }
    
    dini='1960-01-01'
    dend='1990-12-31'
    pr_mon_LTM_1960_1990    = monthlyfunction(window(pr_mon, start = dini, end = dend),FUN=mean)
    pet_mon_LTM_1960_1990   = monthlyfunction(window(pet_mon, start = dini, end = dend),FUN=mean)
    tmin_mon_LTM_1960_1990  = monthlyfunction(window(tmin_mon, start = dini, end = dend),FUN=mean)
    tmax_mon_LTM_1960_1990  = monthlyfunction(window(tmax_mon, start = dini, end = dend),FUN=mean)
    tavg_mon_LTM_1960_1990  = monthlyfunction(window(tavg_mon, start = dini, end = dend),FUN=mean)
    qsim_mon_LTM_1960_1990  = monthlyfunction(window(qsim_mon, start = dini, end = dend),FUN=mean)
    et_mon_LTM_1960_1990    = monthlyfunction(window(et_mon, start = dini, end = dend),FUN=mean)

    if (sum(!is.na(window(qobs_mon, start = dini, end = dend)))/length(is.na(window(qobs_mon, start = dini, end = dend)))>0.8){
      qobs_mon_LTM_1960_1990  = monthlyfunction(window(qobs_mon, start = dini, end = dend),FUN=mean,na.rm=T)
    }else{
      qobs_mon_LTM_1960_1990  = monthlyfunction(window(qobs_mon, start = dini, end = dend),FUN=mean,na.rm=F)
    }
    
    dini='1990-01-01'
    dend='2020-12-31'
    pr_mon_LTM_1990_2020    = monthlyfunction(window(pr_mon, start = dini, end = dend),FUN=mean)
    pet_mon_LTM_1990_2020   = monthlyfunction(window(pet_mon, start = dini, end = dend),FUN=mean)
    tmin_mon_LTM_1990_2020  = monthlyfunction(window(tmin_mon, start = dini, end = dend),FUN=mean)
    tmax_mon_LTM_1990_2020  = monthlyfunction(window(tmax_mon, start = dini, end = dend),FUN=mean)
    tavg_mon_LTM_1990_2020  = monthlyfunction(window(tavg_mon, start = dini, end = dend),FUN=mean)
    qsim_mon_LTM_1990_2020  = monthlyfunction(window(qsim_mon, start = dini, end = dend),FUN=mean)
    et_mon_LTM_1990_2020    = monthlyfunction(window(et_mon, start = dini, end = dend),FUN=mean)

    if (sum(!is.na(window(qobs_mon, start = dini, end = dend)))/length(is.na(window(qobs_mon, start = dini, end = dend)))>0.8){
      qobs_mon_LTM_1990_2020  = monthlyfunction(window(qobs_mon, start = dini, end = dend),FUN=mean,na.rm=T)
    }else{
      qobs_mon_LTM_1990_2020  = monthlyfunction(window(qobs_mon, start = dini, end = dend),FUN=mean,na.rm=F)
    }
    
    pr_mon_LTM_2040_2070    = nazoo
    pet_mon_LTM_2040_2070   = nazoo
    tmin_mon_LTM_2040_2070  = nazoo
    tmax_mon_LTM_2040_2070  = nazoo
    tavg_mon_LTM_2040_2070  = nazoo
    qsim_mon_LTM_2040_2070  = nazoo
    et_mon_LTM_2040_2070    = nazoo
    qobs_mon_LTM_2040_2070  = nazoo
    
    pr_mon_LTM_2070_2100    = nazoo
    pet_mon_LTM_2070_2100   = nazoo
    tmin_mon_LTM_2070_2100  = nazoo
    tmax_mon_LTM_2070_2100  = nazoo
    tavg_mon_LTM_2070_2100  = nazoo
    qsim_mon_LTM_2070_2100  = nazoo
    et_mon_LTM_2070_2100    = nazoo
    qobs_mon_LTM_2070_2100  = nazoo
    
    }else{
      pr_mon_LTM_1960_2020=nazoo
      et_mon_LTM_1960_2020=nazoo
      tmin_mon_LTM_1960_2020=nazoo
      tmax_mon_LTM_1960_2020=nazoo
      tavg_mon_LTM_1960_2020=nazoo
      pet_mon_LTM_1960_2020=nazoo
      qobs_mon_LTM_1960_2020=nazoo
      qsim_mon_LTM_1960_2020=nazoo
      
      pr_mon_LTM_1960_1990=nazoo
      et_mon_LTM_1960_1990=nazoo
      tmin_mon_LTM_1960_1990=nazoo
      tmax_mon_LTM_1960_1990=nazoo
      tavg_mon_LTM_1960_1990=nazoo
      pet_mon_LTM_1960_1990=nazoo
      qobs_mon_LTM_1960_1990=nazoo
      qsim_mon_LTM_1960_1990=nazoo

      pr_mon_LTM_1990_2020=nazoo
      et_mon_LTM_1990_2020=nazoo
      tmin_mon_LTM_1990_2020=nazoo
      tmax_mon_LTM_1990_2020=nazoo
      tavg_mon_LTM_1990_2020=nazoo
      pet_mon_LTM_1990_2020=nazoo
      qobs_mon_LTM_1990_2020=nazoo
      qsim_mon_LTM_1990_2020=nazoo
      
      dini='2040-01-01'
      dend='2070-12-31'
      pr_mon_LTM_2040_2070    = monthlyfunction(window(pr_mon, start = dini, end = dend),FUN=mean)
      pet_mon_LTM_2040_2070   = monthlyfunction(window(pet_mon, start = dini, end = dend),FUN=mean)
      tmin_mon_LTM_2040_2070  = monthlyfunction(window(tmin_mon, start = dini, end = dend),FUN=mean)
      tmax_mon_LTM_2040_2070  = monthlyfunction(window(tmax_mon, start = dini, end = dend),FUN=mean)
      tavg_mon_LTM_2040_2070  = monthlyfunction(window(tavg_mon, start = dini, end = dend),FUN=mean)
      qsim_mon_LTM_2040_2070  = monthlyfunction(window(qsim_mon, start = dini, end = dend),FUN=mean)
      et_mon_LTM_2040_2070    = monthlyfunction(window(et_mon, start = dini, end = dend),FUN=mean)
      qobs_mon_LTM_2040_2070  = nazoo
    
      dini='2070-01-01'
      dend='2099-12-31'
      pr_mon_LTM_2070_2100    = monthlyfunction(window(pr_mon, start = dini, end = dend),FUN=mean)
      pet_mon_LTM_2070_2100   = monthlyfunction(window(pet_mon, start = dini, end = dend),FUN=mean)
      tmin_mon_LTM_2070_2100  = monthlyfunction(window(tmin_mon, start = dini, end = dend),FUN=mean)
      tmax_mon_LTM_2070_2100  = monthlyfunction(window(tmax_mon, start = dini, end = dend),FUN=mean)
      tavg_mon_LTM_2070_2100  = monthlyfunction(window(tavg_mon, start = dini, end = dend),FUN=mean)
      qsim_mon_LTM_2070_2100  = monthlyfunction(window(qsim_mon, start = dini, end = dend),FUN=mean)
      et_mon_LTM_2070_2100    = monthlyfunction(window(et_mon, start = dini, end = dend),FUN=mean)
      qobs_mon_LTM_2070_2100  = nazoo
    
    }
    
    temp_mon          = data.frame(gauge_id=attrib_i$gauge_id,
                                   gauge_name=attrib_i$gauge_name,
                                   area_km2=attrib_i$area_km2,
                                   date=time(pr_mon),
                                   pr_mm_mon=pr_mon,
                                   et_mm_mon=et_mon,
                                   tmin_C_mon=tmin_mon,
                                   tmax_C_mon=tmax_mon,
                                   tavg_C_mon=tavg_mon,
                                   pet_mm_mon=pet_mon,
                                   qobs_m3s=qobs_mon,
                                   qobs_mm_mon=NA,
                                   qsim_m3s=qsim_mon,
                                   qsim_mm_mon=NA,
                                   sim_short_name=sim_short_name,
                                   sim_pr_name=sim_pr_name,
                                   sim_pet_name=sim_pet_name,
                                   sim_txn_name=sim_txn_name,
                                   sim_lai_name=sim_lai_name
    )
    
    m3s_to_mm_mon        = (1/attrib_i$area_km2)*(1/1000000)*(3600*24*30)*1000
    temp_mon$qsim_mm_mon = temp_mon$qsim_m3s*m3s_to_mm_mon
    temp_mon$qobs_mm_mon = temp_mon$qobs_m3s*m3s_to_mm_mon
    
    temp_yr           = data.frame(gauge_id=attrib_i$gauge_id,
                                   gauge_name=attrib_i$gauge_name,
                                   area_km2=attrib_i$area_km2,
                                   date=time(pr_yr),
                                   pr_mm_yr=pr_yr,
                                   et_mm_yr=et_yr,
                                   tmin_C_yr=tmin_yr,
                                   tmax_C_yr=tmax_yr,
                                   tavg_C_yr=tavg_yr,
                                   pet_mm_yr=pet_yr,
                                   qobs_m3s_yr=qobs_yr,
                                   qobs_mm_yr=NA,
                                   qsim_m3s_yr=qsim_yr,
                                   qsim_mm_yr=NA,
                                   sim_short_name=sim_short_name,
                                   sim_pr_name=sim_pr_name,
                                   sim_pet_name=sim_pet_name,
                                   sim_txn_name=sim_txn_name,
                                   sim_lai_name=sim_lai_name
    )
    
    m3s_to_mm_yr       = (1/attrib_i$area_km2)*(1/1000000)*(3600*24*365)*1000
    temp_yr$qsim_mm_yr = temp_yr$qsim_m3s*m3s_to_mm_yr
    temp_yr$qobs_mm_yr = temp_yr$qobs_m3s*m3s_to_mm_yr
    
    
    temp_mon_LTM          = data.frame(gauge_id=attrib_i$gauge_id,
                                       gauge_name=attrib_i$gauge_name,
                                       area_km2=attrib_i$area_km2,
                                       date=time(pr_mon_LTM_1960_2020),
                                       
                                       sim_short_name=sim_short_name,
                                       sim_pr_name=sim_pr_name,
                                       sim_pet_name=sim_pet_name,
                                       sim_txn_name=sim_txn_name,
                                       sim_lai_name=sim_lai_name,
                                       
                                       pr_mm_mon_LTM_1960_2020=pr_mon_LTM_1960_2020,
                                       et_mm_mon_LTM_1960_2020=et_mon_LTM_1960_2020,
                                       tmin_C_mon_LTM_1960_2020=tmin_mon_LTM_1960_2020,
                                       tmax_C_mon_LTM_1960_2020=tmax_mon_LTM_1960_2020,
                                       tavg_C_mon_LTM_1960_2020=tavg_mon_LTM_1960_2020,
                                       pet_mm_mon_LTM_1960_2020=pet_mon_LTM_1960_2020,
                                       qobs_m3s_mon_LTM_1960_2020=qobs_mon_LTM_1960_2020,
                                       qobs_mm_mon_LTM_1960_2020=NA,
                                       qsim_m3s_mon_LTM_1960_2020=qsim_mon_LTM_1960_2020,
                                       qsim_mm_mon_LTM_1960_2020=NA,
                                       
                                       pr_mm_mon_LTM_1960_1990=pr_mon_LTM_1960_1990,
                                       et_mm_mon_LTM_1960_1990=et_mon_LTM_1960_1990,
                                       tmin_C_mon_LTM_1960_1990=tmin_mon_LTM_1960_1990,
                                       tmax_C_mon_LTM_1960_1990=tmax_mon_LTM_1960_1990,
                                       tavg_C_mon_LTM_1960_1990=tavg_mon_LTM_1960_1990,
                                       pet_mm_mon_LTM_1960_1990=pet_mon_LTM_1960_1990,
                                       qobs_m3s_mon_LTM_1960_1990=qobs_mon_LTM_1960_1990,
                                       qobs_mm_mon_LTM_1960_1990=NA,
                                       qsim_m3s_mon_LTM_1960_1990=qsim_mon_LTM_1960_1990,
                                       qsim_mm_mon_LTM_1960_1990=NA,
                                       
                                       pr_mm_mon_LTM_1990_2020=pr_mon_LTM_1990_2020,
                                       et_mm_mon_LTM_1990_2020=et_mon_LTM_1990_2020,
                                       tmin_C_mon_LTM_1990_2020=tmin_mon_LTM_1990_2020,
                                       tmax_C_mon_LTM_1990_2020=tmax_mon_LTM_1990_2020,
                                       tavg_C_mon_LTM_1990_2020=tavg_mon_LTM_1990_2020,
                                       pet_mm_mon_LTM_1990_2020=pet_mon_LTM_1990_2020,
                                       qobs_m3s_mon_LTM_1990_2020=qobs_mon_LTM_1990_2020,
                                       qobs_mm_mon_LTM_1990_2020=NA,
                                       qsim_m3s_mon_LTM_1990_2020=qsim_mon_LTM_1990_2020,
                                       qsim_mm_mon_LTM_1990_2020=NA,

                                       pr_mm_mon_LTM_2040_2070=pr_mon_LTM_2040_2070,
                                       et_mm_mon_LTM_2040_2070=et_mon_LTM_2040_2070,
                                       tmin_C_mon_LTM_2040_2070=tmin_mon_LTM_2040_2070,
                                       tmax_C_mon_LTM_2040_2070=tmax_mon_LTM_2040_2070,
                                       tavg_C_mon_LTM_2040_2070=tavg_mon_LTM_2040_2070,
                                       pet_mm_mon_LTM_2040_2070=pet_mon_LTM_2040_2070,
                                       qobs_m3s_mon_LTM_2040_2070=qobs_mon_LTM_2040_2070,
                                       qobs_mm_mon_LTM_2040_2070=NA,
                                       qsim_m3s_mon_LTM_2040_2070=qsim_mon_LTM_2040_2070,
                                       qsim_mm_mon_LTM_2040_2070=NA,
                                       
                                       pr_mm_mon_LTM_2070_2100=pr_mon_LTM_2070_2100,
                                       et_mm_mon_LTM_2070_2100=et_mon_LTM_2070_2100,
                                       tmin_C_mon_LTM_2070_2100=tmin_mon_LTM_2070_2100,
                                       tmax_C_mon_LTM_2070_2100=tmax_mon_LTM_2070_2100,
                                       tavg_C_mon_LTM_2070_2100=tavg_mon_LTM_2070_2100,
                                       pet_mm_mon_LTM_2070_2100=pet_mon_LTM_2070_2100,
                                       qobs_m3s_mon_LTM_2070_2100=qobs_mon_LTM_2070_2100,
                                       qobs_mm_mon_LTM_2070_2100=NA,
                                       qsim_m3s_mon_LTM_2070_2100=qsim_mon_LTM_2070_2100,
                                       qsim_mm_mon_LTM_2070_2100=NA
                                      

    )
    
    temp_mon_LTM$qsim_mm_mon_LTM_1960_2020 = temp_mon_LTM$qsim_m3s_mon_LTM_1960_2020*m3s_to_mm_mon
    temp_mon_LTM$qobs_mm_mon_LTM_1960_2020 = temp_mon_LTM$qobs_m3s_mon_LTM_1960_2020*m3s_to_mm_mon
    
    temp_mon_LTM$qsim_mm_mon_LTM_1960_1990 = temp_mon_LTM$qsim_m3s_mon_LTM_1960_1990*m3s_to_mm_mon
    temp_mon_LTM$qobs_mm_mon_LTM_1960_1990 = temp_mon_LTM$qobs_m3s_mon_LTM_1960_1990*m3s_to_mm_mon

    temp_mon_LTM$qsim_mm_mon_LTM_1990_2020 = temp_mon_LTM$qsim_m3s_mon_LTM_1990_2020*m3s_to_mm_mon
    temp_mon_LTM$qobs_mm_mon_LTM_1990_2020 = temp_mon_LTM$qobs_m3s_mon_LTM_1990_2020*m3s_to_mm_mon

    temp_mon_LTM$qsim_mm_mon_LTM_2040_2070 = temp_mon_LTM$qsim_m3s_mon_LTM_2040_2070*m3s_to_mm_mon
    temp_mon_LTM$qobs_mm_mon_LTM_2040_2070 = temp_mon_LTM$qobs_m3s_mon_LTM_2040_2070*m3s_to_mm_mon

    temp_mon_LTM$qsim_mm_mon_LTM_2070_2100 = temp_mon_LTM$qsim_m3s_mon_LTM_2070_2100*m3s_to_mm_mon
    temp_mon_LTM$qobs_mm_mon_LTM_2070_2100 = temp_mon_LTM$qobs_m3s_mon_LTM_2070_2100*m3s_to_mm_mon
    
    data_platform_mon = structure(rbind(data_platform_mon,temp_mon),.Names = names(data_platform_mon))
    data_platform_mon_LTM = structure(rbind(data_platform_mon_LTM,temp_mon_LTM),.Names = names(data_platform_mon_LTM))
    data_platform_yr = structure(rbind(data_platform_yr,temp_yr),.Names = names(data_platform_yr))
    
    
    path_figs = paste0('../figs/sim_',sim_short_name,'_qsim_gauge_',i,'.pdf')
    
    # processing historical simulations
    
    if (date_ini<'1970-01-01' & date_end<'2050-01-01'){
      
    pdf(file = path_figs) #
    # dev.new()
    layout(matrix(c(1,1,2,3),nrow=2,ncol = 2, byrow = TRUE))
    # par(mfrow = c(2,1),lwd = 1)
    plot.new()
    plot.window(xlim = range((time(qsim_mon))), ylim = range(qsim_mon, qobs_mon,na.rm=T))
    axis(2, at = axTicks(2), labels =axTicks(2), tck = 0, cex.lab=.9,cex.axis=.8)
    box(cex=.8)
    lines(time(qsim_mon),qsim_mon, col = "dodgerblue", lwd = 1.2, type = "l")
    lines(time(qsim_mon),qobs_mon, col = "grey27", lwd = 1, type = "l")
    title(main = paste0(attrib_i$gauge_name, ', Lat = ', round(attrib_i$gauge_lat,1), ' Mean Elev = ', round(attrib_i$mean_elev,0),'\n',
                        'Sim: ',sim_short_name,'Sim. period: ', date_ini, '-',date_end, ' gauge_id: ',i), cex.main = 0.7)
    title(xlab = "", ylab = "Streamflow (m3/s)", cex.lab = 0.8)
    xticks <- seq((as.Date(start(qsim_mon))), (as.Date(end(qsim_mon))), by = "month")
    legend("topright", legend = c("Sim", "Obs"), col = c("dodgerblue", "grey27"), lwd = 1, cex = 0.8, bty = "n")
    axis(1, at = xticks, labels = format(as.POSIXct(xticks), "%m-%Y"), tck = 0, cex.lab=.9,cex.axis=.8)
    mtext(paste0('daily KGE (1965-2021)= ', as.numeric(gof(qsim_mon,qobs_mon)["KGE",])),cex=.7,line=-1,adj = 0.02)
    
    plot(1:12,temp_mon_LTM$qsim_mm_mon_LTM_1960_1990, col = "red", lwd = 1.2, type = "l",ylab="",xaxt='n',xlab='',
         ylim = range(0,temp_mon_LTM$pr_mm_mon_LTM_1960_1990,temp_mon_LTM$pr_mm_mon_LTM_1990_2020))
    lines(1:12,temp_mon_LTM$qobs_mm_mon_LTM_1960_1990, col = "black", lwd = 1.2, type = "l")
    lines(1:12,temp_mon_LTM$pr_mm_mon_LTM_1960_1990, col = "blue", lwd = 1.2, type = "l")
    lines(1:12,temp_mon_LTM$et_mm_mon_LTM_1960_1990, col = "orange", lwd = 1.2, type = "l")
    title(main = 'Mean monthly flows 1960-1990', cex.main = 0.7)
    title(xlab = "", ylab = "Mean flow (mm)", cex.lab = 0.8)
    xticks <- time(pr_mon_LTM_1960_1990)
    legend("topright", legend = c("Qsim", "Qobs","Pr","ET"), col = c("red", "black","blue","orange"), lwd = 1, cex = 0.8, bty = "n")
    axis(1, at = 1:12, labels = time(pr_mon_LTM_1960_1990), tck = 0, cex.lab=.9,cex.axis=.8)

    plot(1:12,temp_mon_LTM$qsim_mm_mon_LTM_1990_2020, col = "red", lwd = 1.2, type = "l",ylab="",xaxt='n',xlab='',
         ylim = range(0,temp_mon_LTM$pr_mm_mon_LTM_1960_1990,temp_mon_LTM$pr_mm_mon_LTM_1990_2020))
    lines(1:12,temp_mon_LTM$qobs_mm_mon_LTM_1990_2020, col = "black", lwd = 1.2, type = "l")
    lines(1:12,temp_mon_LTM$pr_mm_mon_LTM_1990_2020, col = "blue", lwd = 1.2, type = "l")
    lines(1:12,temp_mon_LTM$et_mm_mon_LTM_1990_2020, col = "orange", lwd = 1.2, type = "l")
    title(main = 'Mean monthly flows 1990-2020', cex.main = 0.7)
    title(xlab = "", ylab = "Mean flow (mm)", cex.lab = 0.8)
    xticks <- time(pr_mon_LTM_1990_2020)
    legend("topright", legend = c("Qsim", "Qobs","Pr","ET"), col = c("red", "black","blue","orange"), lwd = 1, cex = 0.8, bty = "n")
    axis(1, at = 1:12, labels = time(pr_mon_LTM_1990_2020), tck = 0, cex.lab=.9,cex.axis=.8)
    
    dev.off()
    }else{
      pdf(file = path_figs) #
      # dev.new()
      layout(matrix(c(1,1,2,3),nrow=2,ncol = 2, byrow = TRUE))
      # par(mfrow = c(2,1),lwd = 1)
      plot.new()
      plot.window(xlim = range((time(qsim_mon))), ylim = range(qsim_mon, qobs_mon,na.rm=T))
      axis(2, at = axTicks(2), labels =axTicks(2), tck = 0, cex.lab=.9,cex.axis=.8)
      box(cex=.8)
      lines(time(qsim_mon),qsim_mon, col = "dodgerblue", lwd = 1.2, type = "l")
      lines(time(qsim_mon),qobs_mon, col = "grey27", lwd = 1, type = "l")
      title(main = paste0(attrib_i$gauge_name, ', Lat = ', round(attrib_i$gauge_lat,1), ' Mean Elev = ', round(attrib_i$mean_elev,0),'\n',
                          'Sim: ',sim_short_name,'Sim. period: ', date_ini, '-',date_end, ' gauge_id: ',i), cex.main = 0.7)
      title(xlab = "", ylab = "Streamflow (m3/s)", cex.lab = 0.8)
      xticks <- seq((as.Date(start(qsim_mon))), (as.Date(end(qsim_mon))), by = "month")
      legend("topright", legend = c("Sim", "Obs"), col = c("dodgerblue", "grey27"), lwd = 1, cex = 0.8, bty = "n")
      axis(1, at = xticks, labels = format(as.POSIXct(xticks), "%m-%Y"), tck = 0, cex.lab=.9,cex.axis=.8)
      mtext(paste0('daily KGE (1965-2021)= ', as.numeric(gof(qsim_mon,qobs_mon)["KGE",])),cex=.7,line=-1,adj = 0.02)
      
      plot(1:12,temp_mon_LTM$qsim_mm_mon_LTM_2040_2070, col = "red", lwd = 1.2, type = "l",ylab="",xaxt='n',xlab='',
           ylim = range(0,temp_mon_LTM$pr_mm_mon_LTM_2040_2070,temp_mon_LTM$pr_mm_mon_LTM_2070_2100))
      lines(1:12,temp_mon_LTM$qobs_mm_mon_LTM_2040_2070, col = "black", lwd = 1.2, type = "l")
      lines(1:12,temp_mon_LTM$pr_mm_mon_LTM_2040_2070, col = "blue", lwd = 1.2, type = "l")
      lines(1:12,temp_mon_LTM$et_mm_mon_LTM_2040_2070, col = "orange", lwd = 1.2, type = "l")
      title(main = 'Mean monthly flows 2040-2070', cex.main = 0.7)
      title(xlab = "", ylab = "Mean flow (mm)", cex.lab = 0.8)
      xticks <- time(pr_mon_LTM_2040_2070)
      legend("topright", legend = c("Qsim", "Qobs","Pr","ET"), col = c("red", "black","blue","orange"), lwd = 1, cex = 0.8, bty = "n")
      axis(1, at = 1:12, labels = time(pr_mon_LTM_2040_2070), tck = 0, cex.lab=.9,cex.axis=.8)
      
      plot(1:12,temp_mon_LTM$qsim_mm_mon_LTM_2070_2100, col = "red", lwd = 1.2, type = "l",ylab="",xaxt='n',xlab='',
           ylim = range(0,temp_mon_LTM$pr_mm_mon_LTM_2040_2070,temp_mon_LTM$pr_mm_mon_LTM_2070_2100))
      lines(1:12,temp_mon_LTM$qobs_mm_mon_LTM_2070_2100, col = "black", lwd = 1.2, type = "l")
      lines(1:12,temp_mon_LTM$pr_mm_mon_LTM_2070_2100, col = "blue", lwd = 1.2, type = "l")
      lines(1:12,temp_mon_LTM$et_mm_mon_LTM_2070_2100, col = "orange", lwd = 1.2, type = "l")
      title(main = 'Mean monthly flows 2070-2100', cex.main = 0.7)
      title(xlab = "", ylab = "Mean flow (mm)", cex.lab = 0.8)
      xticks <- time(pr_mon_LTM_2070_2100)
      legend("topright", legend = c("Qsim", "Qobs","Pr","ET"), col = c("red", "black","blue","orange"), lwd = 1, cex = 0.8, bty = "n")
      axis(1, at = 1:12, labels = time(pr_mon_LTM_2070_2100), tck = 0, cex.lab=.9,cex.axis=.8)
      
      dev.off()
      
    }
    }
  } 
}

ids_macro_basin_platform= unique(data_platform_mon_LTM$gauge_id)
basin       <- readOGR('~/Dropbox/Research_projects/AA_CAMELS_CL/Updating_CAMELScl/v2021/processed_data/v2021_12/boundaries',layer='catchments_camels_cl_v2021')
platform_macro_basin <- basin[basin$gauge_id %in% ids_macro_basin_platform,]
path        <- "../results"
writeOGR(obj=platform_macro_basin, layer = 'platform_macro_basin', dsn=path , driver="ESRI Shapefile",overwrite_layer = TRUE)


write.csv(data_platform_mon,file=paste0('../results/data_macrobasin_platform_mon.csv'),row.names = FALSE)
write.csv(data_platform_yr,file=paste0('../results/data_macrobasin_platform_yr.csv'),row.names = FALSE)
write.csv(data_platform_mon_LTM,file=paste0('../results/data_macrobasin_platform_mon_LTM.csv'),row.names = FALSE)

