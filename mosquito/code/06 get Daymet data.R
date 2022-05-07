site_df<-read_csv(paste0(path, "/data/NEON_Field_Site_Metadata_20220412.csv")) %>% 
  dplyr::select(Site=field_site_id, lat=field_latitude, lon=field_longitude) %>% 
  filter(Site %in% unique(full.df$Site))

if(!file.exists( "./mosquito/data/daymet.rds")) {
  library(daymetr)
  daymet_df_list<-vector(mode="list", length=nrow(site_df))
  for (r in 1:nrow(site_df)) {
    daymet_df_list[[r]] <- download_daymet(site = site_df$Site[r],
                                           lat = site_df$lat[r],
                                           lon = site_df$lon[r],
                                           start = 2014,
                                           end = 2019,
                                           internal = TRUE,
                                           simplify = TRUE) %>% 
      mutate(date=as.Date(yday, origin = as.Date(paste0(year,"-01-01")))) %>% 
      spread(key="measurement", value="value") %>% 
      dplyr::select(date, tmax=`tmax..deg.c.`, tmin=`tmin..deg.c.`, prcp=`prcp..mm.day.`) %>% 
      mutate(temp=(tmax+ tmin)/2) %>% 
      dplyr::select(date, temp, prcp) %>% 
      mutate(Site=site_df$Site[r])
  }
  daymet_df<-bind_rows(daymet_df_list) %>% 
    mutate(Year=format(date, "%Y") %>% as.numeric()) %>% 
    group_by(Year, Site) %>% 
    summarise(mat=mean(temp),
              tap=sum(prcp))
  write_rds(daymet_df, "./mosquito/data/daymet.rds")
}  else {
  daymet_df<-read_rds("./mosquito/data/daymet.rds")
}
