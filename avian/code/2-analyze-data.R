#######################
# 2 - analyze BBS data
#######################

# load packages -----------------------------------------------------------

p_load(tidyverse)
p_load(sf)
p_load(bbsBayes)
p_load(raster)
p_load(rstanarm)
p_load(raster)
p_load(exactextractr)
p_load(rnaturalearth)
p_load(gridExtra)
p_load(ggpubr)
# load data ----------------------------------------------------------------
path_in<-"./avian/data/processed/"
files<-list.files(path_in, full.names = T)
m_data_list<-vector(mode="list", length=length(files))
for (f in 1:length(files)) {
  file<-files[f]
  m_data_list[[f]]<-read_rds(file)
  print(f)
}
data_p<-bind_rows(m_data_list)
load("./avian/data/bbs_raw_data.RData")
strat_data <- stratify(by = 'bbs_cws', bbs_data = bbs_data)

# time series plot
tplt <- dplyr::filter(data_p, RID == uRID[1])
p_ts <- ggplot(tplt, aes(Year, residuals, color = factor(AOU))) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  theme(legend.pos = 'none') +
  xlab('year') +
  ylab('log(Count) residuals') +
  ggtitle(paste0(tplt$RouteName[1], ' (', round(tplt$Latitude[1], 2), ', ', round(tplt$Longitude[1], 2), ')'))
# stats data --------------------------------------------------------------

#filter for only passerines
spstrat <- strat_data$species_strat
data <- dplyr::left_join(data_p, spstrat, by = c('AOU' = 'sp.bbs')) %>%
  dplyr::filter(order == 'Passeriformes')



# calculate mn corr and mean evenness for each route ---------------------------------------------------

uRID <- data %>% pull(RID) %>% unique() %>% sort()

p_load(foreach)
p_load(doSNOW)
num_cores <- 50
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

out_list<-
foreach (i = 1:length(uRID),
         .packages = c("tidyverse")) %dopar% {
  #i <- 236
  print(paste0('processing route: ', i, ' of ', length(uRID)))
  
  fdat1 <- dplyr::filter(data, RID == uRID[i])
  
  #mean when multiple obs for a given AOU/year
  fdat2 <- suppressMessages(fdat1 %>%
    dplyr::group_by(AOU, Year) %>%
    dplyr::summarize(log_count = mean(log_count), resid = mean(residuals),
                     count = mean(SpeciesTotal)))
  
  #wide format
  fdat_wide_rs <- fdat2 %>%
    dplyr::select(c(-log_count, -count)) %>%
    tidyr::pivot_wider(names_from = AOU, values_from = resid)
  
  #mean of corr matrix - community synchrony
  cc_rs <- cor(fdat_wide_rs[,-1], use = 'pairwise.complete.obs')
  
  out <- data.frame(RID = uRID[i],
                    mn_cc_rs = mean(cc_rs[lower.tri(cc_rs)], na.rm = TRUE),
                    nsp = length(unique(fdat1$AOU)),
                    nyr = length(unique(fdat1$Year)),
                    lat = fdat1$Latitude[1],
                    lon = fdat1$Longitude[1],
                    BCR = fdat1$BCR[1])
  out
}
out<-bind_rows(out_list)

#save RDS
write_rds(out, './avian/data/pro-cc.rds')
# out<-read_rds('./avian/data/pro-cc.rds')


# plot cc on map ----------------------------------------------------------

out2 <- out %>% filter(nsp >= 10, nyr >= 10)
world <- rnaturalearth::ne_coastline(returnclass = 'sf')
# circle <- st_point(x = c(0,0)) %>% st_buffer(dist = 10000000) %>% 
#   st_sfc(crs = 4326)

#unique route IDs
uRID2 <- unique(out2$RID)

#number of routes
length(uRID2)

#number of species
mean(out2$nsp)
range(out2$nsp)

#extent time
mean(out2$nyr)
range(out2$nyr)


#map of synchrony across space

p_map <- ggplot() +
  geom_sf(data = world, 
          #fill = 'grey89') +
          fill = 'white') +
  #geom_point(data = out, inherit.aes = FALSE, aes(lon, lat, col = mn_cc_lc), size = 1.5, alpha = 0.3) +
  geom_point(data = out2 %>% mutate(rank=rank(mn_cc_rs)), inherit.aes = FALSE, aes(lon, lat, col = rank), size = 1.5, alpha = 0.3) +
  labs(color = 'Synchrony (rho)') +
  # scale_color_gradient2(low = '#C7522B', mid = '#FBF2C4', high = '#3C5941',
  #                      #earth
  #                      # scale_fill_gradient2(low = '#A36B2B', mid = '#EDEAC2', high = '#2686A0',
  #                      # scale_fill_gradient2(low = '#2686A0', mid = '#EDEAC2', high = '#A36B2B',
  #                      midpoint = 0) +
  scale_color_viridis_c( breaks = out2 %>% 
                           mutate(rank=rank(mn_cc_rs)) %>%
                           arrange(mn_cc_rs) %>% 
                           mutate(cut=cut(mn_cc_rs, 
                                          breaks=c(min(mn_cc_rs),c(0, 0.1, 0.2, 0.3, 0.6), max(mn_cc_rs)),
                                          include.lowest=T)) %>% 
                           group_by(cut) %>% 
                           summarise(rank=max(rank)+0.5) %>%
                           slice(-n()) %>% 
                           pull(rank), 
                         label=c(0, 0.1, 0.2, 0.3, 0.6))+
  geom_point(data=tplt %>% distinct(Longitude, Latitude), aes(x=Longitude, y=Latitude), col="red", pch=10, cex=3)+
  xlab("longitude")+
  ylab("latitude")+
  theme_minimal() +
  theme(legend.position="bottom")+
  xlim(c(-170, -50)) +
  ylim(c(20, 75))


pdf('./avian/output/sync_map.pdf', width = 12, height = 10)
print(p_map)
dev.off()

#histogram of community synchrony
hist(out2$mn_cc_rs, breaks = 20)
range(out2$mn_cc_rs)
mean(out2$mn_cc_rs)


# human footprint ---------------------------------------------------------

#BBS route 25 miles long - get buffer (40 km) when extracting vars
#source: https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-1993-human-footprint/data-download#close

pts_sp <- sp::SpatialPoints(cbind(out2$lon, out2$lat))
sf_pts <- sf::st_as_sf(pts_sp)
st_crs(sf_pts) <- 4326

hfi <- terra::rast('./avian/data/wildareas-v3-2009-human-footprint-geotiff/wildareas-v3-2009-human-footprint.tif')
sf_pts3 <- sf::st_transform(sf_pts, crs = crs(hfi)) %>%
  st_buffer(dist = 40000)
eval3 <- exactextractr::exact_extract(hfi, sf_pts3, 
                                      fun = 'mean')

#add to previous cc dataframe - HFI values are between 0 and 50
out5 <- dplyr::mutate(out2, HFI = eval3) %>%
  dplyr::filter(HFI <= 50)

plt6 <- ggplot(out5, aes(x = HFI, y = mn_cc_rs, color = factor(BCR))) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  xlab('Human Footprint Index') +
  ylab('Synchrony (rho)') +
  geom_line(stat = 'smooth', method = 'lm', size = 3, alpha = 0.5)

pdf('./avian/output/sync_hfi.pdf',width = 6, height = 4)
print(plt6)
dev.off()

out5$nsp <- as.numeric(out5$nsp)

ff <- rstanarm::stan_gamm4(mn_cc_rs ~ s(nsp) + HFI,
                         data = out5,
                         chains = 4,
                         cores = 4,
                         control = list(adapt_delta = 0.99))
MCMCvis::MCMCsummary(ff, round = 5)


#Contour plot
pred_y <- predict(ff, newdata = data.frame(nsp = mean(out5$nsp), HFI = 0:50))

axis_x <- seq(min(out5$nsp), max(out5$nsp), length.out = 50)
axis_y <- seq(min(out5$HFI), max(out5$HFI), length.out = 50)
pro_lm_surface <- expand.grid(nsp = axis_x, HFI = axis_y, KEEP.OUT.ATTRS = F)
pro_lm_surface$sc_p <- predict(ff, newdata = pro_lm_surface)

p_contour<- ggplot(pro_lm_surface)+
  geom_contour_filled(aes(x=HFI, y=nsp, z=sc_p))+
  geom_contour(aes(x=HFI, y=nsp, z=sc_p))+
  metR::geom_label_contour(aes(x=HFI, y=nsp, z=sc_p), skip = 0)+
  ylab ('Number of species')+
  xlab("Human Footprint Index (HFI)")+
  guides(fill="none")+
  theme_classic()

pdf('./avian/output/sync_nsp_hfi_contour.pdf', width = 4, height = 4)
print(p_contour)
dev.off()

pdf('./avian/output/ts_map_cont.pdf', width = 12, height = 8)
grid.arrange(annotate_figure(p_ts, fig.lab = "A"),
             annotate_figure(p_map, fig.lab = "B"),
             annotate_figure(p_contour, fig.lab = "C"),
             layout_matrix=rbind(c(1,1,1),
                                 c(2,2,3),
                                 c(2,2,3))
             )
dev.off()
