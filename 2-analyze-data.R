#######################
# 2 - analyze BBS data
#######################


# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(bbsBayes)
library(raster)
library(rstanarm)
library(raster)
library(exactextractr)


# settings ----------------------------------------------------------------

data_pro_date <- '2021-10-12'
run_date <- '2021-10-25'


# load data ----------------------------------------------------------------

data_p <- readRDS(paste0('~/Google_Drive/R/BBS_synchrony/Data/pro-data-', data_pro_date, '.rds'))
#fetch_bbs_data()
strat_data <- stratify(by = 'bbs_cws')



# stats data --------------------------------------------------------------

#filter for only passerines
spstrat <- strat_data$species_strat
data <- dplyr::left_join(data_p, spstrat, by = c('AOU' = 'sp.bbs')) %>%
  dplyr::filter(order == 'Passeriformes')



# calculate mn corr and mean evenness for each route ---------------------------------------------------

uRID <- unique(out$RID)
out <- data.frame(RID = uRID,
                  mn_cc_rs = NA,
                  nsp = NA,
                  nyr = NA,
                  lat = NA,
                  lon = NA,
                  BCR = NA)
for (i in 1:length(uRID))
{
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
  out$mn_cc_rs[i] <- mean(cc_rs[lower.tri(cc_rs)], na.rm = TRUE)
  
  out$nyr[i] <- length(unique(fdat1$Year))
  out$nsp[i] <- length(unique(fdat1$AOU))
  out$lat[i] <- fdat1$Latitude[1]
  out$lon[i] <- fdat1$Longitude[1]
  out$BCR[i] <- fdat1$BCR[1]
}

#save RDS
saveRDS(out, paste0('~/Google_Drive/R/BBS_synchrony/Data/pro-cc-', run_date, '.rds'))
#out <- readRDS(paste0('~/Google_Drive/R/BBS_synchrony/Data/pro-cc-', run_date, '.rds'))


# plot cc on map ----------------------------------------------------------

out2 <- dplyr::filter(out, nsp >= 10, nyr >= 10)
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
plt <- ggplot() +
  geom_sf(data = world, 
          #fill = 'grey89') +
          fill = 'white') +
  #geom_point(data = out, inherit.aes = FALSE, aes(lon, lat, col = mn_cc_lc), size = 1.5, alpha = 0.3) +
  geom_point(data = out2, inherit.aes = FALSE, aes(lon, lat, col = mn_cc_rs), size = 1.5, alpha = 0.3) +
  labs(color = 'Synchrony (rho)') +
  scale_color_gradient2(low = '#C7522B', mid = '#FBF2C4', high = '#3C5941',
                       #earth
                       # scale_fill_gradient2(low = '#A36B2B', mid = '#EDEAC2', high = '#2686A0', 
                       # scale_fill_gradient2(low = '#2686A0', mid = '#EDEAC2', high = '#A36B2B',
                       midpoint = 0) +
  theme_bw() +
  xlim(c(-170, -50)) +
  ylim(c(20, 75))


ggsave(plot = plt, filename = paste0('~/Google_Drive/R/BBS_synchrony/Figures/sync_map_', run_date, '.pdf'), 
       width = 6, height = 4)

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

hfi <- terra::rast('Data/wildareas-v3-1993-human-footprint.tif')
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

ggsave(plot = plt6, filename = paste0('~/Google_Drive/R/BBS_synchrony/Figures/sync_hfi_', run_date, '.pdf'),
       width = 6, height = 4)


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
pp <- matrix(pro_lm_surface$sc_p, nrow = 50, ncol = 50)

pdf(paste0('~/Google_Drive/R/BBS_synchrony/Figures/sync_nsp_hfi_contour_', run_date, '.pdf'), width = 4, height = 4)
contour(axis_x, axis_y, pp, xlab = 'Number of species', ylab = 'HFI')
dev.off()

