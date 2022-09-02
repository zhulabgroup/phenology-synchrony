#######################
# 1 - process BBS data
#######################


# load packages -----------------------------------------------------------
# JAGS downloaded from https://sourceforge.net/projects/mcmc-jags/files/
library(bbsBayes)
library(tidyverse)


# get data ----------------------------------------------------------------

# only needs to be run once
# fetch_bbs_data()
path <- "./case study B: avian/"
load(paste0(path, "data/bbs_raw_data.RData"))

strat_data <- stratify(by = "bbs_cws", bbs_data = bbs_data)

# process
bstrat <- strat_data$bird_strat

# routes are unique within each state (i.e., need to create new)
RID_df <- bstrat %>%
  dplyr::distinct(statenum, Route) %>%
  dplyr::arrange(statenum, Route) %>%
  dplyr::mutate(RID = 1:NROW(.))

# route info
route_df <- strat_data$route_strat %>%
  dplyr::select(
    statenum, Route, RouteName,
    Latitude, Longitude
  ) %>%
  dplyr::distinct()

# merge new RID (unique ID for each route, irrespective of state) with full data
bstrat2 <- left_join(bstrat, RID_df, by = c("statenum", "Route")) %>%
  left_join(route_df, by = c("statenum", "Route")) %>%
  dplyr::mutate(log_count = log(SpeciesTotal)) %>%
  dplyr::arrange(RID, AOU, Year) %>%
  dplyr::select(RID, RouteName, AOU, Year, SpeciesTotal, log_count, BCR, Latitude, Longitude)

# #plot
# bstrat2 %>%
#   select(Latitude, Longitude) %>%
#   distinct() %>%
#   ggplot(aes(Longitude, Latitude)) +
#   geom_point() +
#   theme_bw()


# process data -----------------------------------------------------------------

# unique RIDs
uRID <- unique(bstrat2$RID)
file_list <- paste0(path, "data/processed/") %>% list.files(pattern = ".rds", full.names = T)
if (length(file_list) < 54) {
  # RESIDUALS MIGHT BE OFF (for route 5146 at least was 0)
  # INCLUDE ALL SPECIES, NOT JUST COMPLETE TIME SERIES (determine what to do about zeros in time series) - can filter after using mean count, etc.
  # MEAN WHEN MULTIPLE OBS PER YEAR
  # at least 5 species for each RID

  library(foreach)
  library(doSNOW)
  num_cores <- 50
  cl <- makeCluster(num_cores)
  registerDoSNOW(cl)
  m_data_list <-
    foreach(
      i = 1:length(uRID),
      .packages = c("tidyverse")
    ) %dopar% {
      # i <- 770
      print(paste0("processing route: ", i, " of ", length(uRID)))
      # data for single route
      r_data <- dplyr::filter(bstrat2, RID == uRID[i])

      # only species seen in all years
      sp_keep1 <- r_data %>%
        dplyr::group_by(AOU) %>%
        dplyr::summarize(nyrs = length(unique(Year))) %>%
        dplyr::ungroup() %>%
        dplyr::filter(nyrs == length(unique(r_data$Year)))

      # filtered by 'abundant' species
      r_data2 <- r_data %>%
        dplyr::filter(AOU %in% sp_keep1$AOU)

      # sp_keep1 <- data.frame(AOU = unique(r_data$AOU))

      # get mean count of each species across years
      sp_keep2 <- r_data2 %>%
        dplyr::group_by(AOU) %>%
        dplyr::summarize(mn_count = mean(SpeciesTotal)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(AOU)

      r_data3 <- dplyr::left_join(r_data2, sp_keep2, by = "AOU") %>%
        dplyr::arrange(AOU, Year)
      # r_data2 <- r_data


      # get residuals from linear regression (on log count)
      resids <- r_data3 %>%
        group_by(AOU) %>%
        do(model = residuals(lm(log_count ~ Year, data = .)))

      # add to df
      r_data3$residuals <- unlist(resids$model)

      if (NROW(r_data3) >= 5) {
        r_data3
      }
    }
  m_data <- bind_rows(m_data_list)

  # save RDS ----------------------------------------------------------------
  year_list <- m_data %>%
    pull(Year) %>%
    unique()
  path_out <- paste0(path, "data/processed/")
  dir.create(path_out, recursive = T)
  for (year in year_list) {
    write_rds(m_data %>% filter(Year == year), paste0(path_out, year, ".rds"))
    print(year)
  }
}

m_data_list <- vector(mode = "list")
for (file in file_list) {
  m_data_list[[file]] <- read_rds(file)
}
m_data <- bind_rows(m_data_list)


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
library(rnaturalearth)
library(gridExtra)
library(ggpubr)
# load data ----------------------------------------------------------------

# filter for only passerines
spstrat <- strat_data$species_strat
tplt <- dplyr::left_join(m_data, spstrat, by = c("AOU" = "sp.bbs")) %>%
  dplyr::filter(order == "Passeriformes") %>%
  dplyr::filter(RID == uRID[1])

# time series plot
p_ts <- ggplot(tplt, aes(Year, residuals, color = factor(AOU))) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  theme(legend.pos = "none") +
  xlab("Year") +
  ylab("log (count) residuals") #+
p_ts
# stats data --------------------------------------------------------------



# calculate mn corr and mean evenness for each route ---------------------------------------------------

if (!file.exists(paste0(path, "data/pro-cc.rds"))) {
  library(foreach)
  library(doSNOW)
  num_cores <- 50
  cl <- makeCluster(num_cores)
  registerDoSNOW(cl)

  out_list <-
    foreach(
      i = 1:length(uRID),
      .packages = c("tidyverse")
    ) %dopar% {
      # i <- 236
      print(paste0("processing route: ", i, " of ", length(uRID)))

      fdat1 <- dplyr::filter(data, RID == uRID[i])

      # mean when multiple obs for a given AOU/year
      fdat2 <- suppressMessages(fdat1 %>%
        dplyr::group_by(AOU, Year) %>%
        dplyr::summarize(
          log_count = mean(log_count), resid = mean(residuals),
          count = mean(SpeciesTotal)
        ))

      # wide format
      fdat_wide_rs <- fdat2 %>%
        dplyr::select(c(-log_count, -count)) %>%
        tidyr::pivot_wider(names_from = AOU, values_from = resid)

      # mean of corr matrix - community synchrony
      cc_rs <- cor(fdat_wide_rs[, -1], use = "pairwise.complete.obs")

      out <- data.frame(
        RID = uRID[i],
        mn_cc_rs = mean(cc_rs[lower.tri(cc_rs)], na.rm = TRUE),
        nsp = length(unique(fdat1$AOU)),
        nyr = length(unique(fdat1$Year)),
        lat = fdat1$Latitude[1],
        lon = fdat1$Longitude[1],
        BCR = fdat1$BCR[1]
      )
      out
    }
  out <- bind_rows(out_list)

  # save RDS
  write_rds(out, paste0(path, "data/pro-cc.rds"))
} else {
  out <- read_rds(paste0(path, "data/pro-cc.rds"))
}

# plot cc on map ----------------------------------------------------------

out2 <- out %>% filter(nsp >= 10, nyr >= 10)
world <- rnaturalearth::ne_coastline(returnclass = "sf")
world <- as(world, "Spatial")
# circle <- st_point(x = c(0,0)) %>% st_buffer(dist = 10000000) %>%
#   st_sfc(crs = 4326)

# unique route IDs
uRID2 <- unique(out2$RID)

# number of routes
length(uRID2)

# number of species
mean(out2$nsp)
range(out2$nsp)

# extent time
mean(out2$nyr)
range(out2$nyr)


# map of synchrony across space

p_map <- ggplot() +
  geom_path(data = world, aes(x = long, y = lat, group = group)) +
  # geom_point(data = out, inherit.aes = FALSE, aes(lon, lat, col = mn_cc_lc), size = 1.5, alpha = 0.3) +
  geom_point(data = out2 %>% mutate(rank = rank(mn_cc_rs)), inherit.aes = FALSE, aes(lon, lat, col = rank), size = 1.5, alpha = 0.3) +
  labs(color = "Correlation") +
  # scale_color_gradient2(low = '#C7522B', mid = '#FBF2C4', high = '#3C5941',
  #                      #earth
  #                      # scale_fill_gradient2(low = '#A36B2B', mid = '#EDEAC2', high = '#2686A0',
  #                      # scale_fill_gradient2(low = '#2686A0', mid = '#EDEAC2', high = '#A36B2B',
  #                      midpoint = 0) +
  scale_color_viridis_c(
    breaks = out2 %>%
      mutate(rank = rank(mn_cc_rs)) %>%
      arrange(mn_cc_rs) %>%
      mutate(cut = cut(mn_cc_rs,
        breaks = c(min(mn_cc_rs), c(0, 0.1, 0.2, 0.3, 0.6), max(mn_cc_rs)),
        include.lowest = T
      )) %>%
      group_by(cut) %>%
      summarise(rank = max(rank) + 0.5) %>%
      slice(-n()) %>%
      pull(rank),
    label = c(0, 0.1, 0.2, 0.3, 0.6)
  ) +
  geom_point(data = tplt %>% distinct(Longitude, Latitude), aes(x = Longitude, y = Latitude), col = "red", pch = 10, cex = 6, lwd = 2) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_minimal() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(legend.position = "bottom") +
  theme(
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(1.5, "cm")
  ) +
  xlim(c(-170, -50)) +
  ylim(c(20, 75)) +
  coord_map("bonne", lat0 = 50)
p_map


# histogram of community synchrony
hist(out2$mn_cc_rs, breaks = 20)
range(out2$mn_cc_rs)
mean(out2$mn_cc_rs)


# human footprint ---------------------------------------------------------

# BBS route 25 miles long - get buffer (40 km) when extracting vars
# source: https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-1993-human-footprint/data-download#close

pts_sp <- sp::SpatialPoints(cbind(out2$lon, out2$lat))
sf_pts <- sf::st_as_sf(pts_sp)
st_crs(sf_pts) <- 4326

hfi <- terra::rast(paste0(path, "/data/wildareas-v3-2009-human-footprint-geotiff/wildareas-v3-2009-human-footprint.tif"))
sf_pts3 <- sf::st_transform(sf_pts, crs = crs(hfi)) %>%
  st_buffer(dist = 40000)
eval3 <- exactextractr::exact_extract(hfi, sf_pts3,
  fun = "mean"
)

# add to previous cc dataframe - HFI values are between 0 and 50
out5 <- dplyr::mutate(out2, HFI = eval3) %>%
  dplyr::filter(HFI <= 50)

plt6 <- ggplot(out5, aes(x = HFI, y = mn_cc_rs, color = factor(BCR))) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  xlab("HFI") +
  ylab("Mean pairwise correlation coefficient") +
  geom_line(stat = "smooth", method = "lm", size = 3, alpha = 0.5)

# pdf('./avian/output/sync_hfi.pdf',width = 6, height = 4)
# print(plt6)
# dev.off()

out5$nsp <- as.numeric(out5$nsp)

if (!file.exists(paste0(path, "data/mcmc.rds"))) {
  ff <- rstanarm::stan_gamm4(mn_cc_rs ~ s(nsp) + HFI,
    data = out5,
    chains = 4,
    cores = 4,
    control = list(adapt_delta = 0.99)
  )
  MCMCvis::MCMCsummary(ff, round = 5)
  write_rds(ff, paste0(path, "data/mcmc.rds"))
} else {
  ff <- read_rds(paste0(path, "data/mcmc.rds"))
}

# check accuracy
out5 <- out5 %>% mutate(pred = predict(ff))
ggplot(out5) +
  geom_point(aes(x = mn_cc_rs, y = pred), alpha = 0.3) +
  theme_classic()

# calculate R^2
cor(out5$mn_cc_rs, out5$pred)^2

# Contour plot
pred_y <- predict(ff, newdata = data.frame(nsp = mean(out5$nsp), HFI = 0:50))

axis_x <- seq(min(out5$nsp), max(out5$nsp), length.out = 50)
axis_y <- seq(min(out5$HFI), max(out5$HFI), length.out = 50)
pro_lm_surface <- expand.grid(nsp = axis_x, HFI = axis_y, KEEP.OUT.ATTRS = F)
pro_lm_surface$sc_p <- predict(ff, newdata = pro_lm_surface)

p_contour <- ggplot() +
  geom_point(data = out5 %>% mutate(rank = rank(mn_cc_rs)), aes(x = HFI, y = nsp, col = rank), alpha = 0.3) +
  # geom_contour_filled(data=pro_lm_surface, aes(x=HFI, y=nsp, z=sc_p))+
  geom_contour(data = pro_lm_surface, aes(x = HFI, y = nsp, z = sc_p)) +
  metR::geom_label_contour(data = pro_lm_surface, aes(x = HFI, y = nsp, z = sc_p), skip = 0) +
  ylab("Number of species") +
  xlab("HFI") +
  scale_color_viridis_c(
    breaks = out5 %>%
      mutate(rank = rank(mn_cc_rs)) %>%
      arrange(mn_cc_rs) %>%
      mutate(cut = cut(mn_cc_rs,
        breaks = c(min(mn_cc_rs), c(0, 0.1, 0.2, 0.3, 0.6), max(mn_cc_rs)),
        include.lowest = T
      )) %>%
      group_by(cut) %>%
      summarise(rank = max(rank) + 0.5) %>%
      slice(-n()) %>%
      pull(rank),
    label = c(0, 0.1, 0.2, 0.3, 0.6)
  ) +
  guides(col = "none") +
  # labs(color = 'Synchrony (rho)') +
  # theme(legend.position="bottom")+
  # theme(legend.key.height= unit(0.5, 'cm'),
  #       legend.key.width= unit(1, 'cm'))+
  theme_classic()
p_contour

# pdf('./avian/output/sync_nsp_hfi_contour.pdf', width = 4, height = 4)
# print(p_contour)
# dev.off()

pdf(paste0(path, "output/ts_map_cont.pdf"), width = 12, height = 8)
grid.arrange(annotate_figure(p_ts, fig.lab = "(a)"),
  annotate_figure(p_map, fig.lab = "(b)"),
  annotate_figure(p_contour, fig.lab = "(c)"),
  layout_matrix = rbind(
    c(1, 1),
    c(2, 3),
    c(2, 3)
  )
)
dev.off()
