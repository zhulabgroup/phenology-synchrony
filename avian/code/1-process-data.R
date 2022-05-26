library(pacman)
#######################
# 1 - process BBS data
#######################


# load packages -----------------------------------------------------------
# JAGS downloaded from https://sourceforge.net/projects/mcmc-jags/files/
p_load(bbsBayes)
p_load(tidyverse)


# get data ----------------------------------------------------------------

#only needs to be run once
# fetch_bbs_data()
path<-"./avian/"
load(paste0(path,"data/bbs_raw_data.RData"))

strat_data <- stratify(by = 'bbs_cws', bbs_data = bbs_data)

#process
bstrat <- strat_data$bird_strat

#dplyr::filter(route_df, RouteName == 'ALABAMA PT 2')
#dplyr::filter(bstrat, AOU == 580, Year == 2010, statenum == 2, Route == 141)

#routes are unique within each state (i.e., need to create new)
RID_df <- bstrat %>%
  dplyr::distinct(statenum, Route) %>%
  dplyr::arrange(statenum, Route) %>%
  dplyr::mutate(RID = 1:NROW(.))

#route info
route_df <- strat_data$route_strat %>%
  dplyr::select(statenum, Route, RouteName,
                Latitude, Longitude) %>%
  dplyr::distinct()

#merge new RID (unique ID for each route, irrespective of state) with full data
bstrat2 <- left_join(bstrat, RID_df, by = c('statenum', 'Route')) %>%
  left_join(route_df, by = c('statenum', 'Route')) %>%
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

#unique RIDs
uRID <- unique(bstrat2$RID)

#create empty df to fill

#RESIDUALS MIGHT BE OFF (for route 5146 at least was 0)
#INCLUDE ALL SPECIES, NOT JUST COMPLETE TIME SERIES (determine what to do about zeros in time series) - can filter after using mean count, etc.
#MEAN WHEN MULTIPLE OBS PER YEAR
#at least 5 species for each RID

p_load(foreach)
p_load(doSNOW)
num_cores <- 50
cl <- makeCluster(num_cores)
registerDoSNOW(cl)
m_data_list<-
  foreach (i = 1:length(uRID),
           .packages = c("tidyverse")) %dopar% {
             #i <- 770
             print(paste0('processing route: ', i, ' of ', length(uRID)))
             #data for single route
             r_data <- dplyr::filter(bstrat2, RID == uRID[i])
             
             #only species seen in all years
             sp_keep1 <- r_data %>%
               dplyr::group_by(AOU) %>%
               dplyr::summarize(nyrs = length(unique(Year))) %>%
               dplyr::ungroup() %>%
               dplyr::filter(nyrs == length(unique(r_data$Year)))
             
             #filtered by 'abundant' species
             r_data2 <- r_data %>%
               dplyr::filter(AOU %in% sp_keep1$AOU)
             
             # sp_keep1 <- data.frame(AOU = unique(r_data$AOU))
             
             #get mean count of each species across years
             sp_keep2 <- r_data2 %>%
               dplyr::group_by(AOU) %>%
               dplyr::summarize(mn_count = mean(SpeciesTotal)) %>%
               dplyr::ungroup() %>%
               dplyr::arrange(AOU)
             
             r_data3 <- dplyr::left_join(r_data2, sp_keep2, by = 'AOU') %>%
               dplyr::arrange(AOU, Year)
             # r_data2 <- r_data
             
             
             #get residuals from linear regression (on log count)
             resids <- r_data3 %>%
               group_by(AOU) %>%
               do(model = residuals(lm(log_count ~ Year, data = .)))
             
             #add to df
             r_data3$residuals <- unlist(resids$model)
             
             if (NROW(r_data3) >= 5)
             {
               r_data3
             }
           }
m_data<-bind_rows(m_data_list)

# save RDS ----------------------------------------------------------------
year_list<-m_data %>% pull(Year) %>% unique()
path_out<-paste0(path, "data/processed/")
dir.create(path_out, recursive = T)
for (year in year_list) {
  write_rds(m_data %>% filter(Year==year), paste0(path_out, year,'.rds'))
  print(year)
}


# plot time series for single site ----------------------------------------

#line size
SZ <- 2
#text size
TSZ <- 14
#text title size
TTSZ <- 16
#bounding box size
LSZ <- 1.5
#tick size
TISZ <- 1.5

tplt <- dplyr::filter(m_data, RID == uRID[1])
#plot log count ~ year
lc1 <- ggplot(tplt, aes(Year, log_count, color = factor(AOU))) +
  #geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  theme_bw() +
  theme(legend.pos = 'none') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = LSZ),
        axis.ticks = element_line(size = TISZ),
        axis.text.x = element_text(size = TSZ),
        axis.title.x = element_text(size = TTSZ),
        axis.text.y = element_text(size = TSZ),
        axis.title.y = element_text(size = TTSZ),
        legend.position = 'none') +
  xlab('Year') +
  ylab('log(Count)') +
  ggtitle(paste0(tplt$RouteName[1], ' (', round(tplt$Latitude[1], 2), ', ', round(tplt$Longitude[1], 2), ')'))

#plot residuals of log count ~ year
lc2 <- ggplot(tplt, aes(Year, residuals, color = factor(AOU))) +
  #geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  theme_bw() +
  theme(legend.pos = 'none') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = LSZ),
        axis.ticks = element_line(size = TISZ),
        axis.text.x = element_text(size = TSZ),
        axis.title.x = element_text(size = TTSZ),
        axis.text.y = element_text(size = TSZ),
        axis.title.y = element_text(size = TTSZ),
        legend.position = 'none') +
  xlab('Year') +
  ylab('log(Count) residuals') +
  ggtitle(paste0(tplt$RouteName[1], ' (', round(tplt$Latitude[1], 2), ', ', round(tplt$Longitude[1], 2), ')'))

#save figures
dir.create("./avian/output", recursive = T)
pdf( './avian/output/lc_year.pdf', width = 6, height = 4)
print(lc1)
dev.off()

pdf( './avian/output/resid_year.pdf', width = 6, height = 4)
print(lc2)
dev.off()


