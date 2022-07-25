library(pacman)
p_load(tidyverse)
p_load(ggcorrplot)
p_load(ggrepel)
p_load(raster)
p_load(maps)
p_load(minpack.lm)
p_load(sp)
p_load(maptools)
p_load(rgeos)
p_load(gridExtra)
p_load(geosphere)
p_load(tsibble)
p_load(fable)

path<-"./marine/"
data<-read_csv(paste0(path,"data/bethany_bottomtrawl_fall.csv")) %>% 
  dplyr::select(-SEASON, -name, -LAT, -LON) %>% 
  tidyr::gather(key="species", value="CPUE", -YEAR, -STRATUM, -NTOWS, -midlat, -midlon, -group) %>% 
  mutate(abundance=abs(CPUE)) %>% 
  mutate(logabun=log(abundance+1)) %>% 
  mutate(site=paste0(midlon,"_", midlat))
  # group_by (YEAR, group, species) %>% 
  # summarize (abundance=sum(abundance),
  #            lat=mean(midlat),
  #            lon=mean(midlon)) %>% 
  # ungroup()

coord_df<-data %>% 
  distinct(midlon, midlat) %>% 
  arrange(midlat, midlon) %>%
  mutate(id=row_number())

data<-data %>% 
  left_join(coord_df, by=c("midlat", "midlon"))


usa <- map("state", fill = TRUE)
IDs <- sapply(strsplit(usa$names, ":"), function(x) x[1])
usa <- map2SpatialPolygons(usa, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

area<-extent(min(data$midlon)-0.5,max(data$midlon)+0.5, min(data$midlat)-0.5,max(data$midlat)+0.5 ) 
usa_crop <- crop(usa,area)

p_map<-ggplot()+
  geom_polygon( data=usa_crop, aes(x=long, y=lat, group=group),
                color="darkblue", fill="lightblue", size = .1 )+
  geom_point(data=data %>% 
               distinct(midlon, midlat, id),
             aes(x=midlon, y=midlat))+
  # geom_label_repel(data=data %>% 
  #             distinct(midlon, midlat, id),
  #           aes(x=midlon, y=midlat, label=id))+
  theme_minimal()+
  xlab("Longitude")+
  ylab("Latitude")+
  guides(col=guide_legend(title=""))+
  coord_equal()
p_map
pdf("./marine/output/map.pdf", width = 8, height = 8)
print(p_map)
dev.off()

sp_list<-unique(data$species) %>% sort()
decay_df_list<-vector(mode="list", length(sp_list))
for (i in 1:length(sp_list)) {
  sp<-sp_list[i]
  data_sp<-data %>%
    filter(species==sp_list[i]) %>% 
    arrange(midlat, midlon) %>% 
    dplyr::select(-STRATUM, -NTOWS, -group,-midlat, -midlon, -species, -CPUE) %>% 
    group_by(site,id) %>% 
    filter(sum(abundance==0)/n()<0.50) %>% 
    ungroup() 
  
  data_tsbl<-data_sp %>% 
    dplyr::select(-site) %>% 
    as_tsibble(index = YEAR, key=id)
  lm_df <- data_tsbl %>%
    model(TSLM(logabun ~ YEAR)) %>%
    broom::tidy(fit) %>% 
    dplyr::select(id, term, estimate) %>% 
    spread(key="term", value="estimate") %>% 
    rename(intercept=`(Intercept)`, beta=YEAR)
  
  data_sp<-data_sp %>% 
    left_join(lm_df , by=c("id")) %>% 
    mutate(predict=YEAR*beta+intercept) %>% 
    mutate(detrend=logabun-predict)
  
  # data_sp %>% 
  #   ggplot(aes(x=YEAR, y=detrend, col=id %>% as.factor())) +
  #   geom_line() 
  
  see_sites<-data_sp %>% 
    group_by(site, id) %>%
    summarise(sum=sum(abundance)) %>%
    arrange(desc(sum)) %>% 
    head(6) %>% 
    pull(id)
  p_ts<-ggplot(data_sp %>% filter(id %in% see_sites))+
    geom_line(aes(x=YEAR, y=logabun, col=site, group=site), alpha=0.5)+
    guides(col="none")+
    theme_classic()+
    xlab("Year")+
    ylab("log (CPUE + 1)")
    if (i==8) {
    p_ts_1sp<-p_ts
    }
  
  data_mat<- data_sp %>% 
    dplyr::select(YEAR, id, detrend) %>% 
    spread(key = "id", value="detrend") %>% 
    dplyr::select(-YEAR)
  
  res<-cor(data_mat,use="complete.obs")
  res[lower.tri(res)]<-NA

    # variogram
  
  p_corrmat<-ggcorrplot(res, method="circle")
  
  mean_corr<-mean(res-diag(nrow=nrow(res), ncol=ncol(res)), na.rm = T)
  
  corr_df<-as.data.frame(res) %>% 
    rownames_to_column(var="start") %>% 
    tidyr::gather(key="end", value = "correlation", -start) %>%
    mutate(start=as.integer(start),
           end=as.integer(end)) %>% 
    left_join(coord_df, by=c("start"="id")) %>% 
    rename(start_lat=midlat,
           start_lon=midlon) %>% 
    left_join(coord_df, by=c("end"="id")) %>% 
    rename(end_lat=midlat,
           end_lon=midlon) %>% 
    rowwise() %>% 
    mutate(distance=distm(c(start_lon, start_lat), c(end_lon, end_lat), fun = distHaversine) %>% as.numeric() %>% `/`(100000)
           # distance=sqrt((start_lat-end_lat)^2+(start_lon-end_lon)^2)
           ) %>% 
    ungroup() %>% 
    filter(distance>0) %>% 
    drop_na(correlation)
  
  # set.seed(42)
  p_corrmap<-ggplot ()+
    geom_polygon( data=usa_crop, aes(x=long, y=lat, group=group),
                  color="darkblue", fill="lightblue", size = .1 )+
    geom_point(data=data ,aes(x=midlon, y=midlat))+
    geom_segment(data=corr_df %>%
                   filter(distance<=quantile(corr_df$distance, 0.5))
                 ,aes(x=start_lon, y=start_lat, xend=end_lon, yend=end_lat, alpha=abs(correlation), col=correlation))+
    scale_color_viridis_c()+
    guides(alpha="none",
           col=guide_legend(title="Correlation"))+
    theme_minimal()+
    xlab("Longitude")+
    ylab("Latitude")+
    coord_equal()
  
  if (i==8) {
    p_corrmap_1sp<-p_corrmap
  }
  
  loess_model<-loess(correlation ~ distance, data=corr_df)
  smooth<-predict(loess_model) 
  corr_df<-corr_df %>% 
    mutate(smooth=smooth)
  
  model <- nlsLM(correlation ~ SSasymp(distance, yf, y0, log_alpha), data = corr_df ,
                 start=c(y0=0.5, yf=0, log_alpha=log(-log(0.5)/5)),
                 lower = c(0, 0, log(-log(0.5)/10)),
                 upper=c(1,1, log(-log(0.5)/0.1)),
                 trace = F)
  corr_df <- corr_df %>% 
    mutate(fit = predict(model))  
  
  upper<-coef(model)[1]
  lower<-coef(model)[2]
  phi<- exp(coef(model)[3]) 
  range<- -log(0.5)/phi
  
  p_decay<-ggplot(corr_df %>% mutate(bin=distance%/%0.5*0.5*100 ) )+
    geom_boxplot(aes(x=bin, y=correlation, group=bin), alpha=0.25)+
    # geom_point(aes(x=distance*100, y=correlation), alpha=0.25)+
    geom_line(aes(x=distance*100, y=fit), color="blue", lwd=2)+
    # geom_smooth(aes(x=distance, y=correlation),method=loess,col="red")+
    # geom_line(aes(x=distance, y=smooth), col="red")+
    geom_hline(yintercept = lower, lty=2, col="red")+
    geom_hline(yintercept = upper, lty=2, col="red")+
    geom_vline(xintercept = range*100, lty=2, col="red")+
    xlim(-50,12*100+50)+
    # ylim(0,1)+
    ylab ("Pairwise correlation coefficient")+
    xlab ("Distance (km)")+
    theme_classic()
  # p_decay
  if (i==8) {
    p_decay_1sp<-p_decay
  }
  
  decay_df_list[[i]]<-tibble(mean=mean_corr, 
                                 upper=upper, 
                                 lower=lower,
                                 range=range,
                                 n=nrow(corr_df))
  
  dir.create(paste0(path, "output/by_species/"), showWarnings = F)
  pdf(paste0(path,"output/by_species/",sp,".pdf"), width = 16, height = 16)
  grid.arrange(p_ts, p_corrmat, p_corrmap, p_decay, nrow=2, top=sp)
  dev.off()
  
  print(paste0(i, ", ", sp))
}
decay_df<-bind_rows(decay_df_list)

p_ts_1sp
p_corrmap_1sp
p_decay_1sp

sp_list
management=c(0, 1, 0, 1, 0,
          1, 1, 1, 1, 1,
          0, 1, 1, 1, 1,
          0)
# https://media.fisheries.noaa.gov/2021-04/Mid-Atlantic-Managed-Species.pdf
# https://media.fisheries.noaa.gov/2021-04/New-England-Managed-Species.pdf
# https://www.fisheries.noaa.gov/new-england-mid-atlantic/population-assessments/fishery-stock-assessments-new-england-and-mid-atlantic

sync_df<-tibble(species=sp_list, decay_df,
                    management=management) %>% 
  mutate(management=case_when(management==0~"Unmanaged",
                           management==1~"Managed")) %>% 
  arrange(species)



abundance_df<-data %>% 
  group_by(species, YEAR) %>% 
  summarize(abundance=sum(abundance)) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  do(broom::tidy(lm(abundance ~ YEAR, .))) %>%
  filter(term == "YEAR") %>%
  dplyr::select(species, roc = estimate, p = p.value) %>% 
  ungroup() %>% 
  mutate(management=management) %>% 
  mutate(management=case_when(management==0~"Unmanaged",
                              management==1~"Managed")) %>% 
  arrange(species)

summary_df<-full_join(sync_df, abundance_df %>% dplyr::select(-p),by=c("species", "management"))
pdf(paste0(path, "output/summary.pdf"))
p_summary<-ggplot(summary_df %>%
                    mutate(range=range*100) %>% 
                    dplyr::select(species, management, 
                                  # `mean correlation`=mean,
                                  `Distance where correlation decays to half (km)`=range,
                                  `Rate of change in CPUE (per year)`=roc) %>% 
                    tidyr::gather(key="var", value="value", -species, -management)  %>% 
                    mutate(var=as.factor(var)))+
  geom_boxplot(aes(x=management, y=value))+
  geom_point(aes(x=management, y=value), cex=2, col="red", pch=1)+
  # geom_label_repel(aes(x=management, y=value, label=species), cex=3, col="red")+
  theme_classic()+
  xlab("")+
  facet_wrap(.~var, scales = "free_y", labeller = labeller(var = label_wrap_gen(30)))+
  ylab("")
print(p_summary)
dev.off()

pdf(paste0(path, 'output/ts_map_decay_box.pdf'), width = 12, height = 10)
grid.arrange(annotate_figure(p_ts_1sp, fig.lab = "(a)"),
             annotate_figure(p_corrmap_1sp, fig.lab = "(b)"),
             annotate_figure(p_decay_1sp, fig.lab = "(c)"),
             annotate_figure(p_summary, fig.lab = "(d)"),
             layout_matrix=rbind(c(1,2),
                                 c(3, 4))
)
dev.off()

summary_df %>%
  tidyr::gather(key="var", value="value", -species, -management) %>% 
  filter(var!="number of pairs of sites") %>%  
  group_by(management, var) %>% 
  summarise(mean=mean(value),
            sd=sd(value))

t.test(sync_df %>% filter(management=="Managed") %>% pull(range),
       sync_df %>% filter(management=="Unmanaged") %>% pull(range),
       alternative = "less")

t.test(sync_df %>% filter(management=="Managed") %>% pull(mean),
       sync_df %>% filter(management=="Unmanaged") %>% pull(mean),
       alternative = "greater")

t.test(abundance_df %>% filter(management=="Managed") %>% pull(roc),
       abundance_df %>% filter(management=="Unmanaged") %>% pull(roc),
       alternative = "greater")
