library(tidyverse)
library(corrplot)
library(ggrepel)
library(raster)
library(minpack.lm)

data<-read_csv("./marine/data/bethany_bottomtrawl_fall.csv") %>% 
  dplyr::select(-SEASON, -name, -LAT, -LON) %>% 
  gather(key="species", value="CPUE", -YEAR, -STRATUM, -NTOWS, -midlat, -midlon, -group) %>% 
  mutate(abundance=CPUE*NTOWS) %>% 
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

require(maps)
require(sp)
require(maptools)
usa <- map("state", fill = TRUE)
IDs <- sapply(strsplit(usa$names, ":"), function(x) x[1])
usa <- map2SpatialPolygons(usa, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

area<-extent(min(data$midlon)-0.5,max(data$midlon)+0.5, min(data$midlat)-0.5,max(data$midlat)+0.5 ) 
usa_crop <- crop(usa,area)

p<-ggplot()+
  geom_polygon( data=usa_crop, aes(x=long, y=lat, group=group),
                color="darkblue", fill="lightblue", size = .1 )+
  geom_point(data=data %>% 
               distinct(midlon, midlat, id),
             aes(x=midlon, y=midlat))+
  geom_label_repel(data=data %>% 
              distinct(midlon, midlat, id),
            aes(x=midlon, y=midlat, label=id))+
  theme_minimal()+
  xlab("longitude")+
  ylab("latitude")+
  guides(col=guide_legend(title=""))+
  coord_equal()
  
cairo_pdf("./marine/output/map.pdf", width = 8, height = 8)
print(p)
dev.off()
# test<-data %>% filter(species==unique(data$species)[5]) %>% 
#   pull(abundance)  
# mean(test==0) # check with Casey about 0

sp_list<-unique(data$species)
corr_decay_list<-mean_corr_list<-rep(NA, length(sp_list))
for (i in 1:length(sp_list)) {
  sp<-sp_list[i]
  data_sp<-data %>%
    filter(species==sp_list[i]) %>% 
    arrange(midlat, midlon) %>% 
    dplyr::select(-STRATUM, -NTOWS, -group,-midlat, -midlon, -species, -CPUE) %>% 
    group_by(site,id) %>% 
    filter(sum(abundance)>1000) %>% 
    ungroup() 
  
  dir.create("./marine/archive/ts/", recursive = T)
  cairo_pdf(paste0("./marine/archive/ts/",sp, ".pdf"))
  p<-ggplot(data_sp)+
    geom_line(aes(x=YEAR, y=logabun, col=site, group=site))+
    guides(col="none")+
    theme_classic()+
    xlab("year")+
    ylab("log (abundance + 1)")+
    ggtitle(sp)
  print(p)
  dev.off()
  
  data_mat<- data_sp %>% 
    dplyr::select(-site, -abundance) %>% 
    spread(key = "id", value="logabun") %>% 
    dplyr::select(-YEAR)
  
  res<-cor(data_mat,use="complete.obs")
  # variogram
  
  dir.create("./marine/archive/corr_mat/", recursive = T)
  cairo_pdf(paste0("./marine/archive/corr_mat/",sp, ".pdf"))
  corrplot(res, method="circle", title=sp)
  dev.off()
  
  mean_corr_list[i]<-mean(res-diag(nrow=nrow(res), ncol=ncol(res)), na.rm = T)
  
  corr_df<-as.data.frame(res) %>% 
    rownames_to_column(var="start") %>% 
    gather(key="end", value = "correlation", -start) %>%
    mutate(start=as.integer(start),
           end=as.integer(end)) %>% 
    left_join(coord_df, by=c("start"="id")) %>% 
    rename(start_lat=midlat,
           start_lon=midlon) %>% 
    left_join(coord_df, by=c("end"="id")) %>% 
    rename(end_lat=midlat,
           end_lon=midlon) %>% 
    mutate(distance=sqrt((start_lat-end_lat)^2+(start_lon-end_lon)^2)) 
  
  dir.create("./marine/archive/corr_map/", recursive = T)
  cairo_pdf(paste0("./marine/archive/corr_map/",sp, ".pdf"))
  p<-ggplot ()+
    geom_polygon( data=usa_crop, aes(x=long, y=lat, group=group),
                  color="darkblue", fill="lightblue", size = .1 )+
    geom_point(data=data,aes(x=midlon, y=midlat))+
    geom_segment(data=corr_df,aes(x=start_lon, y=start_lat, xend=end_lon, yend=end_lat, col=correlation))+
    scale_color_viridis_c(begin=0, end=1)+
    theme_minimal()+
    xlab("longitude")+
    ylab("latitude")+
    coord_equal()+
    ggtitle(sp)
  
  print(p)
  dev.off()
  
  model <- nlsLM(correlation ~ SSasymp(distance, yf, y0, log_alpha), data = corr_df ,
                 start=c(y0=1, yf=0, log_alpha=log(1)),
                 lower = c(0, 0, 0),
                 upper=c(1,1, 12))
  corr_df <- corr_df %>% 
    mutate(fit = predict(model))  
  
  y0<-coef(model)[1]
  yf<-coef(model)[2]
  phi<- exp(coef(model)[3]) 
  range<- -log(0.05)/phi
  
  dir.create("./marine/archive/decay/", recursive = T)
  cairo_pdf(paste0("./marine/archive/decay/",sp, ".pdf"))
  p<-ggplot(corr_df)+
    geom_point(aes(x=distance, y=correlation))+
    geom_line(aes(x=distance, y=fit), color="blue")+
    geom_smooth(aes(x=distance, y=correlation),method=loess,col="red")+
    # geom_smooth(aes(x=distance, y=correlation),method=lm, col="red")+
    geom_hline(yintercept = yf, lty=2)+
    geom_hline(yintercept = y0, lty=2)+
    geom_vline(xintercept = range, lty=2)+
    xlim(0,12)+
    ylim(0,1)+
    theme_classic()+
    ggtitle(sp)
  print(p)
  dev.off()
  
  # corr_decay_list[i]<-coef(lm(data=corr_df, correlation~distance))["distance"]
  corr_decay_list[i]<-range
  
  print(paste0(i, ", ", sp))
}

management=c(1, 1, 1, 1, 0,
          0, 0, 1, 1, 1,
          1, 0, 1, 0, 1,
          1)
# https://media.fisheries.noaa.gov/2021-04/Mid-Atlantic-Managed-Species.pdf
# https://media.fisheries.noaa.gov/2021-04/New-England-Managed-Species.pdf
# https://www.fisheries.noaa.gov/new-england-mid-atlantic/population-assessments/fishery-stock-assessments-new-england-and-mid-atlantic

sync_df<-data.frame(species=sp_list, corr=mean_corr_list,decay=corr_decay_list,
                    management=management) %>% 
  mutate(management=case_when(management==0~"unmanaged",
                           management==1~"managed")) %>% 
  arrange(species)

cairo_pdf("./marine/output/mean_corr.pdf")
p<-ggplot(sync_df)+
  geom_boxplot(aes(x=management, y=corr))+
  geom_point(aes(x=management, y=corr), cex=2, col="red", pch=1)+
  geom_label_repel(aes(x=management, y=corr, label=species), cex=3, col="red")+
  theme_classic()+
  xlab("")+
  ylab("mean correlation")
print(p)
dev.off()

mean(sync_df %>% filter(management=="managed") %>% pull(corr))
sd(sync_df %>% filter(management=="managed") %>% pull(corr))

mean(sync_df %>% filter(management=="unmanaged") %>% pull(corr))
sd(sync_df %>% filter(management=="unmanaged") %>% pull(corr))

t.test(sync_df %>% filter(management=="managed") %>% pull(corr),
       sync_df %>% filter(management=="unmanaged") %>% pull(corr))


cairo_pdf("./marine/output/corr_decay.pdf")
p<-ggplot(sync_df)+
  geom_boxplot(aes(x=management, y=decay))+
  geom_point(aes(x=management, y=decay), cex=2, col="red", pch=1)+
  geom_label_repel(aes(x=management, y=decay, label=species), cex=3, col="red")+
  theme_classic()+
  xlab("")+
  ylab("log (effective range)")
print(p)
dev.off()

mean(sync_df %>% filter(management=="managed") %>% pull(decay))
sd(sync_df %>% filter(management=="managed") %>% pull(decay))

mean(sync_df %>% filter(management=="unmanaged") %>% pull(decay))
sd(sync_df %>% filter(management=="unmanaged") %>% pull(decay))

t.test(sync_df %>% filter(management=="managed") %>% pull(decay),
       sync_df %>% filter(management=="unmanaged") %>% pull(decay))


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
  mutate(management=case_when(management==0~"unmanaged",
                              management==1~"managed")) %>% 
  arrange(species)

cairo_pdf("./marine/output/abundance_roc.pdf")
p<-ggplot(abundance_df)+
  geom_boxplot(aes(x=management, y=roc))+
  geom_point(aes(x=management, y=roc), cex=2, col="red", pch=1)+
  geom_label_repel(aes(x=management, y=roc, label=species), cex=3, col="red")+
  theme_classic()+
  xlab("")+
  ylab("abundance rate of change")
print(p)
dev.off()

t.test(abundance_df %>% filter(management=="managed") %>% pull(roc),
       abundance_df %>% filter(management=="unmanaged") %>% pull(roc))

# does depth matter