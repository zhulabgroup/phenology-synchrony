library(pacman)
p_load(tidyverse)

path <- "./mosquito/"
pred.df <- readRDS(paste0(path, "data/predict.rds")) %>% 
  filter(Site !="HARV")
pred1.df <- readRDS(paste0(path, "data/predict_year.rds"))
pred3.df <- readRDS(paste0(path, "data/predict_long.rds")) %>% 
  filter(Site != "HARV")

full.df <- bind_rows(pred.df, pred1.df,pred3.df)

usa <- map("state", fill = TRUE)
IDs <- sapply(strsplit(usa$names, ":"), function(x) x[1])
usa <- map2SpatialPolygons(usa, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
p_map<-ggplot()+
  geom_polygon( data=usa, aes(x=long, y=lat, group=group),
                color="darkblue", fill="lightblue", size = .1 )+
  geom_point(data=site_df ,
             aes(x=lon, y=lat))+
  geom_label_repel(data=site_df,
            aes(x=lon, y=lat, label=Site))+
  theme_minimal()+
  xlab("longitude")+
  ylab("latitude")+
  coord_equal()

count.df <- full.df %>%
  group_by(Site, Year, SciName, DOY) %>% 
  summarise(
    Count = mean(Count)#,
    # Q_2.5 = mean(Q_2.5),
    # Q_97.5 = mean(Q_97.5),
  ) %>%
  ungroup() %>% 
  group_by(Site, Year, SciName) %>% 
  mutate(cumCount=cumsum(Count)) %>% 
  mutate(proCount=cumCount/max(cumCount))

count.df %>% 
  ggplot( aes(x=DOY,y=proCount, color=interaction(SciName, Year))) + 
  geom_line(size=2,alpha=.5)+
  facet_wrap(.~Site)+
  theme(legend.position = "none")

p_ts<-ggplot()+
  geom_line(data=count.df %>% filter(Site=="HARV", SciName=="Aedes trivittatus"),
            aes(x=DOY,y=cumCount, col=Year))+
  geom_rect(data=range.df %>% filter(Site=="HARV", SciName=="Aedes trivittatus"),
            aes(xmin=start, xmax=end,ymin=0, ymax=Inf, fill=Year),alpha=0.25)+
  # geom_vline(data=range.df %>% filter(Site=="HARV", SciName=="Aedes trivittatus"),
  #           aes(xintercept=start, col=Year))+
  # geom_vline(data=range.df %>% filter(Site=="HARV", SciName=="Aedes trivittatus"),
  #            aes(xintercept=end, col=Year))+
  theme_classic()+
  facet_wrap(.~Year, ncol=1, scales="free_y")+
  guides(fill="none")+
  xlab("day of year")+
  ylab ("cumulated abundance")+
  ggtitle("HARV")
  
## generating the number of days between the 10% quantile nad 90th quanitile
## Higher values indicate

range.df <- count.df %>%
  filter(proCount >= .1 & proCount <= .9) %>% 
  group_by( SciName, Year, Site) %>%
  summarise(start=min(DOY),
            end=max(DOY),
            TDays = n())

# ggplot(range.df , aes(x=SciName, y= TDays, col=Year)) + geom_point() + 
#   facet_wrap(.~Site)

ggplot(range.df %>%
         # filter(Site != "ORNL" & Site != "SERC" & Site != "TALL" & Site != "WREF" & Site != "YELL") %>% 
         left_join( range.df %>% group_by(SciName, Site) %>% distinct(Year) %>% summarise(n=n()),
                    by=c("Site", "SciName")) %>% 
         filter(n>3) %>% 
         mutate(genus=str_split(SciName," ", simplify = T)[1]),
       aes(x=Year, y= TDays, col=SciName,fill=SciName,group=SciName)) +
  geom_line() + 
  # geom_smooth(method="lm", se=F, alpha=0.3)+
  facet_wrap(.~Site, ncol=1)+
  theme_classic()+
  guides(col="none")

range.df.proc<-range.df %>%
  left_join( range.df %>% group_by(SciName) %>% distinct(Year, Site) %>% summarise(n=n()),
             by=c( "SciName")) %>% 
  filter(n>=10) %>%
  mutate(Year=as.numeric(Year)) %>% 
  left_join(daymet_df, by=c("Site", "Year"))

p_corr<-ggplot(range.df.proc)+
  geom_point(aes(x=mat, y = TDays, col=Site))+
  # geom_smooth(aes(x=mat, y = TDays, col=Site),method="lm", se=F)+
  geom_smooth(aes(x=mat, y = TDays),method="lm", se=T)+
  facet_wrap(.~SciName, scales="free_y")+
  theme_classic()+
  xlab("mean annual temperature (°C)")+
  ylab ("TDays (day)")

p_load(nlme)
reg_df_list<-vector(mode="list")
for (sp in range.df.proc %>% pull(SciName) %>% unique()) {
  lme.fit <- lme(TDays~mat,random=~1|Site,data=range.df.proc %>% filter(SciName==sp))
  res<-summary(lme.fit)
  reg_df_list[[sp]]<-data_frame(SciName=sp,
             estimate=res$tTable[2, 1],
             se=res$tTable[2, 2],
             p=res$tTable[2, 5])
}
reg_df<-bind_rows(reg_df_list)

p_summary<-ggplot()+
  geom_point(data=reg_df, aes(x=SciName, y=estimate))+
  geom_point(data=reg_df%>% filter(p<0.05), aes(x=SciName, y=estimate), col="blue")+
  geom_errorbar(data=reg_df, aes(x=SciName, ymin=estimate-1.96*se, ymax=estimate+1.96*se))+
  geom_errorbar(data=reg_df %>% filter(p<0.05), aes(x=SciName, ymin=estimate-1.96*se, ymax=estimate+1.96*se), col="blue")+
  geom_hline(yintercept = 0, lty=2)+
  theme_classic()+
  xlab("species")+
  ylab ("regression coefficient (day / °C)")+
  coord_flip()+
  scale_x_discrete(limits=rev)
summary(reg_df$estimate>0)
summary(reg_df$p<0.05)

pdf('./mosquito/output/map_ts_corr_sum.pdf', width = 12, height = 8)
grid.arrange(annotate_figure(p_map, fig.lab = "A"),
             annotate_figure(p_ts, fig.lab = "B"),
             annotate_figure(p_corr, fig.lab = "C"),
             annotate_figure(p_summary, fig.lab = "D"),
             layout_matrix=rbind(c(1,2),
                                 c(3, 4))
)
dev.off()
