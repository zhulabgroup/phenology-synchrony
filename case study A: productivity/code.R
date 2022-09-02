library(lubridate)
path<-"./case study A: productivity/"

# read in flux data
flux_df <- bind_rows(read_csv(paste0(path,"data/FLX_US-Ne2_FLUXNET2015_FULLSET_2001-2013_1-4/FLX_US-Ne2_FLUXNET2015_FULLSET_DD_2001-2013_1-4.csv")) %>% 
                       mutate(site="IR"),
                     read_csv(paste0(path,"data/FLX_US-Ne3_FLUXNET2015_FULLSET_2001-2013_1-4/FLX_US-Ne3_FLUXNET2015_FULLSET_DD_2001-2013_1-4.csv")) %>% 
                                mutate(site="RF"))%>% 
  replace(. == "-9999", NA) %>% 
  rename(date=TIMESTAMP) %>% 
  dplyr::select(date, site, GPP_NT_CUT_50, TA_ERA_DAY, P_F) %>%
  mutate(date=ymd(date))  %>%
  filter(date > as.Date("2000-12-31"),
         date < as.Date("2013-01-01")) 

# read in daymet data 
daymet_df <- bind_rows(read_csv(paste0(path,"data/NE_211742_lat_41.1649_lon_-96.4701_2022-02-02_122643.csv"), skip=7) %>% 
                         mutate(site="IR"),
                       read_csv(paste0(path,"data/NE_311742_lat_41.1797_lon_-96.4397_2022-02-02_122805.csv"), skip=7) %>% 
                         mutate(site="RF")) %>% 
  mutate(date= as.Date(paste(year, yday, sep="-"), format="%Y-%j"))%>%
  rename(Daymet_Precip=`prcp (mm/day)`,Daymet_Tmax=`tmax (deg c)`, Daymet_Tmin=`tmin (deg c)`) %>%
  dplyr::select(date, Daymet_Precip, Daymet_Tmax, Daymet_Tmin, site)

daily_df<-left_join(flux_df, daymet_df, by=c("date", "site")) %>% 
  mutate(month=format(date, "%m") %>% as.integer(),
         year=format(date, "%Y") %>% as.integer())

yearly_df<- daily_df%>% 
  filter(month > 4  & month < 10)%>% # growing season
  group_by(year, site) %>%
  summarise(GPP=sum(GPP_NT_CUT_50, na.rm=TRUE), temp_max=max(TA_ERA_DAY), temp_min=min(TA_ERA_DAY), tower_precip=sum(P_F, na.rm=TRUE), daymet_precip=sum(Daymet_Precip, na.rm=TRUE), daymet_tmax=mean(Daymet_Tmax, na.rm=TRUE), daymet_tmin=mean(Daymet_Tmin, na.rm=TRUE))

pdf(paste0(path, 'output/precip.pdf'), width = 8, height = 4)
ggplot(df %>% 
         dplyr::select(year,site, FLUXNET=tower_precip,Daymet=daymet_precip) %>% 
         tidyr::gather(key="data", value="value", -year, -site) )+
  geom_line(aes(x=year, y=value, group=site, col=site, lty=site), lwd=1)+
  scale_color_brewer(palette = "Set1")+
  facet_wrap(.~data, nrow=1)+
  theme_classic()+
  ylab ("Total annual precipitation (mm)")+
  xlab("Year")+
  guides(col=guide_legend(title="Site"),
         lty="none")+
  scale_x_continuous(breaks= scales::pretty_breaks())
dev.off()

# read in EVI data
evi_df <- read_csv(paste0(path, "data/Nebraska-All-MOD13A1-006-results.csv")) %>%
  dplyr::select(site=ID, date=Date, evi="MOD13A1_006__500m_16_days_EVI")%>%
  mutate(date=as.Date(date))%>%
  filter(date >= as.Date("2001-01-01")) %>% 
  mutate(site=case_when(site=="US-NE2"~ "IR",
                        site=="US-NE3"~"RF")) %>% 
  mutate(year=format(date, "%Y") %>% as.integer())

metric_list<-c("SOS", "EOS", "LOS", "MGS","Peak")
metric_df <-evi_df %>% 
  dplyr::select(-date, -year) %>% 
  group_by(site) %>% 
  summarise(ts(evi, frequency = 23, start = c(2001,1) %>% as.numeric()) %>% 
           FillPermanentGaps( min.gapfrac = 0.2, fill=0.15) %>% 
           as.ts() %>% 
           fortify() %>% 
           # filter(Index < 2013) %>%
           dplyr::select(-Index) %>% 
           ts(start=c(2001,1) %>% as.numeric(), frequency=23) %>% 
           Phenology(tsgf="TSGFspline", approach="White") %>% 
           magrittr::extract(4:15) %>% 
           as.data.frame() %>% 
           tidyr::gather(key="metric_low", value="value") %>% 
           group_by(metric_low) %>% 
           summarise(mean=mean(value, na.rm=T),
                     sd=sd(value, na.rm=T)) ) %>% 
  left_join(data.frame(metric=metric_list, metric_low=tolower(metric_list)), by="metric_low") %>% 
  mutate(metric=as.factor(metric)) %>% 
  mutate(metric=fct_relevel(metric, levels=metric_list)) %>% 
  drop_na(metric)
  
p_metric<-ggplot(metric_df)+
  geom_point(aes(x=site, y=mean))+
  geom_errorbar(aes(x=site, ymin=mean-1.95*sd, ymax=mean+1.95*sd), width=0.3)+
  facet_wrap(.~metric, nrow=1, scales = "free_y")+
  theme_classic()+
  ylab("Value of phenological metric")+
  xlab("Site")
p_metric


evi_annual_df <-evi_df %>% 
  dplyr::select(-date, -year) %>% 
  group_by(site) %>% 
  summarise(ts(evi, frequency = 23, start = c(2001,1) %>% as.numeric()) %>% 
              FillPermanentGaps( min.gapfrac = 0.2, fill=0.15) %>% 
              as.ts() %>% 
              fortify() %>% 
              # filter(Index < 2013) %>%
              dplyr::select(-Index) %>% 
              ts(start=c(2001,1) %>% as.numeric(), frequency=23) %>% 
              Phenology(tsgf="TSGFspline", approach="White") %>% 
              magrittr::extract(4:15) %>% 
              as.data.frame() %>% 
              tidyr::gather(key="metric_low", value="value") %>% 
              group_by(metric_low) %>% 
              mutate(year=2001+row_number()-1) %>% 
              ungroup()) %>% 
  left_join(data.frame(metric=metric_list, metric_low=tolower(metric_list)), by="metric_low") %>% 
  mutate(metric=as.factor(metric)) %>% 
  mutate(metric=fct_relevel(metric, levels=metric_list)) %>% 
  drop_na(metric) %>% 
  left_join(daymet_df %>% 
              mutate(month=format(date, "%m") %>% as.integer(),
                     year=format(date, "%Y") %>% as.integer()) %>% 
              group_by(year, site) %>%
              summarise( daymet_precip=sum(Daymet_Precip, na.rm=TRUE), daymet_tmax=mean(Daymet_Tmax, na.rm=TRUE), daymet_tmin=mean(Daymet_Tmin, na.rm=TRUE))
              )

reg_df<-evi_annual_df %>% 
  group_by(site, metric) %>%
  do(broom::tidy(lm(value ~ daymet_precip, .))) %>%
  filter(term %in% c("daymet_precip")) %>%
  dplyr::select( -statistic)


p_reg<-ggplot(evi_annual_df %>% 
                filter(metric %in% c("Peak","MGS")) %>%
                left_join(reg_df, by=c("site", "metric")) %>% 
                mutate(value_plot=case_when(p.value<0.05~value)))+
  geom_point(aes(x=daymet_precip, y=value))+
  geom_smooth(aes(x=daymet_precip, y=value_plot), method="lm")+
  ggpubr::stat_cor(
    aes(x = daymet_precip, y = value, group = interaction( metric, site),
        # col = interaction( metric, site),
        label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
    # p.accuracy = 0.05,
    label.x.npc = "left",
    label.y.npc = "top",
    show.legend=F,
    col="blue"
  )+
  facet_grid(vars(metric), vars(site), scales = "free")+
  theme_classic()+
  ylab("Value of phenological metric")+
  xlab("Total annual precipitation (mm)")
p_reg

dir.create(paste0(path, "output/"))
pdf(paste0(path, 'output/metric_regression.pdf'), width = 8, height = 8)
grid.arrange(annotate_figure(p_metric, fig.lab = "(a)"),
             annotate_figure(p_reg, fig.lab = "(b)"),
             layout_matrix=rbind(c(1),
                                 c(2),
                                 c(2))
)
dev.off()
