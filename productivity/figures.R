NE_df<-data.frame(NE2_table) %>% 
  `colnames<-` (c("metric", "US-Ne2.IR.Mean", "US-Ne2.IR.Stdev", "US-Ne3.RF.Mean", "US-Ne3.RF.Stdev")) %>% 
  tidyr::gather(key="sitestat", value="value", -metric) %>% 
  rowwise () %>% 
  mutate(site=str_split(sitestat,pattern="\\.",simplify = T)[2]) %>% 
  mutate(stat=str_split(sitestat,pattern="\\.",simplify = T)[3]) %>% 
  ungroup() %>% 
  dplyr::select(-sitestat) %>% 
  spread(key="stat", value="value") %>% 
  mutate(Mean=as.numeric(Mean),
         Stdev=as.numeric(Stdev)) %>% 
  filter(!metric %in% c("MSP", "MAU")) %>% 
  mutate(metric=as.factor(metric)) %>% 
  mutate(metric=fct_relevel(metric, levels=c("SOS", "EOS", "LOS", "MGS", "Peak")))
  

p_metric<-ggplot(NE_df)+
  geom_point(aes(x=site, y=Mean))+
  geom_errorbar(aes(x=site, ymin=Mean-1.95*Stdev, ymax=Mean+1.95*Stdev))+
  facet_wrap(.~metric, nrow=1, scales = "free_y")+
  theme_classic()+
  ylab("Value of phenological metric")+
  xlab("Site")


NE_ann_df<-bind_rows(NE2_ann %>% mutate(site="IR"), NE3_ann%>% mutate(site="RF"))

p_reg<-ggplot(NE_ann_df %>% 
         dplyr::select(site, year, daymet_precip, Peak=peak, MGS=mgs) %>% 
         tidyr::gather(key="metric", value="value", -site, -year,-daymet_precip) %>% 
           mutate(value_plot=case_when(site=="RF"~value)))+
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

library(gridExtra)

path<-"./productivity/"
dir.create(paste0(path, "output/"))
pdf(paste0(path, 'output/metric_regression.pdf'), width = 8, height = 8)
grid.arrange(annotate_figure(p_metric, fig.lab = "(a)"),
             annotate_figure(p_reg, fig.lab = "(b)"),
             layout_matrix=rbind(c(1),
                                 c(2),
                                 c(2))
)
dev.off()
