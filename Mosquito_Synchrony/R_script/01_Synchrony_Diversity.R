###### Synchrony Working Group -- Community synchrony  #######

# Travis McDevitt-Galles
# 09/24/2021
# title: 01_Synchrony_Diversity

library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan)
library(here)

# loading clean mosquito data

load( here("Data", "Mosquito_Data_Clean.Rda" ) )

names(complete.df)

# lots of data : lets first break this down to look at patterns of richness

## reassigning full dataset into new dataset
rich.df <- complete.df

# Column for whether a species is present
rich.df$Pres <- NA 
# Count greater than 1 == present
rich.df$Pres[ rich.df$Count > 0 ] <- 1
# Count less than 1 == absent
rich.df$Pres[ rich.df$Count == 0 ] <- 0

## Lets sum across Sample events Site, year 

rich.df <- rich.df %>% 
           group_by( Domain, Site, Plot,
                     Lat, Long, Year, 
                     Elev, TrapEvent ) %>% 
           summarize(
            DOY = min(DOY),
            Rich = sum(Pres)
          ) %>%  ungroup()

rich.df <- rich.df %>% filter(Site== "HARV"|
                     Site == "SERC" |
                     Site == "ORNL" |
                     Site == "TALL")

rich.df$Site <- factor(rich.df$Site ,levels= c("HARV", "SERC", "ORNL", "TALL") )
rich.df %>% filter(Year == "2018")  %>% 
  ggplot(aes(x=DOY, y=Rich, color= Site)) + geom_point(alpha=1) +
  ylab("Alpha richness") + facet_wrap(~Site, scales="free_y")+
  xlab("Day of year") + ylab("Alpha diversity") + 
  theme_classic()+
  theme( legend.position = "none",
         axis.line.x = element_line(color="black") ,
         axis.ticks.y = element_line(color="black"),
         axis.ticks.x = element_line(color="black"),
         axis.title.x = element_text(size = rel(1.8)),
         axis.text.x  = element_text(vjust=0.5, size=16, color = "black"),
         axis.title.y = element_text(size = rel(1.8), angle = 90) ,
         axis.text.y=element_text(vjust=0.5, size=16),
         strip.text.x = element_text(size=20)
  )


 







library(rstanarm)
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)


### Lets extract the sites we care about


pred.df <- readRDS("Data/predict.rds")
pred1.df <- readRDS("Data/predict_year.rds")
pred3.df <- readRDS("Data/predict_long.rds")

pred.df$SiteYearSpp <- paste(pred.df$Site,pred.df$Year,pred.df$SciName)
pred1.df$SiteYearSpp <- paste(pred1.df$Site,pred1.df$Year,pred1.df$SciName)
pred3.df$SiteYearSpp <- paste(pred3.df$Site,pred3.df$Year,pred3.df$SciName)


unique(pred.df$SiteYearSpp)
unique(pred1.df$SiteYearSpp)
unique(pred3.df$SiteYearSpp)

pred.df <- filter(pred.df, Site !="HARV")

pred3.df <- filter(pred3.df, Site != "HARV")


full.df <- rbind.data.frame(pred.df, pred1.df,pred3.df)

full.df$SiteYear <- paste(full.df$Site,full.df$Year)

full.df <- full.df %>%  group_by(Site, Year, SciName, DOY,SiteYear, SiteYearSpp) %>% 
  summarize(
    Count = mean(Count)
  ) %>% ungroup()

full.df <- full.df %>% select(c("Site", "Year", "DOY", "SiteYear","SciName","Count"))

t <- 1
for(i in 1:length(full.df$SiteYear )){
  at1 <- full.df %>% filter(SiteYear == unique(full.df$SiteYear)[i] )
  test.df <- at1 %>% pivot_wider( names_from = SciName, values_from=Count)
  
  nCols <- dim(test.df)[2]
  test.df$shannon <- NA
  
  

  test.df$shannon <- diversity(as.matrix(test.df[,5:nCols]), index = "shannon")

  plot.df <- test.df %>% select(c( "Site", "Year", "DOY", 
                                   "shannon"))
plot.df <- unique(plot.df)

if( t == 1){
    shannon.df <- plot.df
  t <- t +1 
  }else{
    shannon.df <- rbind.data.frame(shannon.df,plot.df)
  }
}

save <-shannon.df
shannon.df$Site <- factor(shannon.df$Site,levels= c("HARV", "SERC", "ORNL", "TALL"))

shannon.df %>% filter( Site !="ORNL" & Site != "SERC" & Site != "TALL" & 
                         Site != "WREF" & Site != "YELL") %>% 
ggplot(aes(x=DOY, y=shannon, color=Year))+geom_line(size=2, alpha=.5)+
  facet_wrap(~Site)+ xlab("Day of year") + ylab("Shannon diversity") + 
  theme_classic()+
  theme(
    legend.key.size = unit(1, "cm"),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    axis.line.x = element_line( color = "black" ) ,
    axis.ticks.y = element_line( color = "black" ),
    axis.ticks.x = element_line( color = "black" ),
    axis.title.x = element_text( size = rel( 1.8 ) ),
    axis.text.x = element_text( vjust = 0.5, color = "black" ),
    axis.text.y = element_text( vjust = 0.5, color = "black" ),
    axis.title.y = element_text( size = rel( 1.8 ), angle = 90 ),
    strip.text.x = element_text( size = 20 ) )

