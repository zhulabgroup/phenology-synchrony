###### Synchrony Working Group -- Community synchrony  #######

# Travis McDevitt-Galles
# 09/24/2021
# title: 01_Synchrony_Diversity

library(tidyverse)
library(vegan)
path <- "./mosquito/"
# loading clean mosquito data

load(paste0(path, "/data/Mosquito_Data_Clean.Rda"))

names(complete.df)

# lots of data : lets first break this down to look at patterns of richness

## reassigning full dataset into new dataset
rich.df <- complete.df

# Column for whether a species is present
rich.df$Pres <- NA
# Count greater than 1 == present
rich.df$Pres[rich.df$Count > 0] <- 1
# Count less than 1 == absent
rich.df$Pres[rich.df$Count == 0] <- 0

## Lets sum across Sample events Site, year

rich.df <- rich.df %>%
  group_by(
    Domain, Site, Plot,
    Lat, Long, Year,
    Elev, TrapEvent
  ) %>%
  summarize(
    DOY = min(DOY),
    Rich = sum(Pres)
  ) %>%
  ungroup()

rich.df <- rich.df %>% filter(Site == "HARV" |
  Site == "SERC" |
  Site == "ORNL" |
  Site == "TALL")

rich.df$Site <- factor(rich.df$Site, levels = c("HARV", "SERC", "ORNL", "TALL"))
rich.df %>%
  filter(Year == "2018") %>%
  ggplot(aes(x = DOY, y = Rich, color = Site)) +
  geom_point(alpha = 1) +
  ylab("Alpha richness") +
  facet_wrap(~Site, scales = "free_y") +
  xlab("Day of year") +
  ylab("Alpha diversity") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_text(size = rel(1.8)),
    axis.text.x = element_text(vjust = 0.5, size = 16, color = "black"),
    axis.title.y = element_text(size = rel(1.8), angle = 90),
    axis.text.y = element_text(vjust = 0.5, size = 16),
    strip.text.x = element_text(size = 20)
  )










library(rstanarm)
library(tidyverse)


### Lets extract the sites we care about


pred.df <- readRDS(paste0(path, "data/predict.rds"))
pred1.df <- readRDS(paste0(path, "data/predict_year.rds"))
pred3.df <- readRDS(paste0(path, "data/predict_long.rds"))

pred.df$SiteYearSpp <- paste(pred.df$Site, pred.df$Year, pred.df$SciName)
pred1.df$SiteYearSpp <- paste(pred1.df$Site, pred1.df$Year, pred1.df$SciName)
pred3.df$SiteYearSpp <- paste(pred3.df$Site, pred3.df$Year, pred3.df$SciName)


unique(pred.df$SiteYearSpp)
unique(pred1.df$SiteYearSpp)
unique(pred3.df$SiteYearSpp)

pred.df <- filter(pred.df, Site != "HARV")

pred3.df <- filter(pred3.df, Site != "HARV")


full.df <- rbind.data.frame(pred.df, pred1.df, pred3.df)

# full.df$SiteYear <- paste(full.df$Site, full.df$Year)

shannon.df <- full.df %>%
  group_by(Site, Year, SciName, DOY) %>%
  summarise(
    Count = mean(Count)
  ) %>%
  ungroup() %>%
  group_by(Site, Year, DOY) %>%
  summarise(shannon = diversity(as.matrix(Count), index = "shannon")) %>%
  mutate(Site = factor(Site, levels = c("HARV", "SERC", "ORNL", "TALL")))

shannon.df %>%
  filter(Site != "ORNL" & Site != "SERC" & Site != "TALL" & Site != "WREF" & Site != "YELL") %>%
  ggplot(aes(x = DOY, y = shannon, color = Year)) +
  geom_line(size = 2, alpha = .5) +
  facet_wrap(~Site) +
  xlab("Day of year") +
  ylab("Shannon diversity") +
  theme_classic() +
  theme(
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    axis.line.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_text(size = rel(1.8)),
    axis.text.x = element_text(vjust = 0.5, color = "black"),
    axis.text.y = element_text(vjust = 0.5, color = "black"),
    axis.title.y = element_text(size = rel(1.8), angle = 90),
    strip.text.x = element_text(size = 20)
  ) +
  facet_wrap(. ~ Site)
