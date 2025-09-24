#### Synchron working group #####

# Travis McDevitt-Galles
# 09/24/2021
# title: 03_Synchrony_Pop_Extract

# Extracting predicted values of mosquito population trends

library(rstanarm)
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)

# loading in the seasonal abundance data and r function
load(here("Data", "Mosquito_Data_Clean.Rda"))
load(here("Data", "DailyPrismMod.Rda"))

source(here("R_Script", "02_Synchrony_GAM_Function.R"))

### Lets extract the sites we care about


names(contigus.df)

site.df <- unique(select(ungroup(complete.df), c("Site", "Plot")))

weather.df <- right_join(site.df, contigus.df, by = "Plot")

weather.df <- weather.df %>% filter(Site == "HARV" |
  Site == "ORNL" |
  Site == "SERC" |
  Site == "TALL" |
  Site == "UNDE" |
  Site == "WOOD" |
  Site == "WREF" |
  Site == "YELL")


weather.df <- weather.df %>%
  group_by(Site, Year, Lat, Long, DOY, Julian) %>%
  summarise(
    PPT = mean(PPT),
    TMAX = mean(TMAX),
    TMean7 = mean(Tmean7),
    PPT14 = mean(PPT14),
    Photoperiod = mean(Photoperiod)
  )


pred.df <- readRDS("Data/predict.rds")
pred1.df <- readRDS("Data/predict_year.rds")
pred3.df <- readRDS("Data/predict_long.rds")

full.df <- rbind.data.frame(pred.df, pred1.df, pred3.df)

full.df$SiteYear <- paste(full.df$Site, full.df$Year)

weather.df$SiteYear <- paste(weather.df$Site, weather.df$Year)

weather.df <- weather.df %>% filter(SiteYear %in% full.df$SiteYear)

length(unique(weather.df$SiteYear))
length(unique(full.df$SiteYear))


weather.df <- weather.df %>% filter(DOY %in% full.df$DOY)

weather.df <- weather.df %>% filter(Site == "HARV" |
  Site == "UNDE" |
  Site == "WOOD" |
  Site == "WREF" |
  Site == "YELL")

weather.df$Photo <- "Photoperiod"
photo <- ggplot(weather.df, aes(x = DOY, y = Photoperiod, color = Site, group = (SiteYear))) +
  geom_line(size = 2, alpha = .5) +
  facet_wrap(~Photo) +
  xlab("DOY") +
  ylab("Day length (hrs)") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.key.size = unit(1, "cm"),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = rel(1.25)),
    axis.line.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_text(size = rel(1.8)),
    axis.text.x = element_text(vjust = 0.5, color = "black"),
    axis.text.y = element_text(vjust = 0.5, color = "black"),
    axis.title.y = element_text(size = rel(1.8), angle = 90),
    strip.text.x = element_text(size = 20)
  )



weather.df$Precip <- "Precipitation"
precip <- ggplot(weather.df, aes(x = DOY, y = PPT14, color = Site, group = SiteYear)) +
  geom_line(size = 2, alpha = .5) +
  facet_wrap(~Precip) +
  xlab("DOY") +
  ylab("Total 14 day precip (mm)") +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.key.size = unit(1, "cm"),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = rel(1.25)),
    axis.line.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_text(size = rel(1.8)),
    axis.text.x = element_text(vjust = 0.5, color = "black"),
    axis.text.y = element_text(vjust = 0.5, color = "black"),
    axis.title.y = element_text(size = rel(1.8), angle = 90),
    strip.text.x = element_text(size = 20)
  )


weather.df$Temp <- "Temperature"
temp <- ggplot(weather.df, aes(x = DOY, y = TMean7, color = Site, group = SiteYear)) +
  geom_line(size = 2, alpha = .5) +
  facet_wrap(~Temp) +
  xlab("DOY") +
  ylab("7 day mean temperature (C)") +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.key.size = unit(1, "cm"),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = rel(1.25)),
    axis.line.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_text(size = rel(1.8)),
    axis.text.x = element_text(vjust = 0.5, color = "black"),
    axis.text.y = element_text(vjust = 0.5, color = "black"),
    axis.title.y = element_text(size = rel(1.8), angle = 90),
    strip.text.x = element_text(size = 20)
  )

library(patchwork)

(precip + temp + photo)

at1 <- left_join(weather.df, full.df, by = c("Site", "Year", "DOY", "SiteYear"))

at1 <- at1[complete.cases(at1), ]

plot.df <- at1 %>% filter(SciName == "Aedes vexans" &
  Site == "HARV" &
  Year == "2016")

ggplot(plot.df, aes(x = DOY, y = Count)) +
  geom_line() +
  geom_line(aes(x = DOY, y = TMean7 / 5), color = "blue")


### Building a simple data set to store the values of population overlap
library(sfsmisc)

t <- 1


for (s in 1:length(unique(at1$Site))) {
  focSite <- unique(at1$Site)[s]
  site.df <- at1 %>% filter(Site == focSite)

  for (y in 1:length(unique(site.df$Year))) {
    focYear <- unique(site.df$Year)[y]

    pred.df <- site.df %>% filter(Year == focYear)

    weather.combo <- as.data.frame(expand.grid(
      unique(pred.df$SciName),
      "Temp"
    ))

    temp.df <- weather.combo

    colnames(temp.df) <- c("Sp1", "Temp")

    temp.df$Site <- focSite
    temp.df$Year <- focYear


    temp.df$Corr <- NA

    for (i in 1:nrow(temp.df)) { # Number of sites (# of combinations?)

      # Isolating species 1 and species 2
      Sp1 <- temp.df$Sp1[i]


      trial1.df <- filter(ungroup(pred.df), SciName == Sp1)

      trial.df <- trial1.df %>% select(c("DOY", "Count", "TMean7"))


      colnames(trial.df)[2] <- "Sp1"

      trial.df[is.na(trial.df)] <- 0


      temp.df$Corr[i] <- cor(trial.df$Sp1, trial.df$TMean7, method = "spearman")
    }
    if (t == 1) {
      overlap.df <- temp.df

      t <- t + 1
    } else {
      overlap.df <- rbind.data.frame(overlap.df, temp.df)
      t <- t + 1
    }
  }
}


hold1.df <- overlap.df
hold2.df <- overlap.df
hold3.df <- overlap.df

colnames(hold1.df)[1:2] <- c("SciName", "Abio")
colnames(hold2.df)[1:2] <- c("SciName", "Abio")
colnames(hold3.df)[1:2] <- c("SciName", "Abio")

overlap.df <- rbind.data.frame(hold1.df, hold2.df, hold3.df)

overlap.df <- overlap.df %>% filter(Site == "HARV" |
  Site == "UNDE" |
  Site == "WOOD" |
  Site == "WREF" |
  Site == "YELL")
ggplot(overlap.df, aes(x = Year, y = Corr, fill = Site)) +
  theme_classic() +
  geom_hline(yintercept = 0, alpha = .5) +
  geom_boxplot() +
  xlab("Abiotic variables") +
  ylab("Spearman correlation") +
  facet_wrap(~Abio, scales = "free_x") +
  theme(
    legend.position = "top",
    legend.key.size = unit(1, "cm"),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = rel(1.25)),
    axis.line.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_text(size = rel(1.8)),
    axis.text.x = element_text(vjust = 0.5, color = "black"),
    axis.text.y = element_text(vjust = 0.5, color = "black"),
    axis.title.y = element_text(size = rel(1.8), angle = 90),
    strip.text.x = element_text(size = 20)
  )




overlap.df %>%
  filter(Site == "UNDE") %>%
  filter(SciName == "Aedes communis" |
    SciName == "Aedes punctor" |
    SciName == "Coquillettidia perturbans") %>%
  ggplot(aes(x = Abio, y = Corr, color = Year)) +
  theme_classic() +
  geom_hline(yintercept = 0, alpha = .5) +
  geom_boxplot(color = "black") +
  geom_point(size = 3) +
  xlab("Abiotic variables") +
  ylab("Spearman correlation") +
  facet_wrap(~SciName, scales = "free_x") +
  theme(
    legend.position = "top",
    legend.key.size = unit(1, "cm"),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = rel(1.25)),
    axis.line.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_text(size = rel(1.8)),
    axis.text.x = element_text(vjust = 0.5, color = "black"),
    axis.text.y = element_text(vjust = 0.5, color = "black"),
    axis.title.y = element_text(size = rel(1.8), angle = 90),
    strip.text.x = element_text(size = 14)
  )




overlap.df %>%
  ggplot(aes(x = Corr, fill = Abio)) +
  geom_density(alpha = .75) +
  facet_wrap(~Site) +
  xlab("Abiotic overlap") +
  ylab("Density") +
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



overlap.df$Site <- factor(overlap.df$Site, levels = c("HARV", "SERC", "ORNL", "TALL"))

overlap.df %>%
  filter(Overlap < 1) %>%
  ggplot(aes(x = reorder(Sp1, -Overlap), y = Overlap)) +
  geom_boxplot() +
  geom_jitter(width = .1) +
  xlab("Species identity") +
  ylab("Overlap") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_text(size = rel(1.8)),
    axis.text.x = element_text(vjust = 0.5, size = 0, color = "black"),
    axis.title.y = element_text(size = rel(1.8), angle = 90),
    axis.text.y = element_text(vjust = 0.5, size = 16),
    strip.text.x = element_text(size = 20)
  )



overlap.df$Comb <- paste(overlap.df$Sp1, overlap.df$Sp2)




overlap.df %>%
  filter(Overlap < 1) %>%
  ggplot(aes(y = reorder(Comb, -Overlap), x = Overlap, color = Year)) +
  geom_point(alpha = .25) +
  xlab("Species identity") +
  ylab("Overlap") +
  theme_classic() +
  facet_wrap(~Site, scales = "free_y") +
  theme(
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_text(size = rel(1.8)),
    axis.text.x = element_text(vjust = 0.5, size = 16, color = "black"),
    axis.title.y = element_text(size = rel(1.8), angle = 90),
    axis.text.y = element_text(vjust = 0.5, size = 0),
    strip.text.x = element_text(size = 20)
  )
