#### Synchron working group #####

# Travis McDevitt-Galles
# 09/24/2021
# title: 03_Synchrony_Pop_Extract

# Extracting predicted values of mosquito population trends

library(rstanarm)
library(ggplot2)
library(tidyr)
library(dplyr)

# loading in the seasonal abundance data and r function
load(paste0(path, "data/Mosquito_Data_Clean.Rda"))

source(here("R_Script", "02_Synchrony_GAM_Function.R"))

### What should be the abundance cut off for seasonal mosquito count?

complete.df <- complete.df %>%
  group_by(SciName, Year, Site) %>%
  mutate(Total_Count = sum(Count))

hist(log10(complete.df$Total_Count))

mos.df <- complete.df %>% filter(SciName != "Wyeomyia sp." &
  SciName != "Psorophora sp." &
  SciName != "Uranotaenia sp." &
  SciName != "Mansonia sp." &
  SciName != "Culiseta sp." &
  SciName != "Culicidae sp." &
  SciName != "Culex sp." &
  SciName != "Anopheles sp." &
  SciName != "Aedes sp.")

at1 <- mos.df %>% filter(Site == "YELL")



unique(at1$SciName)

hist(log10(at1$Total_Count))

at1 %>%
  filter(Year == "2018") %>%
  filter(Total_Count > 10) %>%
  ggplot(aes(x = DOY, y = Count, color = SciName)) +
  geom_point(size = 2, alpha = .5) +
  facet_wrap(~Year, scales = "free_y")



at1 %>%
  filter(Year == "2014") %>%
  filter(Total_Count > 10) %>%
  ggplot(aes(x = DOY, y = log10(Count + 1), color = SciName)) +
  geom_point(size = 2, alpha = .5) +
  facet_wrap(~SciName, scales = "free_y") +
  theme(legend.position = "none")

test.df <- at1 %>%
  filter(Year == "2019") %>%
  filter(Total_Count > 10)


gam.df <- filter(gam.df, TotalWeight > 0)


attempt9 <- pheno_gam(test.df, Site = "UNDE")



attempt9[[2]]

######### year
play.df <- attempt9[[1]]



# long.df <- play.df





long.df <- rbind.data.frame(long.df, play.df)

long.df %>%
  ggplot(aes(x = DOY, y = log10(Count + 1), color = SciName)) +
  geom_line(size = 1, alpha = .5) +
  facet_wrap(~Site, scales = "free_y") +
  theme(legend.position = "none")



play4.df <- attempt4[[1]]
play4.df %>%
  ggplot(aes(x = DOY, y = log10(Count + 1), color = SciName, group = SciName)) +
  geom_line(size = 1, alpha = .5)


play.all.df <- rbind.data.frame(play.all.df, play4.df)


# saveRDS(play.all.df ,"/Data/predict.rds")

# saveRDS(year.df ,"Data/predict_year.rds")
pred.df <- readRDS("Data/predict.rds")
pred1.df <- readRDS("Data/predict_year.rds")
pred3.df <- readRDS("Data/predict_long.rds")

pred.df$SiteYearSpp <- paste(pred.df$Site, pred.df$Year, pred.df$SciName)
pred1.df$SiteYearSpp <- paste(pred1.df$Site, pred1.df$Year, pred1.df$SciName)
pred3.df$SiteYearSpp <- paste(pred3.df$Site, pred3.df$Year, pred3.df$SciName)


unique(pred.df$SiteYearSpp)
unique(pred1.df$SiteYearSpp)
unique(pred3.df$SiteYearSpp)

pred.df <- filter(pred.df, Site != "HARV")

pred3.df <- filter(pred3.df, Site != "HARV")


full.df <- rbind.data.frame(pred.df, pred1.df, pred3.df)


full.df <- full.df %>%
  group_by(Site, Year, SciName, DOY, SiteYearSpp) %>%
  summarize(
    Count = mean(Count),
    Q_2.5 = mean(Q_2.5),
    Q_97.5 = mean(Q_97.5),
  ) %>%
  ungroup()


full.df$SiteYearSpp <- paste(full.df$Year, full.df$Site, full.df$SciName)
full.df %>%
  filter(Year == "2018" & DOY < 300) %>%
  filter(Site != "ORNL" & Site != "SERC" & Site != "TALL") %>%
  ggplot(aes(x = DOY, y = log10(Count + 1), color = SciName)) +
  geom_line(size = 2, alpha = .5) +
  facet_wrap(~Site, scales = "free_y") +
  theme(legend.position = "none") +
  scale_color_ordinal(name = "Species name") +
  xlab("DOY") +
  ylab("Estimated log10 mosquito population") +
  theme_classic() +
  theme(
    legend.position = c(.85, .2),
    axis.line.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_text(size = rel(1.8)),
    axis.text.x = element_text(vjust = 0.5, size = 16, color = "black"),
    axis.title.y = element_text(size = rel(1.8), angle = 90),
    axis.text.y = element_text(vjust = 0.5, size = 16),
    strip.text.x = element_text(size = 20)
  )





at1 <- pred3.df
### Building a simple data set to store the values of population overlap
library(sfsmisc)

t <- 1
p <- 1
over_plot <- list() # empty list to store gam plots



for (s in 1:length(unique(at1$Site))) {
  focSite <- unique(at1$Site)[s]
  site.df <- at1 %>% filter(Site == focSite)

  for (y in 1:length(unique(site.df$Year))) {
    focYear <- unique(site.df$Year)[y]

    pred.df <- site.df %>% filter(Year == focYear)

    spp.combo <- as.data.frame(expand.grid(
      unique(pred.df$SciName),
      unique(pred.df$SciName)
    ))

    temp.df <- spp.combo

    colnames(temp.df) <- c("Sp1", "Sp2")

    temp.df$Site <- focSite
    temp.df$Year <- focYear


    temp.df$Overlap <- NA
    temp.df$Corr <- NA

    for (i in 1:nrow(temp.df)) { # Number of sites (# of combinations?)

      # Isolating species 1 and species 2
      Sp1 <- temp.df$Sp1[i]
      Sp2 <- temp.df$Sp2[i]


      trial1.df <- filter(pred.df, SciName == Sp1)

      trial.df <- trial1.df %>% select(c("DOY", "Count"))

      colnames(trial.df)[2] <- "Sp1"

      trial2.df <- filter(pred.df, SciName == Sp2)

      trial2.df <- trial2.df %>% select(c("DOY", "Count"))

      colnames(trial2.df)[2] <- "Sp2"

      trial.df <- full_join(trial.df, trial2.df, by = "DOY")

      trial.df[is.na(trial.df)] <- 0
      ## How you get the interactions

      trial.df$inter <- pmin(trial.df$Sp1, trial.df$Sp2)

      ## Getting the total
      if (sum(trial.df$Sp1) > sum(trial.df$Sp2)) {
        total <- integrate.xy(trial.df$DOY, trial.df$Sp1)
        # total <- integrate.xy(trial.df$DOY, trial.df$Sp1) +
        # integrate.xy(trial.df$DOY, trial.df$Sp2)
      } else {
        total <- integrate.xy(trial.df$DOY, trial.df$Sp2)
      }

      intersection <- integrate.xy(trial.df$DOY, trial.df$inter)

      temp.df$Overlap[i] <- intersection / total

      temp.df$Corr[i] <- cor(trial.df$Sp1, trial.df$Sp2)
      plot1 <- ggplot() +
        geom_line(data = trial.df, aes(x = DOY, y = Sp1)) +
        geom_ribbon(
          data = trial.df, aes(x = DOY, ymax = Sp1, ymin = 0),
          alpha = .5, fill = "light blue"
        ) +
        geom_line(data = trial.df, aes(x = DOY, y = Sp2)) +
        geom_ribbon(
          data = trial.df, aes(x = DOY, ymax = Sp2, ymin = 0),
          alpha = .5, fill = "red"
        ) +
        xlab("Day of year") +
        theme_classic() +
        ylab(" Estiamted mosquito count") +
        theme(
          axis.line.x = element_line(color = "black"),
          axis.ticks.y = element_line(color = "black"),
          axis.ticks.x = element_line(color = "black"),
          axis.title.x = element_text(size = rel(1.8)),
          axis.text.x = element_text(vjust = 0.5, color = "black"),
          axis.text.y = element_text(vjust = 0.5, color = "black"),
          axis.title.y = element_text(size = rel(1.8), angle = 90),
          strip.text.x = element_text(size = 20)
        )

      over_plot[[p]] <- plot1
      p <- p + 1
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

plot(x = overlap.df$Overlap, y = overlap.df$Corr)
overlap.df %>%
  filter(Overlap < 1) %>%
  ggplot(aes(x = Year, y = Overlap, fill = Site)) +
  geom_boxplot() +
  xlab("Year") +
  ylab("Population overlap") +
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



overlap.df %>%
  filter(Overlap < 1) %>%
  ggplot(aes(x = Overlap, fill = Year)) +
  geom_density(alpha = .75) +
  facet_wrap(~Site) +
  xlab("Population overlap") +
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
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .1, alpha = .25) +
  xlab("Species identity") +
  ylab("Overlap") +
  facet_wrap(~Site, scales = "free_x") +
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
