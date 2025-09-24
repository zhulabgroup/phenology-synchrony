#### Synchron working group #####

# Travis McDevitt-Galles
# 09/24/2021
# title: 03_Synchrony_Pop_Extract

# Extracting predicted values of mosquito population trends

library(rstanarm)
library(tidyverse)

# loading in the seasonal abundance data and r function

source(here("R_Script", "02_Synchrony_GAM_Function.R"))

### Lets extract the sites we care about
path <- "./mosquito/"
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


count.df <- full.df %>%
  group_by(Site, Year, SciName, DOY) %>%
  summarise(
    Count = mean(Count) # ,
    # Q_2.5 = mean(Q_2.5),
    # Q_97.5 = mean(Q_97.5),
  ) %>%
  ungroup() %>%
  group_by(Site, Year, SciName) %>%
  mutate(cumCount = cumsum(Count)) %>%
  mutate(proCount = cumCount / max(cumCount))

count.df %>%
  ggplot(aes(x = DOY, y = proCount, color = interaction(SciName, Year))) +
  geom_line(size = 2, alpha = .5) +
  facet_wrap(. ~ Site) +
  theme(legend.position = "none")

## generating the number of days between the 10% quantile nad 90th quanitile
## Higher values indicate

range.df <- count.df %>%
  filter(proCount >= .1 & proCount <= .9) %>%
  group_by(SciName, Year, Site) %>%
  summarise(TDays = n())

ggplot(range.df, aes(x = SciName, y = TDays, col = Year)) +
  geom_point() +
  facet_wrap(. ~ Site)

ggplot(
  range.df %>%
    filter(Site != "ORNL" & Site != "SERC" & Site != "TALL" & Site != "WREF" & Site != "YELL") %>%
    left_join(range.df %>% group_by(SciName, Site) %>% distinct(Year) %>% summarise(n = n()),
      by = c("Site", "SciName")
    ) %>%
    filter(n > 3) %>%
    mutate(genus = str_split(SciName, " ", simplify = T)[1]),
  aes(x = Year, y = TDays, col = SciName, fill = SciName, group = SciName)
) +
  geom_line() +
  # geom_smooth(method="lm", se=F, alpha=0.3)+
  facet_wrap(. ~ Site, ncol = 1)

range.df.proc <- range.df %>%
  left_join(range.df %>% group_by(SciName) %>% distinct(Year, Site) %>% summarise(n = n()),
    by = c("SciName")
  ) %>%
  filter(n >= 5) %>%
  mutate(Year = as.numeric(Year)) %>%
  left_join(daymet_df, by = c("Site", "Year"))
ggplot(range.df.proc, aes(x = mat, y = TDays, col = SciName)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, alpha = 0.3) +
  facet_wrap(. ~ SciName) +
  theme_classic() +
  guides(col = F)

reg.df <- range.df.proc %>%
  group_by(SciName) %>%
  do(broom::tidy(lm(TDays ~ mat, .))) %>%
  filter(term == "mat") %>%
  dplyr::select(-term, -statistic)
ggplot(reg.df) +
  geom_point(aes(x = SciName, y = estimate)) +
  geom_errorbar(aes(x = SciName, ymin = estimate - 1.95 * std.error, ymax = estimate + 1.95 * std.error)) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic()
summary(reg.df$estimate > 0)

## getting phenometrics

pheno.df <- unique(select(new.df, c("Site", "Year", "SciName", "SiteYearSpp")))

pheno.df$First <- NA
pheno.df$Last <- NA
pheno.df$Duration <- NA

for (i in 1:length(unique(new.df$SiteYearSpp))) {
  foc_spp <- unique(new.df$SiteYearSpp)[i]

  foc.df <- new.df %>% filter(SiteYearSpp == foc_spp)

  pheno.df$First[pheno.df$SiteYearSpp == foc_spp] <- foc.df$DOY[which.min(foc.df$ProCount[foc.df$ProCount > .0499])]

  pheno.df$Last[pheno.df$SiteYearSpp == foc_spp] <- foc.df$DOY[which.max(foc.df$ProCount[foc.df$ProCount < .951])]
}

pheno.df$Duration <- pheno.df$Last - pheno.df$First

hist(pheno.df$Duration)
hist(pheno.df$First)


pheno.df %>%
  ggplot(aes(x = Duration, fill = Site)) +
  geom_density(alpha = .25) +
  facet_wrap(~SciName, scales = "free") +
  theme(legend.position = "None")

pheno.df %>%
  ggplot(aes(x = Site, y = Duration, fill = Year)) +
  geom_boxplot(alpha = .25) +
  # geom_jitter(aes(x=Site, y =Duration, color= SciName ),
  # width = .1, size=2,alpha=.4) +
  theme(legend.position = "None")



pheno.df %>%
  filter(Site != "ORNL" & Site != "SERC" & Site != "TALL") %>%
  ggplot(aes(x = Site, y = SynchPro, fill = Year)) +
  geom_boxplot(alpha = .25) +
  xlab("Site") +
  theme_classic() +
  ylab("Synchrony (1/sd)")




full.df <- full.df %>%
  group_by(Site, Year, SciName, SiteYearSpp) %>%
  mutate(MaxCount = sum(Count))


full.df$ProCount <- full.df$Count / full.df$MaxCount

synch.df <- full.df %>%
  group_by(Site, Year, SciName, SiteYearSpp) %>%
  summarize(
    synchrony = 1 / sd(Count),
    SynchPro = 1 / sd(ProCount)
  ) %>%
  ungroup()




pheno.df$Synchrony <- synch.df$synchrony

pheno.df$SynchPro <- synch.df$SynchPro

pairs(pheno.df[, 5:9])

hist(pheno.df$Synchrony)

pheno.df[which.min(pheno.df$SynchPro), ]


plot(x = pheno.df$SynchPro, y = (pheno.df$Synchrony))


library(lme4)

test.df <- pheno.df %>% filter(Site != "ORNL" & Site != "SERC" & Site != "TALL")

at1 <- lmer((Synchrony) ~ (1 | SciName) + (1 | Site) + (1 | Year), data = test.df)
summary(at1)

merge.df <- select(pheno.df, c(
  "SiteYearSpp", "First", "Last", "Duration",
  "Synchrony", "SynchPro"
))

full.df <- left_join(full.df, merge.df, by = "SiteYearSpp")


hist(full.df$SynchPro)

full.df %>%
  filter(log10(SynchPro) < 2.5) %>%
  ggplot(aes(x = DOY, y = Count, color = Site)) +
  geom_line() +
  facet_wrap(~SiteYearSpp,
    scales = "free_y"
  )

pheno.df %>%
  filter(Site != "ORNL" & Site != "SERC" & Site != "TALL") %>%
  ggplot(aes(x = Site, y = SynchPro, fill = Site)) +
  geom_boxplot(alpha = .25) +
  xlab("Site") +
  theme_classic() +
  ylab("Synchrony (1/sd)") +
  facet_wrap(~Year, scales = "free_x") +
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
