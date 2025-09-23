library(tidyverse)
library(lubridate)
path <- "./case study A flux/"

# read in flux data
list.files(str_c(path, "data"), ".zip", full.names = T) %>%
  purrr::walk(~ unzip(.x, exdir = str_c(path, "data/", str_remove(.x, ".zip"))))

df_flux <- bind_rows(
  list.files(str_c(path, "data/AMF_US-Ne2_BASE-BADM_18-5/"), ".csv", full.names = T) %>%
    read_csv(skip = 2, na = c("-9999")) %>%
    mutate(site = "IR"),
  list.files(str_c(path, "data/AMF_US-Ne3_BASE-BADM_19-5/"), ".csv", full.names = T) %>%
    read_csv(skip = 2, na = c("-9999")) %>%
    mutate(site = "RF")
) %>%
  mutate(time = lubridate::ymd_hm(TIMESTAMP_START)) %>%
  arrange(time) %>%
  select(site, time, gpp = GPP_PI_F_1_1_1, reco = RECO_PI_F_1_1_1, le = LE_PI_F_1_1_1, h = H_PI_F_1_1_1) %>%
  mutate(date = lubridate::date(time)) %>%
  group_by(site, date) %>%
  summarise(across(c(gpp, reco, le, h), ~ mean(.x, na.rm = TRUE))) %>%
  drop_na() %>%
  ungroup() %>%
  mutate(
    year = lubridate::year(date),
    month = lubridate::month(date)
  ) %>%
  pivot_longer(cols = c(gpp, reco, le, h), names_to = "variable", values_to = "value") %>%
  # rolling means
  group_by(site, variable) %>%
  arrange(date) %>%
  mutate(value = zoo::rollmean(value, 14, fill = NA, align = "right")) %>%
  ungroup()

df_daymet <- bind_rows(
  read_csv(paste0(path, "data/NE_211742_lat_41.1649_lon_-96.4701_2022-02-02_122643.csv"), skip = 7) %>%
    mutate(site = "IR"),
  read_csv(paste0(path, "data/NE_311742_lat_41.1797_lon_-96.4397_2022-02-02_122805.csv"), skip = 7) %>%
    mutate(site = "RF")
) %>%
  mutate(date = as.Date(paste(year, yday, sep = "-"), format = "%Y-%j")) %>%
  rename(prcp = `prcp (mm/day)`, tmax = `tmax (deg c)`, tmin = `tmin (deg c)`) %>%
  mutate(temp = (tmax + tmin) / 2) %>%
  select(site, date, prcp, temp) %>%
  mutate(
    year = lubridate::year(date),
    month = lubridate::month(date)
  ) %>%
  pivot_longer(cols = c(prcp, temp), names_to = "variable", values_to = "value") %>%
  # rolling means
  group_by(site) %>%
  arrange(date) %>%
  mutate(
    value = zoo::rollmean(value, 14, fill = NA, align = "right")
  ) %>%
  ungroup()

df_combined <- bind_rows(df_flux, df_daymet) %>%
  filter(year >= 2002, year <= 2015) %>%
  mutate(variable_unit = case_when(
    variable == "gpp" ~ "Gross~primary~productivity~(\u00B5*mol~CO[2]~m^-2~s^-1)",
    variable == "le" ~ "Latent~heat~flux~(W~m^-2)",
    variable == "prcp" ~ "Precipitation~(mm)",
    variable == "temp" ~ "Temperature~(degree*C)"
  )) %>%
  mutate(variable = factor(variable, levels = c("temp", "prcp", "le", "gpp"))) %>%
  arrange(variable) %>%
  mutate(variable_unit = factor(variable_unit, levels = unique(variable_unit)))

# first time crossing thresholds

# gpp threshold: 5
# le threshold: 50
# temp threshold: 15
# prcp threshold: 10

df_SOS <- df_combined %>%
  group_by(site, year, variable, variable_unit) %>%
  filter((variable == "gpp" & value >= 0.5 * max(value, na.rm = T)) |
    (variable == "le" & value >= 0.5 * max(value, na.rm = T)) |
    (variable == "temp" & value >= 15) |
    (variable == "prcp" & value >= 10)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(doy = lubridate::yday(date)) %>%
  select(site, year, variable, variable_unit, doy, date)

p_ts <- df_combined %>%
  filter(variable %in% c("gpp", "le", "temp")) %>%
  filter(year %in% 2012:2013) %>%
  ggplot() +
  geom_line(aes(x = date, y = value, col = site)) +
  geom_vline(
    data = df_SOS %>% filter(variable %in% c("gpp", "le", "temp")),
    aes(xintercept = date, col = site),
    show.legend = F
  ) +
  scale_color_viridis_d(option = "D", begin = 0.25, end = 0.75) +
  facet_wrap(. ~ variable_unit, scales = "free", ncol = 1, labeller = label_parsed) +
  theme_classic() +
  labs(x = "Time", y = "", col = "Site")

p_corr1 <- df_SOS %>%
  select(-date, -variable_unit) %>%
  pivot_wider(names_from = variable, values_from = doy) %>%
  ggplot(aes(x = temp, y = le, col = site)) +
  geom_point() +
  geom_smooth(method = MASS::rlm, se = T) +
  ggpubr::stat_cor(aes(x = temp, y = le), label.x.npc = 0.1, label.y.npc = 0.15, show.legend = F, inherit.aes = F) +
  scale_color_viridis_d(option = "D", begin = 0.25, end = 0.75) +
  facet_wrap(. ~ site) +
  theme_classic() +
  labs(
    x = parse(text = "Days~to~reach~15~degree*C~temperature"),
    y = "Days to reach 50% maximum\nlatent heat flux",
    col = "Site"
  ) +
  guides(col = "none")

p_corr2 <- df_SOS %>%
  select(-date, -variable_unit) %>%
  pivot_wider(names_from = variable, values_from = doy) %>%
  ggplot(aes(x = le, y = gpp, col = site)) +
  geom_point() +
  geom_smooth(method = MASS::rlm, se = T) +
  ggpubr::stat_cor(aes(x = le, y = gpp), label.x.npc = 0.1, label.y.npc = 0.05, show.legend = F, inherit.aes = F) +
  facet_wrap(. ~ site) +
  scale_color_viridis_d(option = "D", begin = 0.25, end = 0.75) +
  theme_classic() +
  labs(
    x = "Days to reach 50% maximum latent heat flux",
    y = "Days to reach 50% maximum\ngross primary productivity",
    col = "Site"
  ) +
  guides(col = "none")

dir.create(paste0(path, "output/"))
pdf(paste0(path, "output/ts_corr.pdf"), width = 8, height = 8)
gridExtra::grid.arrange(ggpubr::annotate_figure(p_ts, fig.lab = "(a)"),
  ggpubr::annotate_figure(p_corr1, fig.lab = "(b)"),
  ggpubr::annotate_figure(p_corr2, fig.lab = "(c)"),
  layout_matrix = rbind(
    c(1, 1),
    c(2, 3)
  ),
  heights = c(1.5, 1)
)
dev.off()
