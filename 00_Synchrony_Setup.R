#Synchrony Functions
library(corrplot)
library(cowplot)
library(lubridate)
library(dplyr)
library(DT)
# library(basemaps)
# library(forecast)
library(ggfortify)
library(gghighlight)
library(ggplot2)
library(ggplotify)
library(ggpubr)
#install.packages("greenbrown", repos="http://R-Forge.R-project.org")
# library(greenbrown)
library(grid)
library(Kendall)
library(ncf)
library(RColorBrewer)
library(reshape2)
library(sf)
library(tidyverse)
# library(timetk)
library(tseries)
library(TTR)
library(xts)
library(zoo)
#Some sources: 
#https://towardsdatascience.com/four-ways-to-quantify-synchrony-between-time-series-data-b99136c4a9c9
#https://docs.google.com/document/d/1f0FQB0Kl82hMYax6CAU4UW6ZXkprNiU_V_1PipQVrEM/edit?disco=AAAAQeymShM&usp_dm=false&pli=1

ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y =names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("R² = ",round(signif(summary(fit)$r.squared, 5), 2),
                       "Intercept =",round(signif(fit$coef[[1]],5 ), 2),
                       " Slope =",round(signif(fit$coef[[2]], 5), 2),
                       " p ≈",round(signif(summary(fit)$coef[2,4], 5),2)))+
    theme_pubr(base_size=12)
}

check_corrs <- function(x,y,z, window){
  merged <- merge(x, y, by=z)
  merged <- merged%>%
    dplyr::filter(.[[1]] < as.Date("2013-01-01"))
  
  names(merged) <- c("time", "NE2_GPP", "NE2_Tmx", "NE2_Tmn","NE2_Tpp", "NE2_Dpp", "NE2_DTx", "NE2_DTn","NE3_GPP", "NE3_Tmx", "NE3_Tmn","NE3_Tpp", "NE3_Dpp", "NE3_DTx", "NE3_DTn")
  
  merged$Pdf <- merged$NE2_Tpp - merged$NE3_Tpp
  merged$Ddf <- merged$NE2_Dpp - merged$NE3_Dpp
  merged$Gdf <- merged$NE2_GPP - merged$NE3_GPP
  melted <- melt(merged, id="time")
  melted$site <- str_sub(melted$variable, 1,3)
  melted$site[melted$site==NE_2] <- "US-Ne2.IR"
  melted$site[melted$site==NE_3] <- "US-Ne3.RF"
  melted$var <- str_sub(melted$variable, -3,-1)
  
  merge_gs <- merged %>%
    mutate(month2=as.integer(month(time))) %>%
    filter(month2 > 4  & month2 < 10)
  
  x_scale <- merged %>%
    filter(complete.cases(merged))  %>%
    mutate_at(c("NE2_GPP", "NE3_GPP", "NE2_DTx", "NE2_Dpp"), scale)
  
  x_ggplot <- melt(x_scale, id.vars="time", measure.vars=c("NE2_GPP", "NE3_GPP", "NE2_DTx", "NE2_Dpp")) 
  
  c1 <- round(cor(merged$NE2_GPP, merged$NE2_Dpp, use="complete.obs"),2)
  c2 <- round(cor(merged$NE3_GPP, merged$NE2_Dpp, use="complete.obs"),2)
  c3 <- round(cor(merged$NE2_GPP, merged$NE2_DTx, use="complete.obs"),2)
  c4 <- round(cor(merged$NE3_GPP, merged$NE2_DTx, use="complete.obs"),2)
  
  rolling <- cbind.data.frame(merged$time, runCor(merged$NE2_GPP, merged$NE2_Dpp, n=window, use="complete.obs"))
  colnames(rolling)[1:2] <- c("date", "rolling")
  rolling$month <- month(as.Date(rolling$date))
  rolling$year <- year(as.Date(rolling$date))
  
  rolling2 <- cbind.data.frame(merged$time, runCor(merged$NE3_GPP, merged$NE2_Dpp, n=window, use="complete.obs"))
  colnames(rolling2)[1:2] <- c("date", "rolling")
  rolling2$month <- month(as.Date(rolling2$date))
  rolling2$year <- year(as.Date(rolling2$date))
  
  p1 <- ggplot(subset(melted, var=="Tpp"), aes(x=time, y=value, group=site, color=site))+
    geom_line()+
    theme_minimal()+
    ylab("Unscaled tower precip")+
    xlab("time")+
    theme_pubr(base_size=10)
  
  p2 <- ggplot(subset(melted, var=="Dpp"), aes(x=time, y=value, group=site, color=site))+
    geom_line()+
    theme_minimal()+
    ylab("Unscaled daymet precip")+
    xlab("time")+
    theme_pubr(base_size=10)
  
  
  p3 <- ggplot(subset(melted, var=="Pdf"), aes(x=time, y=value))+
    geom_line()+
    theme_minimal()+
    geom_hline(yintercept=0)+
    ylab("Abs Diff Tower Precip: ")+
    xlab("time")+
    theme_pubr(base_size=10)
  
  
  p4 <- ggplot(subset(melted, var=="Ddf"), aes(x=time, y=value))+
    geom_line()+
    theme_minimal()+
    geom_hline(yintercept=0)+
    ylab("Abs Diff Daymet Precip: ")+
    xlab("time")+
    theme_pubr(base_size=10)
  
  
  p5 <- ggplot(subset(melted, var=="Gdf"), aes(x=time, y=value))+
    geom_line()+
    theme_minimal()+
    geom_hline(yintercept=0)+
    ylab("Abs Diff GPP: ")+
    xlab("time")+
    theme_pubr(base_size=10)
  
  col_order <- c("time", "NE2_GPP", "NE3_GPP", "NE2_Dpp", "NE2_DTx", "NE2_DTn", "NE3_Dpp", "NE3_DTx", "NE3_DTn","NE2_Tmx", "NE2_Tmn","NE2_Tpp", "NE3_Tmx", "NE3_Tmn","NE3_Tpp")
  
  merge_cor <- merged[ , col_order]
  
  #Quiet corplot
  #p6 <- corrplot(cor(merge_cor[,-1], method = "pearson", use = "complete.obs"), method="color", diag=FALSE, type="upper", insig="blank", addCoef.col = "black", mar=c(0,0,1,0))
  
  
  p7 <- ggplot(rolling, aes(x=date, y=rolling))+
    geom_line()+ 
    theme_pubr()+ 
    ylab("Correlation - r")+
    xlab("Time")+
    theme_pubr(base_size=10)+
    ggtitle(paste0("Rolling Correlation NE-2"))
  
  p8 <- ggplot(rolling2, aes(x=date, y=rolling))+
    geom_line()+ 
    theme_pubr()+ 
    ylab("Correlation - r")+
    xlab("Time")+
    theme_pubr(base_size=10)+
    ggtitle(paste0("Rolling Correlation NE-3"))
  
  f2 <- lm(Pdf~NE2_DTx, data=merged)
  p9 <- ggplotRegression(f2)

  col_order <- c("time", "NE2_GPP", "NE3_GPP", "NE2_Dpp", "NE2_DTx", "NE2_DTn", "NE3_Dpp", "NE3_DTx", "NE3_DTn","NE2_Tmx", "NE2_Tmn","NE2_Tpp", "NE3_Tmx", "NE3_Tmn","NE3_Tpp")
  
  #merge_corgs <- merge_gs[ , col_order]  
  #p9 <- corrplot(cor(merge_corgs[,-c(1)], method = "pearson", use = "complete.obs"),method="color", diag=FALSE, type="upper", insig="blank", addCoef.col = "black")
  
  p10 <- ggCcf(merged$NE2_GPP, merged$NE2_Dpp, ylab="Cross-correlation NE2", na.action=na.pass)+ ggtitle(NULL)
  
  p11 <- ggCcf(merged$NE3_GPP, merged$NE3_Dpp, ylab="Cross-correlation NE3", na.action=na.pass)+ ggtitle(NULL)
  
  p12 <- ggplot(subset(x_ggplot, variable !="NE2_DTx"), aes(x=time, y=value, group=variable, color=variable))+
    geom_line(size=1)+
    scale_color_manual(values=c("#FC4D08", "#00AFBB", "#E7B802"))+
    theme_pubr()+ 
    ylab("Scaled Monthly Values")+
    xlab("Time")+
    ylim(-3, 3)+
    ggtitle(paste0("Correlation (r) NE-2:Precip =  ", c1, "   NE-3:Precip = ", c2))+
    theme(legend.position = "bottom")
  
  
  p13 <- ggplot(subset(x_ggplot, variable !="NE2_Dpp"), aes(x=time, y=value, group=variable, color=variable))+
    geom_line(size=1)+
    scale_color_manual(values=c("#FC4D08", "#00AFBB", "#E7B802"))+
    theme_pubr()+ 
    ylab("Scaled Monthly Values")+
    xlab("Time")+
    ylim(-3, 3)+
    ggtitle(paste0("Correlation (r) NE-2:Tmax =  ", c3, "   NE-3:Tmax = ", c4))+
    theme(legend.position = "bottom")
  
  gridExtra::grid.arrange(p1,p2,p3,p4,p12,p13, top = textGrob("NE2 = Irrigated, NE-3 = Rainfed"))
  #gridExtra::grid.arrange(p7,p8)
  #gridExtra::grid.arrange(p10,p11)
}

check_corrs_year <- function(x,y,z,window){
  merged <- merge(x, y, by=z)
  merged <- merged%>%
    transform(year = seq(2001, 2020, by=1))%>%
    dplyr::filter(year < 2013) %>%
    transform(year=as.integer(year))
  
  names(merged) <- c("time", "NE2_GPP", "NE2_Tmx", "NE2_Tmn","NE2_Tpp", "NE2_Dpp", "NE2_DTx", "NE2_DTn","NE3_GPP", "NE3_Tmx", "NE3_Tmn","NE3_Tpp", "NE3_Dpp", "NE3_DTx", "NE3_DTn")
  
  #PDF - difference between tower precip between sites
  merged$Pdf <- merged$NE2_Tpp - merged$NE3_Tpp
  
  #DDF - difference between daymet precip between sites
  merged$Ddf <- merged$NE2_Dpp - merged$NE3_Dpp
  
  #GDF - Difference between GPP between sites
  merged$Gdf <- merged$NE2_GPP - merged$NE3_GPP
  
  
  x_scale <- merged %>%
    filter(complete.cases(merged))  %>%
    mutate_at(c("NE2_GPP", "NE3_GPP", "NE2_DTx", "NE2_Dpp"), scale)
  
  x_ggplot <- melt(x_scale, id.vars="time", measure.vars=c("NE2_GPP", "NE3_GPP", "NE2_DTx", "NE2_Dpp")) 
  
  c1 <- round(cor(merged$NE2_GPP, merged$NE2_Dpp, use="complete.obs"),2)
  c2 <- round(cor(merged$NE3_GPP, merged$NE2_Dpp, use="complete.obs"),2)
  c3 <- round(cor(merged$NE2_GPP, merged$NE2_DTx, use="complete.obs"),2)
  c4 <- round(cor(merged$NE3_GPP, merged$NE2_DTx, use="complete.obs"),2)
  
  #Melt -> mostly for plotting
  melted <- melt(merged, id="time")
  melted$site <- str_sub(melted$variable, 1,3)
  melted$site[melted$site==NE_2] <- "US-Ne2.IR"
  melted$site[melted$site==NE_3] <- "US-Ne3.RF"
  melted$var <- str_sub(melted$variable, -3,-1)
  melted$time <- as.Date(paste(melted$time, "01", "01", sep="-"))

  #Rolling correlation for NE2 using the merged data 
  rolling <- cbind.data.frame(merged$time, runCor(merged$NE2_GPP, merged$NE2_Dpp, n=window, use="complete.obs"))
  colnames(rolling)[1:2] <- c("date", "rolling")
  rolling$month <- month(as.Date(rolling$date))
  rolling$year <- year(as.Date(rolling$date))
  
  #Rolling correlation for NE3 using the merged data 
  rolling2 <- cbind.data.frame(merged$time, runCor(merged$NE3_GPP, merged$NE3_Dpp, n=window, use="complete.obs"))
  colnames(rolling2)[1:2] <- c("date", "rolling")
  rolling2$month <- month(as.Date(rolling2$date))
  rolling2$year <- year(as.Date(rolling2$date))
  
  p1 <- ggplot(subset(melted, var=="Tpp"), aes(x=time, y=value, group=site, color=site))+
    geom_line()+
    theme_minimal()+
    ylab("Unscaled tower precip")+
    scale_colour_manual(values=c("#FC4D08", "#00AFBB"))+
    xlab("time")+
    theme_pubr(base_size=10)
  
  p2 <- ggplot(subset(melted, var=="Dpp"), aes(x=time, y=value, group=site, color=site))+
    geom_line()+
    theme_minimal()+
    ylab("Unscaled daymet precip")+
    scale_colour_manual(values=c("#FC4D08", "#00AFBB"))+
    xlab("time")+
    theme_pubr(base_size=10)
  
  p3 <- ggplot(subset(melted, var=="Pdf"), aes(x=time, y=value))+
    geom_line()+
    theme_minimal()+
    geom_hline(yintercept=0)+
    ylab("Abs Diff Tower Precip: ")+
    xlab("time")+
    ylim(-55,450)+
    theme_pubr(base_size=10)
  
  p4 <- ggplot(subset(melted, var=="Ddf"), aes(x=time, y=value))+
    geom_line()+
    theme_minimal()+
    geom_hline(yintercept=0)+
    ylab("Abs Diff Daymet Precip: ")+
    xlab("time")+
    ylim(-55,450)+
    theme_pubr(base_size=10)
  
  p5 <- ggplot(subset(melted, var=="Gdf"), aes(x=time, y=value))+
    geom_line()+
    theme_minimal()+
    geom_hline(yintercept=0)+
    ylab("Abs Diff GPP:")+
    xlab("time")+
    theme_pubr(base_size=10)
  
  col_order <- c("time", "NE2_GPP", "NE3_GPP", "NE2_Dpp", "NE2_DTx", "NE2_DTn", "NE3_Dpp", "NE3_DTx", "NE3_DTn","NE2_Tmx", "NE2_Tmn","NE2_Tpp", "NE3_Tmx", "NE3_Tmn","NE3_Tpp")
  
  merge_cor <- merged[ , col_order]
  #Quiet Corplot
  #p6 <- corrplot(cor(merge_cor[,-c(1)], method = "pearson", use =   #"complete.obs"),method="color", diag=FALSE, type="upper", insig="blank", addCoef.col = #"black")
  
  p7 <- ggplot(rolling, aes(x=date, y=rolling))+
    geom_line()+ 
    theme_pubr()+ 
    ylab("Correlation - r")+
    xlab("Time")+
    theme_pubr(base_size=10)+
    ggtitle(paste0("Rolling Correlation NE-2"))
  
  p8 <- ggplot(rolling2, aes(x=date, y=rolling))+
    geom_line()+ 
    theme_pubr()+ 
    ylab("Correlation - r")+
    xlab("Time")+
    theme_pubr(base_size=10)+
    ggtitle(paste0("Rolling Correlation NE-3"))
  
  p9 <- ggCcf(merged$NE2_GPP, merged$NE2_Dpp, ylab="Cross-correlation NE2", na.action=na.pass)+ ggtitle(NULL)
  
  p10 <- ggCcf(merged$NE3_GPP, merged$NE3_Dpp, ylab="Cross-correlation NE3", na.action=na.pass)+ ggtitle(NULL)
  

  p11 <- ggplotRegression( lm(Pdf~NE2_DTx, data=merged))
  p12 <- ggplotRegression( lm(Pdf~NE2_Dpp, data=merged))
  
  p13 <- ggplot(subset(x_ggplot, variable !="NE2_DTx"), aes(x=time, y=value, group=variable, color=variable))+
    geom_line(size=1)+
    scale_color_manual(values=c("#FC4D08", "#00AFBB", "#E7B802"))+
    theme_pubr()+ 
    ylab("Scaled Annual Values")+
    xlab("Time")+
    ylim(-3, 3)+
    ggtitle(paste0("Correlation (r) US-Ne2 IR:Precip =  ", c1, "   US-Ne3 RF:Precip = ", c2))+
    theme(legend.position = "bottom")
  
  print(str(x_ggplot))
  p14 <- ggplot(subset(x_ggplot, variable !="NE2_Dpp"), aes(x=time, y=value, group=variable, color=variable))+
    geom_line(size=1)+
    scale_color_manual(values=c("#FC4D08", "#00AFBB", "#E7B802"))+
    theme_pubr()+ 
    ylab("Scaled Annual Values")+
    xlab("Time")+
    ylim(-3, 3)+
    ggtitle(paste0("Correlation (r) NE-2:Tmax = ", c3, "   NE-3:Tmax =  ", c4))+
    theme(legend.position = "bottom")
  

  
  gridExtra::grid.arrange(p1,p2,p3,p4, top = textGrob("US-Ne2 IR - US-Ne3 RF"))
  gridExtra::grid.arrange(p11,p12, ncol=2)
  gridExtra::grid.arrange(p13,p14, ncol=2)
}

get_pheno_NE2<- function(x){
  y <- as.data.frame(t((PhenoDeriv(TsPP(NE2_ts, interpolate=TRUE)[(times$start[x]):(times$end[x])]))))
  y$year <- as.Date(paste((x + 2000), "01", "01", sep="-"))
  return(y)
}

get_pheno_NE3<- function(x){
  y <- as.data.frame(t((PhenoDeriv(TsPP(NE3_ts, interpolate=TRUE)[(times$start[x]):(times$end[x])]))))
  y$year <- as.Date(paste((x + 2000), "01", "01", sep="-"))
  return(y)
}

get_fluxpheno_NE2<- function(x){
  y <- as.data.frame(t((PhenoDeriv(TsPP(NE2_ts_flux, interpolate=TRUE)[(times$start[x]):(times$end[x])]))))
  y$year <- as.Date(paste((x + 2000), "01", "01", sep="-"))
  return(y)
}

get_fluxpheno_NE3<- function(x){
  y <- as.data.frame(t((PhenoDeriv(TsPP(NE3_ts_flux, interpolate=TRUE)[(times$start[x]):(times$end[x])]))))
  y$year <- as.Date(paste((x + 2000), "01", "01", sep="-"))
  return(y)
}


PlotPhenCycle <- structure(function(
  ##title<< 
  ## Plot a easonal cycle with phenology metrics
  ##description<<
  ## This function plots a seasonal cycle with phenology metrics.
  
  x, 
  ### values of one year
  
  xpred = NULL,
  ### smoothed/predicted values 
  
  metrics,
  ### vector of pheology metrics
  
  xlab="DOY",
  ### label for x-axis
  
  ylab="EVI",
  ### label for y-axis
  
  trs = NULL,
  ### threshold for threshold methods
  
  main="",
  ### title
  
  ...
  ### further arguments (currently not used)
  
  ##seealso<<
  ## \code{\link{Phenology}}
  
) {
  
  sos <- metrics[1]
  eos <- metrics[2]
  los <- metrics[3]
  pop <- metrics[4]
  pot <- metrics[5]
  mgs <- metrics[6]
  rsp <- metrics[7]
  rau <- metrics[8]
  peak <- metrics[9]
  trough <- metrics[10]
  msp <- metrics[11]
  mau <- metrics[12]
  
  cols <- c("blue", "darkgreen", "darkorange")
  ylim <- range(c(x, xpred), na.rm=TRUE)
  ylim[2] <- 0.75
  # plot time series
  t <- 1:length(x)
  plot(t, x, xlab=xlab, ylab=ylab, type="l", col="darkgrey", ylim=ylim, main=main, ...)
  if (!is.null(xpred)) lines(1:length(x), xpred, col="red")
  
  # plot SOS, EOS
  abline(v=c(sos, eos), col=cols[1:2], lty=2)
  text(x=c(sos, eos), y=min(ylim), paste(c("SOS", "EOS"), "=", round(c(sos, eos), 0)), col=cols[1:2], pos=3) 
  
  #plot POP
  abline(v=pop, col=cols[3], lty=2)
  #text(x=pop, y=min(ylim)+diff(ylim)*0.1, paste("POP", "=", round(pop, 0)), col=cols[3], pos=3) 
  text(x=pop, y=(peak+0.07), paste("PEAK", "=", signif(peak, 2)), pos=1)
  
  # plot POT
  #abline(v=pot, col=cols[3], lty=2)
  #text(x=pot, y=min(ylim)+diff(ylim)*0.1, paste("POT", "=", round(pot, 0)), col=cols[3], pos=3) 
  #text(x=pot, y=trough, paste("TROUGH", "=", signif(trough, 2)), pos=3)
  
  # # plot treshold
  # if (!is.null(trs)) {
  # abline(h=trs, lty=2)
  # text(x=1, y=trs, paste("trs =", signif(trs, 2)), pos=4)
  # }
  
  # plot LOS
  text(x=pop, y=msp*0.7, paste(c("LOS"), "=", round(los, 0)), pos=1, col=cols[3]) 
  if (!is.na(eos) & !is.na(sos)) {
    if (sos < eos) segments(x0=sos, y0=msp*0.7, x1=eos, y1=msp*0.7, col=cols[3])
    if (sos > eos) segments(x0=c(0, sos), y0=c(mgs, mgs), x1=c(eos, length(t)), y1=c(mgs, mgs), col=cols[3])
  }
  
  # plot MSP and MAU averages
  segments(x0=c(sos-10, eos-10), y0=c(msp, mau), x1=c(sos+10, eos+10), y1=c(msp, mau), col=cols[1:2])	
  text(x=c(sos, eos, pop), y=c(msp, mau, msp*0.7), paste(c("MSP", "MAU", "MGS"), "=", signif(c(msp, mau, mgs), 2)), pos=3, col=cols)
  
  # plot spring and autumn rates
  #if (!is.na(rsp) | !is.na(rau)) {
  #text(x=c(sos, eos), y=c(msp, mau), paste(c("RSP", "RAU"), "=", signif(c(rsp, rau), 2)), col=cols, pos=1) 
  
  #}
  
}, ex=function() {
  
  # perform time series preprocessing for first year of data
  x <- TsPP(ndvi, interpolate=TRUE)[1:365]
  plot(x)
  
  # calculate phenology metrics for first year
  metrics <- PhenoTrs(x, approach="White")
  PlotPhenCycle(x, metrics=metrics)
  
})

PlotPhenCycle2 <- structure(function(
  ##title<< 
  ## Plot a easonal cycle with phenology metrics
  ##description<<
  ## This function plots a seasonal cycle with phenology metrics.
  
  x, 
  ### values of one year
  
  xpred = NULL,
  ### smoothed/predicted values 
  
  metrics,
  ### vector of pheology metrics
  
  xlab="DOY",
  ### label for x-axis
  
  ylab="GPP",
  ### label for y-axis
  
  trs = NULL,
  ### threshold for threshold methods
  
  main="",
  ### title
  
  ...
  ### further arguments (currently not used)
  
  ##seealso<<
  ## \code{\link{Phenology}}
  
) {
  
  sos <- metrics[1]
  eos <- metrics[2]
  los <- metrics[3]
  pop <- metrics[4]
  pot <- metrics[5]
  mgs <- metrics[6]
  rsp <- metrics[7]
  rau <- metrics[8]
  peak <- metrics[9]
  trough <- metrics[10]
  msp <- metrics[11]
  mau <- metrics[12]
  
  cols <- c("blue", "darkgreen", "darkorange")
  ylim <- range(c(x, xpred), na.rm=TRUE)
  ylim[2] <- 27
  # plot time series
  t <- 1:length(x)
  plot(t, x, xlab=xlab, ylab=ylab, type="l", col="darkgrey", ylim=ylim, main=main, ...)
  if (!is.null(xpred)) lines(1:length(x), xpred, col="red")
  
  # plot SOS, EOS
  abline(v=c(sos, eos), col=cols[1:2], lty=2)
  text(x=c(sos, eos), y=min(ylim), paste(c("SOS", "EOS"), "=", round(c(sos, eos), 0)), col=cols[1:2], pos=3) 
  
  #plot POP
  abline(v=pop, col=cols[3], lty=2)
  #text(x=pop, y=min(ylim)+diff(ylim)*0.1, paste("POP", "=", round(pop, 0)), col=cols[3], pos=3) 
  text(x=pop, y=max(peak+2), paste("PEAK", "=", signif(peak, 2)), pos=1)
  
  # plot POT
  #abline(v=pot, col=cols[3], lty=2)
  #text(x=pot, y=min(ylim)+diff(ylim)*0.1, paste("POT", "=", round(pot, 0)), col=cols[3], pos=3) 
  #text(x=pot, y=trough, paste("TROUGH", "=", signif(trough, 2)), pos=3)
  
  # # plot treshold
  # if (!is.null(trs)) {
  # abline(h=trs, lty=2)
  # text(x=1, y=trs, paste("trs =", signif(trs, 2)), pos=4)
  # }
  
  # plot LOS
  text(x=pop, y=mgs, paste(c("LOS"), "=", round(los, 0)), pos=1, col=cols[3]) 
  if (!is.na(eos) & !is.na(sos)) {
    if (sos < eos) segments(x0=sos, y0=mgs, x1=eos, y1=mgs, col=cols[3])
    if (sos > eos) segments(x0=c(0, sos), y0=c(mgs, mgs), x1=c(eos, length(t)), y1=c(mgs, mgs), col=cols[3])
  }
  
  # plot MSP and MAU averages
  segments(x0=c(sos-10, eos-10), y0=c(msp, mau), x1=c(sos+10, eos+10), y1=c(msp, mau), col=cols[1:2])	
  text(x=c(sos, eos, pop), y=c(msp, mau, mgs), paste(c("MSP", "MAU", "MGS"), "=", signif(c(msp, mau, mgs), 2)), pos=3, col=cols)
  
  # plot spring and autumn rates
  #if (!is.na(rsp) | !is.na(rau)) {
  #text(x=c(sos, eos), y=c(msp, mau), paste(c("RSP", "RAU"), "=", signif(c(rsp, rau), 2)), col=cols, pos=1) 
  
  #}
  
}, ex=function() {
  
  # perform time series preprocessing for first year of data
  x <- TsPP(ndvi, interpolate=TRUE)[1:365]
  plot(x)
  
  # calculate phenology metrics for first year
  metrics <- PhenoTrs(x, approach="White")
  PlotPhenCycle2(x, metrics=metrics)
  
})
Pheno_wrap <- function(x,y,z){
  q <- as.numeric(y)-2000
  x2 <- TsPP(x, interpolate=TRUE)[(times$start[q]):(times$end[q])]
  PlotPhenCycle(x2, metrics=PhenoDeriv(x2, approach="White"), main=paste(y, z, sep = "-"))

}

Pheno_wrap2 <- function(x,y){
  q <- as.numeric(y)-2000
  x2 <- TsPP(x, interpolate=TRUE)[(times$start[q]):(times$end[q])]
  PlotPhenCycle2(x2, metrics=PhenoDeriv(x2, approach="White"))

}

RobustMax <- function(x) {if (length(x)>0) max(x) else -Inf}

RobustMin <- function(x) {if (length(x)>0) min(x) else -Inf}

ggplotRegression2 <- function (x,y,a,b,c) {
  require(ggplot2)
  fit <- lm(y~x)
  #print(summary(fit))

regplot <- ggplot(fit$model, aes_string(x = names(fit$model)[2], y =names(fit$model)[1])) +
  geom_point() +
  theme_pubr(base_size=12)+
  theme(plot.title = element_text(size=8))+
  xlab(a)+
  ylab(b)+
  ylim(min(fit$model$y)*0.95, max(fit$model$y)*1.05)+
  xlim(min(fit$model$x)*0.95, max(fit$model$x)*1.05)+
  annotate("text", -Inf, Inf, label = c, hjust = 0, vjust = 1)+
  ggtitle(paste("R² = ",round(signif(summary(fit)$r.squared, 5), 2),
                "Intercept =",round(signif(fit$coef[[1]],5 ), 1),
                " P ≈",round(signif(summary(fit)$coef[2,4], 5),2)))

if(summary(fit)$coef[2,4] < 0.05){
  regplot <- regplot + stat_smooth(method = "lm", col = "red")}
if(summary(fit)$coef[2,4] > 0.05 & summary(fit)$coef[2,4]<0.1){
  regplot <- regplot + stat_smooth(method = "lm", col = "red", linetype="dashed")}
else{regplot <- regplot
}
regplot    
}

ggregression_wrap <- function(x,y,z,a,b){
  t <- NE2_ann %>%
    dplyr::select(x,y)
  
  u <- NE3_ann %>%
    dplyr::select(x,y)
  
  ggarrange(ggplotRegression2(t[,2], t[,1], b=b,a=a,paste("a) NE-2", z, sep="-")), ggplotRegression2(u[,2], u[,1], b=b,a=a,paste("b) NE-3", z, sep="-")))
}

ggregression_wrap2 <- function(x,y,z,a,b){
  t <- NE2_gs %>%
    dplyr::select(x,y)
  
  u <- NE3_gs %>%
    dplyr::select(x,y)
  
  ggarrange(ggplotRegression2(t[,2], t[,1], b=b,a=a,paste("NE2 Gs", z, sep="-")), ggplotRegression2(u[,2], u[,1], b=x,a=y,paste("NE3 Gs", z, sep="-")))
  
}

ggregression_wrap3 <- function(x,y,z){
  t <- NE2_ann %>%
    dplyr::select(x,y)
  
  u <- NE3_ann %>%
    dplyr::select(x,y)
  
  ggarrange(ggplotRegression2(t[,2], t[,1], b=x,a=y,paste("a) NE-2 Annual", z, sep="-")), ggplotRegression2(u[,2], u[,1], b=x,a=y,paste("NE3 Annual", z, sep="-")))
}

ggregression_wrap4 <- function(x,y,z){
  t <- NE2_gs_flux %>%
    dplyr::select(x,y)
  
  u <- NE3_gs_flux %>%
    dplyr::select(x,y)
  
  ggarrange(ggplotRegression2(t[,2], t[,1], b=x,a=y,paste("NE2 Gs", z, sep="-")), ggplotRegression2(u[,2], u[,1], b=x,a=y,paste("NE3 Gs", z, sep="-")))
  
}

ggregression_wrap5 <- function(x,y,z,a,b){
  t <- Corrob_rs_NE2 %>%
    dplyr::select(x,y)
  
  u <- Corrob_rs_NE3 %>%
    dplyr::select(x,y)
  
  ggarrange(ggplotRegression2(t[,2], t[,1], b=b,a=a,paste("a) US-Ne2 IR", z, sep="-")), ggplotRegression2(u[,2], u[,1], b=b,a=a, paste("b) US-Ne3 RF", z, sep="-")))
  
}

Process_RS <- function(x,y,date){
  #x=filename of RS timeseries
  #y=variable name 
  Neb <- read.csv(x) %>%
    dplyr::select(ID, Date, y)%>%
    mutate(Date=as.Date(Date))%>%
    filter(Date > date)%>%
    rename("y"=y)
  
  NE2_NDVI <- Neb[Neb$ID=="US-NE2",]
  NE3_NDVI <- Neb[Neb$ID=="US-NE3",]
  
  NE2_NDVI$year <- as.character(year(as.Date(NE2_NDVI$Date)))
  NE3_NDVI$year <- as.character(year(as.Date(NE3_NDVI$Date)))

  return(list(NE2_NDVI, NE3_NDVI))}

f <- function(bar) bar[which.max(bar[,1]),2]
