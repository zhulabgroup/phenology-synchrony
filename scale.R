library(tidyverse)

scale_dict<-read_csv("./synchrony scale.csv", skip=0,n_max = 5)

scale_df<-read_csv("./synchrony scale.csv", skip=8, col_select = 1:4) %>% 
  mutate(ss=ss*2,
         ts=ts*2)

scale_sum<-scale_df %>% 
  group_by(Name) %>% 
  summarise(ss_r=max(ss)-min(ss)+2,
            ts_r=max(ts)-min(ts)+2,
            loo_r=max(loo)-min(loo)+1,
            ss_m=mean(ss),
            ts_m=mean(ts),
            loo_m=mean(loo)
            ) %>% 
  ungroup()
scale_sum

theta <- seq(-pi/2, pi/2, by=0.1)
phi <- seq(0, 2*pi, by=0.1)
mgrd <- meshgrid(phi, theta)
phi <- mgrd$X
theta <-  mgrd$Y
ell_df_list<-vector(mode="list")
for (i in 1:nrow(scale_sum)) {
  x <- cos(theta) * cos(phi) * scale_sum$ss_r[i]/2 +scale_sum$ss_m[i]
  dim(x)<-NULL
  y <- cos(theta) * sin(phi) * scale_sum$ts_r[i]/2 +scale_sum$ts_m[i]
  dim(y)<-NULL
  z <- sin(theta) *  scale_sum$loo_r[i]/2 +scale_sum$loo_m[i]
  dim(z)<-NULL
  ell_df_list[[i]]<-data_frame(x=x, y=y, z=z, name=scale_sum$Name[i])
}
ell_df<-bind_rows(ell_df_list)

library(plotly)
# https://plotly.com/r/3d-axes/
p<-plot_ly() %>%
  add_mesh(data=ell_df,
            x = ~x, y = ~y, z = ~z, color=~as.factor(name),
            opacity=0.01) %>% 
  add_text(data=scale_sum,
    x = ~ss_m, y = ~ts_m, z = ~loo_m-0.25, text=~as.factor(Name),
    opacity=1) %>% 
  layout(scene = list(xaxis=list(
    title = "·<b><br><br><br>spatial scale</b>",
    ticketmode = 'array',
    ticktext = c("Single-<br>location", "Multiple-<br>locations"),
    tickvals = c(2,4),
    range = c(0,5)
  ),
  yaxis=list(
    title = "·<b><br><br><br>temporal scale</b>",
    ticketmode = 'array',
    ticktext = c("Intra-<br>annual", "Inter-<br>annual"),
    tickvals = c(2,4),
    range = c(0,5)
  ),
  zaxis=list(
    title = "·<b><br><br><br>level of organization</b>",
    ticketmode = 'array',
    ticktext = c("Population", "Meta-<br>population", "Community", "Ecosystem"),
    tickvals = c(1, 2, 3, 4),
    range = c(0,5)#,
    # titlefont = list(size = 15),
    # tickfont = list(size = 15),
    # tickangle = -30
  )
  )
  )
print(p)


data(iris)
head(iris)

# x, y and z coordinates
x <- sep.l <- iris$Sepal.Length
y <- pet.l <- iris$Petal.Length
z <- sep.w <- iris$Sepal.Width

library("plot3D")
scatter3D(scale_df$ss, scale_df$ts, scale_df$loo, bty = "g", pch = 18, 
          col.var = scale_df$Name %>% as.factor() %>% as.integer(),
          phi=20, theta=30,
          cex=3
          # col = c("#1B9E77", "#D95F02", "#7570B3"),
          # pch = 18#,
          # colkey = list(at = c(1, 2, 3, 4, 5), side = 1,
          #               addlines = TRUE, length = 0.5, width = 0.5,
          #               labels = scale_df$Name %>% unique() %>% sort())
          )

text3D(scale_df$ss, scale_df$ts, scale_df$loo,  labels = scale_df$Name,
       phi=20, theta=30,
       add = TRUE, colkey = FALSE, cex = 0.5)

p <- plot_ly() %>%
  # the scatter plot of the data points 
  add_trace(x=scale_df$ss, y=scale_df$ts, z=scale_df$loo,
            type="scatter3d", mode="markers",
            # marker = list(color=scale_df$ss, 
            #               colorscale = c("#FFE1A1", "#683531"), 
            #               opacity = 0.7, size=2),
            name=scale_df$Name) %>% 
  add_trace(type = 'scatter3d', size = 1, 
            x = ellipse$vb[1,], y = ellipse$vb[2,], z = ellipse$vb[3,], 
            opacity=0.01) %>% 
  add_text(x=scale_df$ss, y=scale_df$ts, z=scale_df$loo,
            type="scatter3d", mode="markers",
            # marker = list(color=scale_df$ss, 
            #               colorscale = c("#FFE1A1", "#683531"), 
            #               opacity = 0.7, size=2),
            text=scale_df$Name) 
print(p)
