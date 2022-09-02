#### Synchron working group #####

# Travis McDevitt-Galles
# 09/24/2021
# title: 02_Synchrony_GAM_Function

# Functions to run to iterate through species, year and sites to extract
# fitted values to quantify phenological overlap
# Modified from Mosquito GAM project to extract phenometrics
 
pheno_gam <- function(data, Site){
  
  foc.df <- filter( data , Site == Site  )
  
  t <- 1

  plot_gam <- list() # empty list to store gam plots
  
  for(y in 1:length( unique( foc.df$Year ) ) ){
    
    year.df <- foc.df %>% filter( Year == unique( foc.df$Year )[y])

    for( j in 1:length( unique( year.df$SciName ) ) ){
      
      gam.df <- year.df %>%  filter(  SciName == unique( year.df$SciName )[j])
      
      gam.df <- gam.df %>% group_by( SciName, Site, Plot, DOY, Year)  %>% 
        summarise(
          Count =  sum( Count ),
          TrapHours = sum( TrapHours ),
          TotalWeight = sum( TotalWeight ),
          SubsetWeight = sum( SubsetWeight )
        ) %>% ungroup()
      
      offset <- gam.df$SubsetWeight/gam.df$TotalWeight * (gam.df$TrapHours/24)
      
      if( length( unique( gam.df$DOY ) ) > 6 ){
        
        
        bayGam <- stan_gamm4( Count ~ s( DOY, bs = "cc", k =6) +
                              offset( log( offset ) ) ,
                              #random= ~(1|Plot) ,
                              data = gam.df , family = 'poisson',
                              chains = 4, iter = 3000 ,
                              control = list( adapt_delta = 0.99 )
                            )
        
        ## creating a new data frame
        DOY <- seq(range(gam.df$DOY)[1],range(gam.df$DOY)[2], by = 1 )
        
        new.data <- as.data.frame(DOY)
      
        new.data$offset <- offset1 <- rep(.75, nrow(new.data))
        
        
        
        dfit  <- posterior_epred(bayGam, newdata = new.data, 
                                 offset = log(offset1))
        
        plot.df <- as.data.frame(matrixStats::colMedians(dfit))
        
        plot.df$DOY <- new.data$DOY
        
        hold <- matrixStats::colQuantiles(dfit, probs=c(0.025, 0.975))
        
        plot.df$Q_2.5 <- hold[,1]
        
        plot.df$Q_97.5 <- hold[,2]
        
        colnames(plot.df)[1] <- "Count"
        
        ## Plotting best fit ##
        p1 <- ggplot() +
          geom_line( data = plot.df, aes( x = DOY, y = Count ) )+
          geom_ribbon( data = plot.df, aes( x = DOY, ymin = Q_2.5,
                                            ymax = Q_97.5 ),
                       alpha = .5, color = "grey") + theme_classic() +
          geom_point( data = gam.df, aes( x = DOY, y = Count*offs ) ) +
          xlab( paste( "DOY — ", unique( gam.df$Year[1] ), "—", 
                     unique( gam.df$Site )[1] ) ) +
          ylab( paste( "Count — ", unique( gam.df$SciName[1] ) ) ) +
          theme(
            axis.line.x = element_line( color = "black" ) ,
            axis.ticks.y = element_line( color = "black" ),
            axis.ticks.x = element_line( color = "black" ),
            axis.title.x = element_text( size = rel( 1.8 ) ),
            axis.text.x = element_text( vjust = 0.5, color = "black" ),
            axis.text.y = element_text( vjust = 0.5, color = "black" ),
            axis.title.y = element_text( size = rel( 1.8 ), angle = 90 ),
            strip.text.x = element_text( size = 20 ) )
        
        plot_gam[[t]] <- p1  # add each plot into plot list
        
        ## Saving predited values 
        pred.df <- plot.df 
        # Adding new columns
        pred.df$Year <- unique( gam.df$Year )
        pred.df$Site <- unique( gam.df$Site )
        pred.df$SciName <- unique( gam.df$SciName )
        
        if( t == 1 ){
          
          out.df <- pred.df
          
          t <- t + 1
          
        }else{
          
          out.df <- rbind.data.frame( out.df, pred.df )
          
          t <- t + 1
          
        }
      }
    }
  }
  return( list( out.df, plot_gam ) )
}

