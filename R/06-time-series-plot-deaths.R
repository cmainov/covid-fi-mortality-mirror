###-----------------------------------------------------------------------------
###   06-TIME SERIES PLOT OF DEATHS STRATIFIED ON DATASET AND URBAN-RURAL STATUS
###-----------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we generate a time series plot of deaths across both datasets
# employed in this analysis (JHU and CDC) and across urban-rural status.
#
# INPUT DATA FILE: "03-Data-Rodeo/02-long-form-analytic-data.rds"
#
# OUPUT DATA FILE: "04-Tables-Figures/time-series-stratified.tiff"
#
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

library( tidyverse)
library( ggrepel )
library( latex2exp )

### 0.0 Data Import and Time-Series Data Wrangle ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


## (0.1) Import raw long-form/time-series data ##
d.time <- readRDS( "03-Data-Rodeo/02-long-form-analytic-data.rds")


## (0.2) Summarize/aggregate time series data by week and by urban-rural status ##
d.time.urban <- d.time %>%
  filter( !is.na( urb.cat.code ) ) %>%  # remove counties with missing urban-rural code (these are all counties in Puerto Rico--PR)
  group_by( week, urb.cat.code ) %>%
  mutate( week.deaths.cdc = sum( crude.deaths.cdc, na.rm = T ),  # note we are using crude and not age adjusted deaths counts for this plot
          week.deaths.jhu = sum( deaths.jhu, na.rm = T  ) ) %>%
  ungroup() %>%
  distinct( week, urb.cat.code, week.deaths.cdc, week.deaths.jhu, week.start, week.end ) %>%
  pivot_longer( cols = week.deaths.cdc:week.deaths.jhu,
                names_to = "source",
                values_to = "week.deaths" )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (1.0) Time Series Plot ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (1.1) Prepare data for plotting ##
label.data <- d.time.urban %>%
  
  # isolate max and min values of CDC curve and assign label (data) so that we can use geom_label_repel below
  mutate( line.label = ifelse( week == min( .$week ) & urb.cat.code == "Metropolitan" & source == "week.deaths.cdc",
                               paste0( week.start ),
                               ifelse( week == max( .$week ) & urb.cat.code == "Metropolitan" & source == "week.deaths.cdc",
                                       paste0(week.end), NA ) ),
          
          # nudge values for those labels
          nudge.y = ifelse( !is.na( line.label ) & week == min( .$week ),
                            300000,
                            ifelse( !is.na( line.label ) & week == max( .$week ),
                                    0, NA ) ),
          nudge.x = ifelse( !is.na( line.label ) & week == min( .$week ),
                            0,
                            ifelse( !is.na( line.label ) & week == max( .$week ),
                                    -40, NA ) ) ) %>% data.frame()


## (1.2) Custom X-axis break marks/ticks ##

levels.week <- length( levels( factor( label.data$week ) ) )

brks <- c( min( label.data$week ),
           round( levels.week / 4 ),
           round( levels.week / 2) ,
           round( levels.week * 0.75 ),
           max( label.data$week ) )


## (1.3) Values for line ranges showing differences in curves ##

# metropolitan curves
met.high.y.cdc <- label.data[ which( label.data$week == max( label.data$week ) &
                                       label.data$urb.cat.code == "Metropolitan" &
                                       label.data$source == "week.deaths.cdc"), "week.deaths"]

met.high.y.jhu <- label.data[ which( label.data$week == max( label.data$week ) &
                                       label.data$urb.cat.code == "Metropolitan" &
                                       label.data$source == "week.deaths.jhu"), "week.deaths"]

# non-metropolitan curves
nonmet.high.y.cdc <- label.data[ which( label.data$week == max( label.data$week ) &
                                          label.data$urb.cat.code == "Non-metropolitan" &
                                          label.data$source == "week.deaths.cdc"), "week.deaths"]


nonmet.high.y.jhu <- label.data[ which( label.data$week == max( label.data$week ) &
                                          label.data$urb.cat.code == "Non-metropolitan" &
                                          label.data$source == "week.deaths.jhu"), "week.deaths"]

# y-values for labeling differences between metropolitan curves
hi.pt.met <- max( met.high.y.cdc, met.high.y.jhu )
lo.pt.met <- min( met.high.y.cdc, met.high.y.jhu )

# y-values for labeling differences between non-metropolitan curves
hi.pt.nonmet <- max( nonmet.high.y.cdc, nonmet.high.y.jhu )
lo.pt.nonmet <- min( nonmet.high.y.cdc, nonmet.high.y.jhu )

# y-value position of text showing differences in deaths (taken at the midpoint between high and low values)
met.txt.pos.y <- ( (hi.pt.met-lo.pt.met) / 2 ) + lo.pt.met
nonmet.txt.pos.y <- ( (hi.pt.nonmet-lo.pt.nonmet) / 2 ) + lo.pt.nonmet

# x parameter and horizontal nudge values for death differences text
x.par <- max( label.data$week ) + 3
x.nudge <- 16


## (1.4) Generate the Plot ##

ggplot( data = label.data,
        mapping = aes( x = week, 
                       y = week.deaths, 
                       col = source, 
                       lty = urb.cat.code ) ) +
  scale_linetype_manual( name = "Urban-Rural Status", 
                         values = c( 1, 2 ) ) +
  scale_color_manual( name = "Data Source", 
                      values = c( "navyblue", "firebrick2" ), 
                      breaks = c( "week.deaths.jhu", "week.deaths.cdc" ),
                      labels = c( "Johns Hopkins CRC", "CDC" ) )+
  geom_line() +
  theme_classic() +
  scale_y_continuous( labels = ~ format( .x, big.mark = ",", scientific = FALSE ) ) + # remove scientific notation on y-axis labels and add comma to thousands place (`big.mark` option)
  ylab( "Cumulative COVID-19 Deaths" ) +
  xlab( "Week" ) +
  
  # text sizes
  theme( axis.title.x = element_text( size = 14),
         axis.title.y = element_text( size = 14 ),
         axis.text.x = element_text( size = 12),
         axis.text.y = element_text( size = 12 ),
         legend.title = element_text( size = 14 ),
         legend.text = element_text( size = 12 ),) +
  
  # add date labels at tails of the lines
  geom_label_repel( label = label.data$line.label,
                    nudge_x = label.data$nudge.x,
                    nudge_y = label.data$nudge.y,
                    segment.color = "grey50",
                    colour = "black",
                    family = "Avenir",
                    size = 3.4,
                    segment.size = 0.3 ) +
  scale_x_continuous( breaks = brks ) +
  
  # vertical different segment metro lines
  geom_linerange(x = x.par, 
                 ymin = hi.pt.met, 
                 ymax = lo.pt.met,
                 size = 0.01,
                 show.legend = F) +
  geom_segment( # upper left-pointing segment of vertical range bar
    x = x.par, xend = x.par -1,
    y = hi.pt.met, yend = hi.pt.met,
    size = 0.01,
    show.legend = F ) +
  geom_segment( # lower left-pointing segment of vertical range bar
    x = x.par, xend = x.par -1,
    y = lo.pt.met, yend = lo.pt.met,
    size = 0.011,
    show.legend = F ) +
  
  # vertical different segment non-metro lines
  geom_linerange(x = x.par, 
                 ymin = hi.pt.nonmet, 
                 ymax = lo.pt.nonmet,
                 size = 0.01,
                 show.legend = F) +
  geom_segment( # upper left-pointing segment of vertical range bar
    x = x.par, xend = x.par -1,
    y = hi.pt.nonmet, yend = hi.pt.nonmet,
    size = 0.011,
    show.legend = F ) +
  geom_segment( # lower left-pointing segment of vertical range bar
    x = x.par, xend = x.par -1,
    y = lo.pt.nonmet, yend = lo.pt.nonmet,
    size = 0.011,
    show.legend = F ) + 
  
  # text label of differences in deaths in metro counties
  annotate( "text", label = unname( TeX( paste0(  "$\\overset{\\Delta = ",
                                                  hi.pt.met-lo.pt.met, "}",
                                                  "{Deaths}$" ) ) ), 
            x = x.par + x.nudge, 
            y = met.txt.pos.y,
            size = 4,
            family = "Avenir") +
  
  # text label of differences in deaths in metro counties
  annotate( "text", label = unname( TeX( paste0(  "$\\overset{\\Delta = ",
                                                  hi.pt.nonmet-lo.pt.nonmet, "}",
                                                  "{Deaths}$" ) ) ), 
            x = x.par + x.nudge, 
            y = nonmet.txt.pos.y,
            size = 4,
            family = "Avenir") +
  coord_cartesian(
    xlim = c( 0, max( label.data$week ) + 28 ) ) + # controls how much text gets into (does not bleed off of) plot
  theme( text = element_text( family = "Avenir" ))


## (1.5) Save Plot  ##

# this plot will be combined with the chloropleths generated in the descriptives
# see: "R/03-descriptives.R"
ggsave( width = 7.79, height = 4.68,
        "04-Tables-Figures/00-descriptives/time-series-stratified.tiff" )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------




