###------------------------------------------------------------
###   03-DESCRIPTIVES
###------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we will generate table 1 (epidemiological characteristics of the 
# study samples--which are counties in this analysis), a correlation matrix
# and some maps. We will use the observations employed in the analysis file
# ("03-Data-Rodeo/01-analytic-data.rds").
#
# INPUT DATA FILES: 
# "03-Data-Rodeo/01-analytic-data.rds"
#
# OUPUT FILES:
# i. "01-cor-matrix-b.txt"
# ii. "02-cor-matrix-sub-b.txt"
# iii. "03-table-1.txt"
# iv. "04-new-bivariate-map-continent-all.png"
#
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

library( radiant.data ) # for weighted means
library( tidyverse )
library( sf )           # map-making
library( tigris )
library( ggpubr )       # for arranging
library( biscale )      # for bivariate maps
library( cowplot )
library( sjstats )      # for weighted wilcoxon rank sum test
library( latex2exp )    # for latex
library( spdep )        # for neighborhood matrices
library( survey )
library( ggspatial )   # for nothern pointing arrow and map scale

# read in helper functions
source( "R/utils.R")


### (1.0) Data Import ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# read in data with observations inlcuded in analysis
d.desc <- readRDS( "03-Data-Rodeo/01-analytic-data.rds" ) %>%
  mutate( jhu.age.adj.mort.rate = ( deaths.jhu.adj / pop ) * 100000 )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (2.0) Correlation Matrix ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (2.1) Columns to Use for Matrix ##

these.columns <- c( "pop", "pop.density", "fi.perc.20", "jhu.age.adj.mort.rate", "age.adj.mort.rate", "median.age","perc.female", 
                    "perc.native","perc.hisp","perc.black","perc.asian",
                     "perc.nh.white", "perc.fb", "pct.emp.trade","pct.emp.trans","ed.1less.than.hspct",
                     "ed.5college.plus.pct","poverty.rate", "disability", "no.health.insur", "no.vehic", "med.hhinc", "avg.hhsize", "gini.index", "elec.2020.margin", 
                    "inc.prop.jhu","perc.vaccinated","ratio.pop.edp","health.index" )

## ---o--- ##


## (2.2) Generate Weighted Pearson Correlation Matrix ##

cor.fct <- cov.wt( d.desc[, these.columns ], wt = d.desc$pop, cor = TRUE )

# round entries
cor.mat <- cor.fct$cor %>% round( digits = 2)

# subset of columns
# keep only columns for FI and mort. rate
keep <- which( colnames( cor.mat ) %in% c( "fi.perc.20", "jhu.age.adj.mort.rate", "age.adj.mort.rate" ) )

cor.mat.sub <- cor.mat[, keep ] # matrix with obnly the three columns listed above in line 66

## ---o--- ##


## (2.3) Clean up the Table Aesthetically ##

# keep lower triangular entries only
cor.mat[ upper.tri( cor.mat ) ] <- ""

# replace entries in the diagonal
diag( cor.mat ) <- "--"

# clean up significant digits
for( i in 1:ncol( cor.mat ) ){
  cor.mat[ , i ] <- str_replace( cor.mat[ , i ], "(\\.\\d)$", "\\10" )
  
  cor.mat[ , i ] <- str_replace( cor.mat[ , i ], "^0$", "0.00" )
}

for( i in 1:ncol( cor.mat.sub ) ){
  cor.mat.sub[ , i ] <- str_replace( cor.mat.sub[ , i ], "(\\.\\d)$", "\\10" )
  
  cor.mat.sub[ , i ] <- str_replace( cor.mat.sub[ , i ], "^0$", "0.00" )
 
  cor.mat.sub[ , i ] <- str_replace( cor.mat.sub[ , i ], "^1$", "--" )

}

# row and column names
nice.names <- c( "Population", "Population Density", "Food Insecurity (%)", "JHU COVID-19 Mortality Rate (Deaths/100k)", "CDC COVID-19 Mortality Rate (Deaths/100k)", "Median Age","% Female", "% Native American", "% Hispanic","% Black","% Asian",
                 "% Non-Hispanic White", "% Foreign Born", "% Employed in Trade/Retail","% Employed in Transportation","% w/ < HS Diploma",
                 "% at least 4-yr college degree","Poverty Rate (%)","% Disabled", "% w/ Health Insurance", "% w/ No Vehicle Access",
                 "Median Household (HH) Income", "Average HH Size", "Gini Index", 
                "2020 General Election Vote Differential", "COVID-19 Incidence (CDC)","% Vaccinated",
                 "Population:Emergency Physicians","Health Index")

rownames( cor.mat ) <- nice.names 
rownames( cor.mat.sub ) <- nice.names 
colnames( cor.mat ) <- nice.names 
colnames( cor.mat.sub ) <- nice.names[ keep ]

## ---o--- ##


## (2.4) Save Tables ##
write.table( cor.mat, "04-Tables-Figures/00-descriptives/01-cor-matrix-b.txt", sep = "," )
write.table( cor.mat.sub, "04-Tables-Figures/00-descriptives/02-cor-matrix-sub-b.txt", sep = "," )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (3.0) Table 1  (Reported as Table 2 in the Manuscript) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (3.1) Loop through Same Variables in Corr Matrix to Make Table 1 ##

# initialize lists to hold results
tabsg <- list( ) # store column for combined sample
tabsg.met <- list( ) # store column for metropolitan counties
tabsg.nonmet <- list( ) # store column for non-metropolitan counties


## Loop through variables to generate table

for ( i in 1:length( these.columns ) ){
  
  # we will not weight population and population density variables by population size
  if ( i %in% c( 1, 2, 4 ) ){
    
    # combined sample
    tabsg[[i]] <- tab1.var.mean( var.name = these.columns[i], data = d.desc, weight = NULL, table.var.name = these.columns[i] )
    
    # metro
    tabsg.met[[i]] <- tab1.var.mean( var.name = these.columns[i], data = d.desc,
                                  weight = NULL, table.var.name = these.columns[i],
                                  strata.var = 'urb.cat.code', strata.level = 'Metropolitan' )
    # non-metro
    tabsg.nonmet[[i]] <- tab1.var.mean( var.name = these.columns[i], data = d.desc,
                                     weight = NULL, table.var.name = these.columns[i],
                                     strata.var = 'urb.cat.code', strata.level = 'Non-metropolitan' )  
  }
  
  # the rest of the variables will be weighted by population size
  if ( i %notin% c( 1, 2, 4 ) ){
    
    # combined sample
    tabsg[[i]] <- tab1.var.mean( var.name = these.columns[i], data = d.desc, weight = 'pop', table.var.name = these.columns[i] )
    
    # metro
    tabsg.met[[i]] <- tab1.var.mean( var.name = these.columns[i], data = d.desc,
                                     weight = 'pop', table.var.name = these.columns[i],
                                     strata.var = 'urb.cat.code', strata.level = 'Metropolitan' )
    # non-metro
    tabsg.nonmet[[i]] <- tab1.var.mean( var.name = these.columns[i], data = d.desc,
                                        weight = 'pop', table.var.name = these.columns[i],
                                        strata.var = 'urb.cat.code', strata.level = 'Non-metropolitan' )  
  }
}

## ---o--- ##


## (3.2) Assemble Table and Make Nice ##

tab.1 <- bind_cols( list_it( tabsg ),
                    list_it( tabsg.met ),
                    list_it( tabsg.nonmet ) )[ , c( 1, 2, 4, 6 ) ]


# column names and sample sizes

comb.n <- nrow( d.desc )
met.n <- nrow( d.desc[ d.desc$urb.cat.code == "Metropolitan", ] )
nonmet.n <- nrow( d.desc[ d.desc$urb.cat.code == "Non-metropolitan", ] )

colnames( tab.1 ) <- c( "Characteristic",
                        paste0( "Combined Sample (n = ", comb.n ),
                        paste0( "Metropolitan (n = ", met.n ),
                        paste0( "Non-Metropolitan (n = ", nonmet.n ) )


# make rownames presentable

tab.1$Characteristic <- nice.names

## ---o--- ##


## (3.3) Wilcoxon Rank-Sum test P-values ##

# store Wilcoxon Rank Sum test p values in vector initialized below
p.out <- vector()

# loop through variables
for( i in 1:length( these.columns ) ){
  
  # population, population density, and mortality rate not weighted
  if ( i %in% c( 1, 2, 4 ) ){
    
    # formula
    mw.form <- formula( paste0( these.columns[i], "~ urb.cat.code"))
    
    p.out[i] <- wilcox.test( mw.form, data = d.desc )$p.value
  }
  
  # remainder of variables are weighted by population size
  if (i %notin% c( 1, 2, 4 ) ){
    
    # weighted WRS test uses the function from the sjstats package
    p.out[i] <- eval( parse( text = paste0( "sjstats::weighted_mannwhitney( x = ", these.columns[i], " , grp = urb.cat.code, 
                      weights = pop, data = d.desc )" ) ) )$p.value
  
  }
}



# clean up vector
p.out <- ifelse( p.out < 0.01, "< 0.01**",
                 ifelse( p.out >= 0.01 & p.out < 0.05, paste0( round( p.out, 2 ), "*" ), round( p.out, 2 ) ) ) 

# append to table 1
tab.1$p <- p.out

# clean up significant digits
for ( i in c( 2:4 ) ){
  
  tab.1[ , i ] <- str_replace( tab.1[ , i ], "(\\.\\d)\\s\\(", "\\10 (")
  tab.1[ , i ] <- str_replace( tab.1[ , i ], "(\\.\\d)\\)", "\\10)" )
  
}

## ---o--- ##


## (3.4) Save Table ##

write.table( tab.1, "04-Tables-Figures//00-descriptives/03-table-1.txt", sep = "," )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (4.0) USA County Chloropleths ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (4.0) Bring in Shapefiles ##

# pull county-level shapefiles from `tigris`
t <- counties( cb = TRUE ) # USA Census tract shapefiles


# extract state FIPS codes by extracting first two characters in GEOID string
state.codes <- str_extract( t$GEOID, "^.{2}")

# remove FIPS codes for Mariana Islands, American Samoa, Guam, Virgin Islands
# Alaska, PR,  and HI (`tidycensus` won't return data for these)
state.codes <- state.codes[ !state.codes %in% c( "69", "60", "66", 
                                                 "78", "15", "02",
                                                 "72") ]

## ---o--- ##


## (4.1) Food Insecurity Chloropleth ##

# breaks in continuous scale (for the legend)
qs.fi <- quantile( d.desc$fi.perc.20, (0:9)/9 )
qs.fi <- setNames( qs.fi, c( round( min(d.desc$fi.perc.20), 2 ), rep( "", 8 ), round( max(d.desc$fi.perc.20), 2 ) ) )

# plot
 continent.fi <- ( 
    
    # wrangle/join shapefile data before plotting
    t %>%
    mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
    left_join( ., d.desc %>% select( -state.code ), by = c( "GEOID" = "fips") ) %>%
    dplyr::filter( state.code %in% state.codes ) %>% # filter only mainland states
    
    # plot 
    ggplot(  ) +
    geom_sf(           
      mapping = aes( fill = fi.perc.20),
      color = NA, size = 0.5 ) +
    coord_sf( crs = st_crs( 2163 ), expand = T ) +
      # add north-poiting compass rose
      ggspatial::annotation_north_arrow( location = "tl", which_north = "true", 
                                         pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                                         style = north_arrow_fancy_orienteering, width = unit(0.8, "cm"), 
                                         height = unit(1, "cm") ) +
      annotation_scale( height = unit(0.14, "cm"), style = "bar" ) + # adds map scale bar
      scale_fill_fermenter( palette = 7, direction = 1, guide = "colourbar", # `scale_fill_fermenter` for discretizing a continuous scale
                            breaks = qs.cdc, na.value = "black" ) +
   scale_fill_fermenter( palette = 13, direction = 1, guide = "colourbar", # `scale_fill_fermenter` for discretizing a continuous scale
                         breaks = qs.fi, na.value = "black" ) +
  theme_minimal( ) +
    theme( text = element_text( family = "Avenir" ),
           legend.title = element_blank(),
           axis.text = element_text( size = 12 ),
           plot.title = element_text( size = 19, hjust = -0.5 ),
           axis.title = element_text( size = 17 ),
           legend.text = element_text( size = 14 ),
           plot.background = element_rect(color = "white")) +
  ggtitle( "2020 County Food Insecurity (%)")  +
  xlab( "Longitude" ) +
  ylab( "Latitude" ) +
   guides( fill = guide_colorbar( ticks.colour = NA,
                                  frame.colour =  "black",
                                  barwidth = 0.6,
                                  barheight = 9 ) )
) %>% 
  ggarrange(., nrow = 1, ncol = 1 ) + # need to pass to `ggarrange` so that label is positioned properly
  # add data source label
  cowplot::draw_plot_label( label = "Data Source: Feeding America",
                              x = 1, y = 1, hjust = 1.04, 
                              vjust = 22.1,
                              family = "Avenir", color = "grey44",
                              size = 12 ) 
 
 ## ---o--- ##
 
 
## (4.2) COVID-19 Mortality (CDC) Chloropleth ##

 # breaks in continuous scale (for the legend)
 qs.cdc <- quantile( d.desc$age.adj.mort.rate, (0:9)/9 )
 qs.cdc <- setNames( qs.cdc, c( round( min(d.desc$age.adj.mort.rate), 2 ), rep( "", 8 ), round( max(d.desc$age.adj.mort.rate), 2 ) ) )
 
 # plot
 continent.cdc <- ( 
   
   # wrangle/join shapefile data before plotting
   t %>%
     mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
     left_join( ., d.desc %>% select( -state.code ), by = c( "GEOID" = "fips") ) %>%
     dplyr::filter( state.code %in% state.codes ) %>% # filter only mainland states
     
     # plot 
     ggplot(  ) +
     geom_sf(           
       mapping = aes( fill = age.adj.mort.rate ),
       color = NA, size = 0.5 ) +
     coord_sf( crs = st_crs( 2163 ), expand = T ) +
     # add north-poiting compass rose
     ggspatial::annotation_north_arrow( location = "tl", which_north = "true", 
                                        pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                                        style = north_arrow_fancy_orienteering, width = unit(0.8, "cm"), 
                                        height = unit(1, "cm") ) +
     scale_fill_fermenter( palette = 7, direction = 1, guide = "colourbar", # `scale_fill_fermenter` for discretizing a continuous scale
                           breaks = qs.cdc, na.value = "black" ) +
     theme_minimal( ) +
     theme( text = element_text( family = "Avenir" ),
            legend.title = element_blank(),
            axis.text = element_text( size = 12 ),
            plot.title = element_text( size = 19, hjust = -1.8 ),
            axis.title = element_text( size = 17 ),
            legend.text = element_text( size = 14 ),
            plot.background = element_rect(color = "white")) +
     ggtitle( "COVID-19 Mortality Rate (Deaths/100k)")  +
     xlab( "" ) +
     ylab( "" ) +
     guides( fill = guide_colorbar( ticks.colour = NA,
                                    frame.colour =  "black",
                                    barwidth = 0.6,
                                    barheight = 9 ) )
 ) %>% 
   ggarrange(., nrow = 1, ncol = 1 ) + # need to pass to `ggarrange` so that label is postitioned properl
   # add data source label
   cowplot::draw_plot_label( label = "Data Source: CDC",
                             x = 1, y = 1, hjust = 1.2, 
                             vjust = 26.5,
                             family = "Avenir", color = "grey44",
                             size = 12 )  
 
 ## ---o--- ##
  

## (4.3) COVID-19 Mortality (JHU) Chloropleth ##

# breaks in continuous scale (for the legend)
qs.jhu <- quantile( d.desc$jhu.age.adj.mort.rate, (0:9)/9 )
qs.jhu <- setNames( qs.jhu, c( round( min(d.desc$jhu.age.adj.mort.rate), 2 ), rep( "", 8 ), round( max(d.desc$jhu.age.adj.mort.rate), 2 ) ) )

 continent.jhu <- ( 
      
      # wrangle/join shapefile data before plotting
      t %>%
      mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
      left_join( ., d.desc %>% select( -state.code ), by = c( "GEOID" = "fips") ) %>%
      dplyr::filter( state.code %in% state.codes ) %>% # filter only mainland states
      
      # plot 
      ggplot(  ) +
      geom_sf(           
        mapping = aes( fill = jhu.age.adj.mort.rate ),
        color = NA, size = 0.5 ) +
      coord_sf( crs = st_crs( 2163 ), expand = T ) +
        # add north-poiting compass rose
        ggspatial::annotation_north_arrow( location = "tl", which_north = "true", 
                                           pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                                           style = north_arrow_fancy_orienteering, width = unit(0.8, "cm"), 
                                           height = unit(1, "cm") ) +
        annotation_scale( height = unit(0.14, "cm"), style = "bar" ) + # adds map scale bar
        scale_fill_fermenter( palette = 7, direction = 1, guide = "colourbar", # `scale_fill_fermenter` for discretizing a continuous scale
                              breaks = qs.cdc, na.value = "black" ) +
      scale_fill_fermenter( palette = 7, direction = 1, guide = "colourbar", # `scale_fill_fermenter` for discretizing a continuous scale
                            breaks = qs.jhu, na.value = "black" ) +
      theme_minimal( ) +
        theme( text = element_text( family = "Avenir" ),
               legend.title = element_blank(),
               axis.text = element_text( size = 12 ),
               plot.title = element_text( size = 19, hjust = -1.8 ),
               axis.title = element_text( size = 17 ),
               legend.text = element_text( size = 14 ),
               plot.background = element_rect(color = "white")) +  #unit(c(top, right, bottom, left), units))  
      ggtitle( "COVID-19 Mortality Rate (Deaths/100k)")  +
      xlab( "" ) +
      ylab( "Latitude" ) +
      guides( fill = guide_colorbar( ticks.colour = NA,
                                     frame.colour =  "black",
                                     barwidth = 0.6,
                                     barheight = 9 ) )
) %>% 
    ggarrange(., nrow = 1, ncol = 1 ) + # need to pass to `ggarrange` so that label is positioned properly
    # add data source label
    cowplot::draw_plot_label( label = "Data Source: JHU",
                              x = 1, y = 1, hjust = 1.2, 
                              vjust = 26.5,
                              family = "Avenir", color = "grey44",
                              size = 12 )  

 ## ---o--- ##
 
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (5.0) Bivariate Map ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (5.1) Prepare Data ##
 
# use `biscale` package to make bivariate map variable (`bi_class`)
# NOTE: we make the bivariate map with the JHU COVID-19 mortality rates as the Y variable and not the CDC
 
mainland.biv <-
  
  # wrangle/join shapefile data before plotting
  t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  left_join( ., d.desc %>% select( -state.code ), by = c( "GEOID" = "fips") ) %>%
  bi_class(., x = fi.perc.20, y = jhu.age.adj.mort.rate, style = "quantile", dim = 3) %>%
  filter( state.code %in% state.codes )  # filter only continental states


# create pallette for color scale
# source: https://observablehq.com/@d3/bivariate-choropleth 
custom_pal3 <- c(
  "1-1" = "#e8e8e8", # low x, low y
  "2-1" = "#ace4e4",
  "3-1" = "#5ac8c8", # high x, low y
  "1-2" = "#dfb0d6",
  "2-2" = "#a5add3", # medium x, medium y
  "3-2" = "#5698b9",
  "1-3" = "#be64ac", # low x, high y
  "2-3" = "#8c62aa",
  "3-3" = "#3b4994" # high x, high y
)

## ---o--- ##


## (5.2) Plot Continental U.S. ##

biv.plot <- mainland.biv %>%
ggplot(  ) +
  geom_sf(           
    mapping = aes( fill = bi_class ),
    color = NA, size = 0.5 ) +
  coord_sf( crs = st_crs( 2163 ), expand = T ) +
  # add north-poiting compass rose
  ggspatial::annotation_north_arrow( location = "tl", which_north = "true", 
                                     pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                                     style = north_arrow_fancy_orienteering, width = unit(0.8, "cm"), 
                                     height = unit(1, "cm") ) +
scale_fill_manual( values = custom_pal3, na.value = "black" ) +
  theme_minimal( ) +
  theme( text = element_text( family = "Avenir" ),
         legend.title = element_blank(),
         axis.text = element_text( size = 12 ),
         plot.title = element_text( size = 19, hjust= -0.1 ),
         axis.title = element_text( size = 17 ),
         legend.text = element_text( size = 14 ),
         plot.background = element_rect(color = "white"),
         legend.position = "none" ) +
  ggtitle( TeX( "Bivariate Map$^a$" ) ) +
  xlab( "Longitude" ) +
  ylab( "" ) 

## ---o--- ##


## (5.3) Generate Legend Grob ##

legend <- bi_legend(  pal = custom_pal3,
                      dim = 3,
                      ylab = "Food Insecurity (%)",
                      xlab = "COVID-19 Deaths/100k (JHU)",
                      size = 14 ) + 
  ylab("Food Insecurity (%)\n(tertile)")+
  xlab("COVID-19\nDeaths/100k\n(tertile)")+
  theme(plot.margin=unit(c( 0,0,0.1,0), "cm"),
                                         text = element_text( family = "Avenir", face = "bold" ),
                                         axis.text.y = element_text( angle = 0, size = 10,color = "gray44" ),
                                         axis.text.x = element_text( size = 10.5, color = "gray44", margin = unit(c(0.2, 0, 0, 0), "mm") ),
                                         axis.title.x = element_text( size = 12,
                                                                      margin = unit(c(3, 0, 0, 0), "mm")),
                                         axis.title.y = element_text(angle = 90, vjust = 0.4,
                                                                     size = 12,
                                                                     margin = unit(c(0, 3, 0, 0), "mm")),
        plot.background = element_rect(fill='transparent')) # make legend background transparent

## ---o--- ##


### (5.4) Create Missing Data Legend Grob ###

leg.p <- generate_NA_legend( "Missing Data", "black")


## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (6.0) Assemble Chloropleths and Bivariate Map ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# arrange and set margins accordingly so that all maps fit into final image
this.p <-  ggarrange( 
  
  ggarrange( continent.jhu + theme( plot.margin = unit(c(-10,.5,-15,-0.35), "lines")),
             continent.fi + theme( plot.margin = unit(c(-10,1,-10,-0.35), "lines")), nrow = 2, ncol  = 1 ),
  
  ggarrange( continent.cdc+theme(plot.margin = unit(c(-10,0,-15,-0.35), "lines" ) ), 
             biv.plot+theme(plot.margin = unit(c(-10,8,-10,-1.8), "lines")), nrow = 2, ncol  = 1),
  nrow = 1, ncol = 2
  
  )


# add legend layer 
this.p.leg <- ggdraw( this.p ) + # increase bottom margin otherwise the map bleeds off the page
  draw_plot( legend, 0.80, 0.06, 0.27, 0.27, scale = 0.63 ) +
  draw_plot( leg.p, 0.914, 0.428, 0.05, 0.05, scale = 0.05 )

# save
ggsave( "04-Tables-Figures/00-descriptives/04-new-bivariate-map-continent-all.tiff" ,
        units = "px", width = (4124), height = (2808) , plot = this.p.leg )  

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (7.0) Mean and Total Mortality Counts and Rates ### 
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# mean mortality rate (age=adjusted--jhu)
mean( d.desc$jhu.age.adj.mort.rate )

# mean mortality rate (age=adjusted--cdc)
mean( d.desc$age.adj.mort.rate )

# total crude COVID-19 deaths (jhu)
sum( d.desc$deaths.jhu )

# total crude COVID-19 deaths (cdc)
sum( d.desc$raw.deaths.cdc )

# no. counties with zero deaths
sum( d.desc$deaths.jhu == 0 )

# highest count totals by state
d.desc %>%
  group_by( state ) %>%
  mutate( state.deaths = sum( deaths.jhu ) ) %>%
  ungroup() %>%
  arrange( desc( state.deaths ) ) %>%
  select( state, state.deaths ) %>%
  distinct() 

# highest average mortality rates by state
d.desc %>%
  group_by( state ) %>%
  mutate( state.deaths = mean( jhu.age.adj.mort.rate ) ) %>%
  ungroup() %>%
  arrange( desc( state.deaths ) ) %>%
  select( state, state.deaths ) %>%
  distinct() %>% data.frame()

# highest average mortality counts by county
d.desc %>%
  group_by( fips ) %>%
  mutate( county.deaths = sum( deaths.jhu ) ) %>%
  ungroup() %>%
  arrange( desc( county.deaths ) ) %>%
  select( state, county, county.deaths ) %>%
  distinct() %>% data.frame()

# highest average mortality rates by county
d.desc %>%
  group_by( fips ) %>%
  mutate( county.deaths = mean( jhu.age.adj.mort.rate ) ) %>%
  ungroup() %>%
  arrange( desc( county.deaths ) ) %>%
  select( state, county, county.deaths ) %>%
  distinct() %>% data.frame()

# highest average FI rates by state
d.desc %>%
  group_by( state ) %>%
  mutate( state.fi = mean( fi.perc.20 ) ) %>%
  ungroup() %>%
  arrange( desc( state.fi ) ) %>%
  select( state, state.fi ) %>%
  distinct() %>% data.frame()

# highest average FI rates by county
d.desc %>%
  group_by( fips ) %>%
  mutate( county.fi = mean( fi.perc.20) ) %>%
  ungroup() %>%
  arrange( desc( county.fi ) ) %>%
  select( state, county, county.fi ) %>%
  distinct() %>% data.frame()


# total COVID-19 deaths in the time period (3/25/2020-12/25/2021)
sum( d.desc$deaths.jhu, na.rm = T ) # 773154 in JHU dataset
sum( d.desc$deaths.cdc, na.rm = T ) # 495350 in CDC dataset

# differences in deaths across datasets and across urban-rural status
met.deaths.jhu <- sum( d.desc[ d.desc$urb.cat.code=="Metropolitan",]$deaths.jhu, na.rm = T ) # 773154 in JHU dataset
nonmet.deaths.jhu <- sum( d.desc[ d.desc$urb.cat.code=="Non-metropolitan",]$deaths.jhu, na.rm = T ) # 495350 in CDC dataset

met.deaths.cdc <- sum( d.desc[ d.desc$urb.cat.code=="Metropolitan",]$deaths.cdc, na.rm = T ) # 773154 in cdc dataset
nonmet.deaths.cdc <- sum( d.desc[ d.desc$urb.cat.code=="Non-metropolitan",]$deaths.cdc, na.rm = T ) # 495350 in CDC dataset

( met.deaths.jhu - met.deaths.cdc ) / met.deaths.jhu
( nonmet.deaths.jhu - nonmet.deaths.cdc ) / nonmet.deaths.jhu

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (8.0) Moran's I ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# example code is found: https://mgimond.github.io/simple_moransI_example/

d.poly <- t %>%
  left_join( d.desc, ., by = c( "fips" = "GEOID") ) %>%
  dplyr::filter( state.code %in% state.codes ) %>% # filter only mainland states
  st_as_sf()
  
## (8.1) Generate Neighborhood Matrix ##
( nb <- poly2nb( d.poly ) )

# looks like we have an empty neighbor set at row 417 (San Juan islands in Washington state) so we will remove them for this analysis only
nb.min <- poly2nb( d.poly[-2941,] )

# the x and y vectors minus that of San Juan
x <- d.poly[-2941,]$fi.perc.20
y <- d.poly[-2941,]$jhu.age.adj.mort.rate

# assign equal weights for computing neighboring variables
lw <- nb2listw( nb.min, style = "W", zero.policy = TRUE )

## ---o--- ##


## (8.2) Global Bivariate Moran's I ##

# (see: https://geodacenter.github.io/workbook/5b_global_adv/lab5b.html)
moran_bv( x, y, lw, nsim = 500 )

## ---o--- ##


## (8.3) Univariate Moran's I for the Dependent and Independent Variables ##

# obtain statistics using `spdep::moran()`
moran( x, lw, n = length(x), Szero( lw ) )
moran( y, lw, n = length(y), Szero( lw ) )

# we can also plot and obtain Moran's I statistic by obtaining the slope of the lags regressed on the variables
space.lag.x <- lag.listw(lw, x )
space.lag.y <- lag.listw(lw, y )

dev.copy( png,'04-Tables-Figures/08-other-supplementary-files/02-moran-fi-plot.png' )
plot( space.lag.x ~ d.poly[-2941,]$fi.perc.20, pch = 16, asp = 1, 
      xlab = "Food Insecurity",
      ylab = "Spatial Lag" )
m.x <- lm( space.lag.x ~ d.poly[-2941,]$fi.perc.20 )
abline(m.x, col="blue")
coef( m.x )[2]
text( x = 5, y = 20, label = unname( TeX( paste0( "$\\beta$ = ", round( coef( m.x )[2], 2) ) ) ) ) 
dev.off()

dev.copy( png,'04-Tables-Figures/08-other-supplementary-files/03a-moran-covid-plot-jhu.png' )
plot( space.lag.y ~ d.poly[-2941,]$jhu.age.adj.mort.rate, pch = 16, asp = 1,
      xlab = "COVID-19 Mortality Rate (JHU)",
      ylab = "Spatial Lag" )
m.y <- lm( space.lag.y ~ d.poly[-2941,]$jhu.age.adj.mort.rate )
abline(m.y, col="blue")
coef( m.y )[2]
text( x = 200, y = -200, label = unname( TeX( paste0( "$\\beta$ = ", round( coef( m.y )[2], 2) ) ) ) ) 
dev.off()

dev.copy( png,'04-Tables-Figures/08-other-supplementary-files/03b-moran-covid-plot-cdc.png' )
plot( space.lag.y ~ d.poly[-2941,]$age.adj.mort.rate, pch = 16, asp = 1,
      xlab = "COVID-19 Mortality Rate (CDC)",
      ylab = "Spatial Lag" )
m.y <- lm( space.lag.y ~ d.poly[-2941,]$age.adj.mort.rate )
abline(m.y, col="blue")
coef( m.y )[2]
text( x = 1500, y = 1000, label = unname( TeX( paste0( "$\\beta$ = ", round( coef( m.y )[2], 2) ) ) ) ) 
dev.off()

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (9.0) Correlation of No. of Neighbors to Other Covariates ### 
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# bring in neighborhood object from above in section 8.1
nnbs <- card( nb ) 

d.nb <- d.desc %>%
  mutate( nnbs = nnbs )

cor( d.nb$nnbs, d.nb$elec.2020.margin )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



