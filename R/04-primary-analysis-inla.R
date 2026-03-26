###-------------------------------------------------------------
###   04-FIT CAR MODELS WITH INLA 
###------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we fit the conditional autoregressive models in INLA using two functions written for this purpose.
# The `car_inla_analysis` function carries out the entire analysis.
# The models we fit include two random effects: a spatially structured random effect at the county level and
# an unstructured random effect at the state level. NOTE: the `car_inla_analysis` fct stipulates the default log-
# gamma distributions for the hyperparameters of the models we fit (we later carry out a sensitivity analysis)
# with different priors specified). All internal functions are found in "utils.R".
#
# The `car_inla_analysis` fct also returns a number of maps, plots, and metadata that are part of the 
# results and are otherwise useful. A number of those maps/plots also rely on internal functions housed in 
# "utils.R". i.e., `inla_var_explained`-- for computing the proportion of variance explained by the random
# effects, `plot_rr` -- for plotting chloropleths/maps, `f_plot_inla` --  for generating a forest plot,
# `inla_results_save` --  for saving the results in a local directory.
# 
# 
#
# INPUT DATA FILE: "03-Data-Rodeo/analytic-data.rds"
#
#
#
# Resources (Accessed 3 August 2023): 
# i. Pub explaining and using the map decomposition technique: https://doi.org/10.1111/j.1538-4632.2004.tb01132.x
# ii. Pub this analysis is modeled after: https://doi.org/10.1093/aje/kwy026
# iii. Pub with more detailed description of CAR model, its implementation in INLA, and examples: https://doi.org/10.1016/j.sste.2013.07.003
# iv. the R-INLA Project: https://www.r-inla.org/
#
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

library( INLA )
library( tidyverse )
library( spdep ) # for creating neighborhood matrix
library( ggpubr ) # for arranging plots
library( tigris ) # shapefiles
library( latex2exp ) # latex in plots
library( ggpattern )
library( viridis ) # color palletes for ggplot
library( reshape2 ) # for melting data (needed in `car_inla_analysis` fct)
library( scales ) # for number formatting in ggplot (i.e., numbers past the decimal mark)
library( cowplot ) # for figure labels
library( ggspatial )

# helper functions
source( "R/utils.R" )

### (1.0) Data Import and Preparation ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# read-in data and omit states not included in analysis
d <- readRDS( "03-Data-Rodeo/01-analytic-data.rds" ) 

# pull county-level shapefiles from `tigris`
t <- counties( cb = TRUE ) # USA Census tract shapefiles

# merge geometries to our dataset
d.poly <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  dplyr::select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d, by = c( "fips" ) ) %>%
  filter( fips %in% d$fips ) %>%
  mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
          E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using CDC data)
          exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
          exp.cases.cdc = ceiling( mean( .$inc.prop.cdc, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using CDC data)
          sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
          sir.cdc = cases.cdc / exp.cases.cdc, # SIR computation, county-level (CDC data)
          exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
          exp.cases.cdc.state = ceiling( mean( .$inc.prop.cdc.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using CDC data)
          sir.jhu.state = cases.jhu.state / exp.cases.jhu.state, # SIR computation, state-level (JHU data)
          sir.cdc.state= cases.cdc.state / exp.cases.cdc.state )# SIR computation, state-level (JHU data)
## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
## for each dataset we feed through `car_inla_analysis`

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (2.0) Formula Specifications ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (2.1) Models with Full Set of Covariates ##

# model specs
source( "R/model-specs.R" )

## ---o--- ##


# ---------------------------------------------------------------------------------------------------------------------------------------------------------


### (3.0) INLA Analysis: National (Nearly All U.S. Counties) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (3.1) Call INLA ##

# states to remove from maps

states.fips.out <- unique( c( "69", "60", "66", "78", "15", "02", "72" ) )

# jhu analysis
res.jhu.1 <- car_inla_analysis( d.jhu = d.poly, d.cdc = NULL,
                                formula.jhu = f.jhu.model.5,
                                which.model = "jhu", E = "E_d.j",
                                term = "fi", null.model.formula = f.jhu.model.null,
                                model.2.formula = f.jhu.model.2, 
                                model.3.formula = f.jhu.model.3,
                                model.4.formula = f.jhu.model.4,
                                state.codes.omit = states.fips.out )

# plot state-level residuals to identify potential outliers
ags <- aggregate( d$deaths.jhu.adj, FUN = "mean", by = list( d$state) )

res.df <- left_join( res.jhu.1$model.fit$summary.random$state,
                     ags, by = c("ID" = "Group.1" ) )

ggplot( data = res.df ) +
  geom_point( mapping = aes( x = ID, y = mean ) ) +
  ggtitle( "state-level random effects" )

ggsave( "/Users/mainovieytesca/Documents/GitHub/COVID-FI-Mortality/04-Tables-Figures/01-main-analysis/jhu-results/state-level-random-effects.tiff",
        width = 12, height = 4 )


# cdc analysis
res.cdc.1 <- car_inla_analysis( d.cdc = d.poly, 
                                formula.cdc = f.cdc.model.5,
                                which.model = "cdc", E = "E_d.c",
                                term = "fi", null.model.formula = f.cdc.model.null,
                                model.2.formula = f.cdc.model.2, 
                                model.3.formula = f.cdc.model.3,
                                model.4.formula = f.cdc.model.4,
                                state.codes.omit = states.fips.out )

## ---o--- ##


## (3.2) Save Results ##

# dir.create( "04-Tables-Figures/01-main-analysis" )
# dir.create( "04-Tables-Figures/01-main-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/01-main-analysis/cdc-results" )

inla_results_save( res.jhu.1, path = paste0( "04-Tables-Figures/01-main-analysis/jhu-results/" ),
                   tag = "jhu-main", map.decomp.height = 2060 )

inla_results_save( res.cdc.1, path = paste0( "04-Tables-Figures/01-main-analysis/cdc-results/" ),
                   tag = "cdc-main", map.decomp.height = 2060  )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (4.0) INLA Analysis: Northeast Census Region ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (4.1) Subset Data Based on Census Region ##

d.ne <- d %>%
  filter( state %in% c( "ME", "NH", "VT", "MA", "RI", "CT", "NY", "NJ", "PA" ) )

states.out.ne <- unlist( c( unique( d[ d$state %notin%  c( "ME", "NH", "VT", "MA", "RI", "CT", "NY", "NJ", "PA" ), "state.code" ] ),
                            "69", "60", "66", "78", "15", "02", "72" ) )

# wrangle/join shapefile data before plotting
( d.poly.ne <- t %>%
    mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
    select( -state.code ) %>%
    dplyr::rename( fips = GEOID ) %>%
    left_join( ., d.ne, by = c( "fips" ) ) %>%
    filter( fips %in% d.ne$fips ) %>%
    mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop  ), # expected death count for offset term (using JHU data)
            E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using CDC data)
            exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
            exp.cases.cdc = ceiling( mean( .$inc.prop.cdc, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using CDC data)
            sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
            sir.cdc = cases.cdc / exp.cases.cdc, # SIR computation, county-level (CDC data)
            exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
            exp.cases.cdc.state = ceiling( mean( .$inc.prop.cdc.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using CDC data)
            sir.jhu.state = cases.jhu.state / exp.cases.jhu.state, # SIR computation, state-level (JHU data)
            sir.cdc.state= cases.cdc.state / exp.cases.cdc.state ) ) # SIR computation, state-level (JHU data)
  ## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
  ## for each dataset we feed through `car_inla_analysis`
  

## ---o--- ##


## (4.2) Call INLA ##

# jhu
res.jhu.ne <- car_inla_analysis( d.jhu = d.poly.ne, d.cdc = NULL,
                                 formula.jhu = f.jhu.model.5,
                                 which.model = "jhu", E = "E_d.j",
                                 term = "fi", null.model.formula = f.jhu.model.null,
                                 model.2.formula = f.jhu.model.2, 
                                 model.3.formula = f.jhu.model.3,
                                 model.4.formula = f.jhu.model.4,
                                 state.codes.omit = states.out.ne  )

# cdc
res.cdc.ne <- car_inla_analysis( d.jhu = NULL, d.cdc = d.poly.ne,
                                 formula.cdc = f.cdc.model.5,
                                 which.model = "cdc", E = "E_d.c",
                                 term = "fi", null.model.formula = f.cdc.model.null,
                                 model.2.formula = f.cdc.model.2, 
                                 model.3.formula = f.cdc.model.3,
                                 model.4.formula = f.cdc.model.4,
                                 state.codes.omit = states.out.ne )

## ---o--- ##


## (4.3) Save Results ##

# dir.create( "04-Tables-Figures/02-northeast-analysis" )
# dir.create( "04-Tables-Figures/02-northeast-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/02-northeast-analysis/cdc-results" )

inla_results_save( res.jhu.ne, path = paste0( "04-Tables-Figures/02-northeast-analysis/jhu-results/" ),
                   tag = "jhu-ne", map.decomp.height = 2060  )

inla_results_save( res.cdc.ne, path = paste0( "04-Tables-Figures/02-northeast-analysis/cdc-results/" ),
                   tag = "cdc-ne", map.decomp.height = 2060  )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (5.0) INLA Analysis: South Census Region ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (5.1) Subset Data Based on Census Region ##

d.so <- d %>%
  filter( state %in% c( "MD", "DE", "DC", "WV", "VA", "NC", "KY", "TN", 
                        "SC","GA", "FL", "AL", "MS", "AR", "LA", "OK", "TX" ))

states.out.so <- unlist( c( unique( d[ d$state %notin%  c( "MD", "DE", "DC", "WV", "VA", "NC", "KY", "TN", 
                                                           "SC","GA", "FL", "AL", "MS", "AR", "LA", "OK", "TX" ), "state.code" ] ),
                            "69", "60", "66", "78", "15", "02", "72" ) )

# wrangle/join shapefile data before plotting
( d.poly.so <- t %>%
    mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
    select( -state.code ) %>%
    dplyr::rename( fips = GEOID ) %>%
    left_join( ., d.so, by = c( "fips" ) ) %>%
    filter( fips %in% d.so$fips ) %>%
    mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop  ), # expected death count for offset term (using JHU data)
            E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using CDC data)
            exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
            exp.cases.cdc = ceiling( mean( .$inc.prop.cdc, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using CDC data)
            sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
            sir.cdc = cases.cdc / exp.cases.cdc, # SIR computation, county-level (CDC data)
            exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
            exp.cases.cdc.state = ceiling( mean( .$inc.prop.cdc.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using CDC data)
            sir.jhu.state = cases.jhu.state / exp.cases.jhu.state, # SIR computation, state-level (JHU data)
            sir.cdc.state= cases.cdc.state / exp.cases.cdc.state ) ) # SIR computation, state-level (JHU data)
  ## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
  ## for each dataset we feed through `car_inla_analysis`
  

## ---o--- ##


## (5.2) Call INLA ##

# jhu
res.jhu.so <- car_inla_analysis( d.jhu = d.poly.so, d.cdc = NULL,
                                 formula.jhu = f.jhu.model.5,
                                 which.model = "jhu", E = "E_d.j",
                                 term = "fi", null.model.formula = f.jhu.model.null,
                                 model.2.formula = f.jhu.model.2, 
                                 model.3.formula = f.jhu.model.3,
                                 model.4.formula = f.jhu.model.4,
                                 state.codes.omit = states.out.so  )
# cdc
res.cdc.so <- car_inla_analysis( d.jhu = NULL, d.cdc = d.poly.so,
                                 formula.cdc = f.cdc.model.5,
                                 which.model = "cdc", E = "E_d.c",
                                 term = "fi", null.model.formula = f.cdc.model.null,
                                 model.2.formula = f.cdc.model.2, 
                                 model.3.formula = f.cdc.model.3,
                                 model.4.formula = f.cdc.model.4,
                                 state.codes.omit = states.out.so )

## ---o--- ##


## (5.3) Save Results ##

# dir.create( "04-Tables-Figures/03-south-analysis" )
# dir.create( "04-Tables-Figures/03-south-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/03-south-analysis/cdc-results" )

inla_results_save( res.jhu.so, path = paste0( "04-Tables-Figures/03-south-analysis/jhu-results/" ),
                   tag = "jhu-so", map.decomp.height = 2060  )

inla_results_save( res.cdc.so, path = paste0( "04-Tables-Figures/03-south-analysis/cdc-results/" ),
                   tag = "cdc-so", map.decomp.height = 2060  )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (6.0) INLA Analysis: Midwest Census Region ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (6.1) Subset Data Based on Census Region ##

d.mw <- d %>%
  filter( state %in% c( "OH", "IN", "MI", "IL", "WI", "MN", "IA", "MO","ND", "SD","NE", "KS"))

states.out.mw <- unlist( c( unique( d[ d$state %notin%  c( "OH", "IN", "MI", "IL", "WI", "MN", "IA", "MO","ND", "SD","NE", "KS"), "state.code" ] ),
                            "69", "60", "66", "78", "15", "02", "72" ) )

# wrangle/join shapefile data before plotting
( d.poly.mw <- t %>%
    mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
    select( -state.code ) %>%
    dplyr::rename( fips = GEOID ) %>%
    left_join( ., d.mw, by = c( "fips" ) ) %>%
    filter( fips %in% d.mw$fips ) %>%
    mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
            E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using CDC data)
            exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
            exp.cases.cdc = ceiling( mean( .$inc.prop.cdc, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using CDC data)
            sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
            sir.cdc = cases.cdc / exp.cases.cdc, # SIR computation, county-level (CDC data)
            exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
            exp.cases.cdc.state = ceiling( mean( .$inc.prop.cdc.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using CDC data)
            sir.jhu.state = cases.jhu.state / exp.cases.jhu.state, # SIR computation, state-level (JHU data)
            sir.cdc.state= cases.cdc.state / exp.cases.cdc.state ) ) # SIR computation, state-level (JHU data)
  ## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
  ## for each dataset we feed through `car_inla_analysis`
  
## ---o--- ##


## (6.2) Call INLA ##

# jhu
res.jhu.mw <- car_inla_analysis( d.jhu = d.poly.mw, d.cdc = d.poly.mw,
                                 formula.jhu = f.jhu.model.5,
                                 which.model = "jhu", E = "E_d.j",
                                 term = "fi", null.model.formula = f.jhu.model.null,
                                 model.2.formula = f.jhu.model.2, 
                                 model.3.formula = f.jhu.model.3,
                                 model.4.formula = f.jhu.model.4,
                                 state.codes.omit = states.out.mw )

res.cdc.mw <- car_inla_analysis( d.jhu = NULL, d.cdc = d.poly.mw,
                                 formula.cdc = f.cdc.model.5,
                                 which.model = "cdc", E = "E_d.c",
                                 term = "fi", null.model.formula = f.cdc.model.null,
                                 model.2.formula = f.cdc.model.2, 
                                 model.3.formula = f.cdc.model.3,
                                 model.4.formula = f.cdc.model.4,
                                 state.codes.omit = states.out.mw )

## ---o--- ##


## (6.3) Save Results ##

# dir.create( "04-Tables-Figures/04-midwest-analysis" )
# dir.create( "04-Tables-Figures/04-midwest-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/04-midwest-analysis/cdc-results" )

inla_results_save( res.jhu.mw, path = paste0( "04-Tables-Figures/04-midwest-analysis/jhu-results/" ),
                   tag = "jhu-mw", map.decomp.height = 2060  )

inla_results_save( res.cdc.mw, path = paste0( "04-Tables-Figures/04-midwest-analysis/cdc-results/" ),
                   tag = "cdc-mw", map.decomp.height = 2060  )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (7.0) INLA Analysis: West Census Region ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (7.1) Subset Data Based on Census Region ##

d.w <- d %>%
  filter( state %in% c( "NM", "CO", "WY", "MT", "ID", "UT", "AZ", "NV", "CA", "OR", "WA"))

states.out.w <- unlist( c( unique( d[ d$state %notin%  c( "NM", "CO", "WY", "MT", "ID", "UT", "AZ", "NV", "CA", "OR", "WA"), "state.code" ] ),
                           "69", "60", "66", "78", "15", "02", "72" ) )

# wrangle/join shapefile data before plotting
( d.poly.w <- t %>%
    mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
    select( -state.code ) %>%
    dplyr::rename( fips = GEOID ) %>%
    left_join( ., d.w, by = c( "fips" ) ) %>%
    filter( fips %in% d.w$fips ) %>%
    mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
            E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using CDC data)
            exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
            exp.cases.cdc = ceiling( mean( .$inc.prop.cdc, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using CDC data)
            sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
            sir.cdc = cases.cdc / exp.cases.cdc, # SIR computation, county-level (CDC data)
            exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
            exp.cases.cdc.state = ceiling( mean( .$inc.prop.cdc.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using CDC data)
            sir.jhu.state = cases.jhu.state / exp.cases.jhu.state, # SIR computation, state-level (JHU data)
            sir.cdc.state= cases.cdc.state / exp.cases.cdc.state ) ) # SIR computation, state-level (JHU data)
  ## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
  ## for each dataset we feed through `car_inla_analysis`
  

## ---o--- ##


## (7.2) Call INLA ##

# jhu
res.jhu.w <- car_inla_analysis( d.jhu = d.poly.w, d.cdc = d.poly.w,
                                formula.jhu = f.jhu.model.5,
                                which.model = "jhu", E = "E_d.j",
                                term = "fi", null.model.formula = f.jhu.model.null,
                                model.2.formula = f.jhu.model.2, 
                                model.3.formula = f.jhu.model.3,
                                model.4.formula = f.jhu.model.4,
                                state.codes.omit = states.out.w )
# cdc
res.cdc.w <- car_inla_analysis( d.jhu = NULL, d.cdc = d.poly.w,
                                formula.cdc = f.cdc.model.5,
                                which.model = "cdc", E = "E_d.c",
                                term = "fi", null.model.formula = f.cdc.model.null,
                                model.2.formula = f.cdc.model.2, 
                                model.3.formula = f.cdc.model.3,
                                model.4.formula = f.cdc.model.4,
                                state.codes.omit = states.out.w )

## ---o--- ##


## (6.3) Save Results ##

# dir.create( "04-Tables-Figures/05-west-analysis" )
# dir.create( "04-Tables-Figures/05-west-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/05-west-analysis/cdc-results" )

inla_results_save( res.jhu.w, path = paste0( "04-Tables-Figures/05-west-analysis/jhu-results/" ),
                   tag = "jhu-w", map.decomp.height = 2060  )

inla_results_save( res.cdc.w, path = paste0( "04-Tables-Figures/05-west-analysis/cdc-results/" ),
                   tag = "cdc-w", map.decomp.height = 2060  )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (8.0) INLA Analysis: Generate Forest Plot Plotting All Results From Models Above ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# prepare lists with the INLA model objects
list.mods.jhu <- list(  res.jhu.ne$model.fit,
                        res.jhu.so$model.fit,
                        res.jhu.mw$model.fit,
                        res.jhu.w$model.fit,
                        res.jhu.1$model.fit )

list.mods.cdc <- list(  res.cdc.ne$model.fit,
                        res.cdc.so$model.fit,
                        res.cdc.mw$model.fit,
                        res.cdc.w$model.fit,
                        res.cdc.1$model.fit )

# names for each of the regions to appear on plot
region.char <- c("Census Region: Northeast", "Census Region: South", 
                 "Census Region: Midwest", "Census Region: West", "National")


# forest plot jhu
( f.plot.jhu <- f_plot_inla( models = list.mods.jhu,
                             term = "fi", 
                             x.label = unname( TeX( "Standardized Mortality Ratio (SMR)" ) ),
                             title = "Data Source: Johns Hopkins",
                             x.names = region.char,
                             x.axis.limits = 1.5,
                             limits.direction = "right") +
    theme( plot.margin=unit(c( 0.1,0.3,0.5,0.1), "cm" ) ) + # unit(c(top, right, bottom, left), units)
    theme( plot.title = element_text( size = 15.2, color = "grey44", hjust = 1 ) ) ) # color title text and right-align it

# forest plot cdc
( f.plot.cdc <- f_plot_inla( models = list.mods.cdc,
                             term = "fi", 
                             x.label = unname( TeX( "Standardized Mortality Ratio (SMR)" ) ),
                             title = "Data Source: CDC",
                             x.names = region.char,
                             x.axis.limits = 1.3,
                             limits.direction = "right" ) +
    theme( plot.margin=unit(c( 0.1,0.3,0.1,0.1), "cm" ) ) +
    theme( plot.title = element_text( size = 15.2, color = "grey44", hjust = 1 ) ) )

# arrange the forest plots and save
ggarrange( f.plot.jhu + theme( axis.title.x = element_blank() ), f.plot.cdc, nrow = 2 )  

ggsave( "04-Tables-Figures/01-main-analysis/forest-plot-both-all-analyses.png",
        width = 6.5, height = 8.5)

# save individual forest plots
ggsave( "04-Tables-Figures/01-main-analysis/jhu-results/forest-plot-jhu-main-analysis.png",
        width = 6.5, height = 8.5, plot = f.plot.jhu )

ggsave( "04-Tables-Figures/01-main-analysis/cdc-results/forest-plot-cdc-main-analysis.png",
        width = 6.5, height = 8.5, plot = f.plot.cdc )


# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (9.0) Check Correlation of Spatial Rand. Eff. w/ Covariates ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


marg.spatial <- res.jhu.1$model.fit$marginals.random$re.s[1:nrow( d.poly )]

spatial.mean <- sapply( marg.spatial, function(x) inla.emarginal( exp, x ) )

d.spatial.check <- cbind( d, spatial.mean )

## (9.1) Columns to Use for Matrix ##

these.columns <- c( "spatial.mean", "pop", "pop.density", "fi.perc.20", 
                    "median.age","perc.female", 
                    "perc.native","perc.hisp","perc.black","perc.asian",
                    "perc.nh.white", "perc.fb", "pct.emp.trade","pct.emp.trans",
                    "ed.1less.than.hspct","ed.5college.plus.pct","poverty.rate",
                    "med.hhinc", "avg.hhsize", "gini.index", "elec.2020.margin", 
                    "inc.prop.jhu","perc.vaccinated","ratio.pop.edp",
                    "health.index" )

## ---o--- ##


## (9.2) Generate Pearson Correlation Matrix ##

cor.fct <- cov.wt( d.spatial.check[, these.columns ], cor = TRUE )

# round entries
cor.mat <- cor.fct$cor %>% round( digits = 2)

# subset of columns
# keep only columns for FI and mort. rate
keep <- which( colnames( cor.mat ) %in% c( "spatial.mean", "pop.density" ) )

cor.mat.sub <- cor.mat[, keep ] # matrix with only spatial random effect and pop density

# spatially structured random effect appears to be orthogonal

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

