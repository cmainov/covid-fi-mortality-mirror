##--------------------------------------------------------------------
###   SENSITIVITY-01: REMOVE UNDERREPORTING STATES IN CDC AND JHU DATA
###-------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we fit the conditional autoregressive models in INLA (with the same specification as in
# "R/04-primary-analysis-inla.R") but we remove outlier states. To identify underreporting states,
# we implement an algorithm whereby the states with COVID-19 mortality rates (based on CDC data) that are
# more than 2x smaller than the COVID-19 mortality rates based on the JHU data are identified and removed 
#
# In addition, we also remove Utah, Nebraska, and Florida counties from the analysis 
# using JHU data as they appeared suspicious in the chloropleths and state-level residuals
# (see R/03-descriptives) and R/04-primary-analysis-inla).
# 
#
# INPUT DATA FILE: "03-Data-Rodeo/01-analytic-data.rds"
#
#
#
# Resources: see "R/04-primary-analysis-inla.R"
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
d <- readRDS( "03-Data-Rodeo/01-analytic-data.rds" ) %>%
  mutate( mort.rate.jhu.adj = amr*100000 )

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


### (3.0) Remove States That Appear to be Underreporting (CDC data)###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# threshold for flagging states that might be underreporting deaths 
f <- 2

# evaluate for under-reporting
( d.states <- d %>%
    group_by( state.code ) %>%
    mutate( 
      crude.cdc.mort.rate = ( raw.deaths.cdc / pop )*100000,
      crude.jhu.mort.rate = ( raw.deaths.cdc / pop )*100000,
      cdc.state.avg = mean( crude.cdc.mort.rate, na.rm = T ), # compute mean state mortality rate using CDC data
      jhu.state.avg = mean( crude.jhu.mort.rate, na.rm = T ),     # compute mean state mortality rate using JHU data
      state.avg.diff = jhu.state.avg - cdc.state.avg,       # compute difference of state mortality rates
      diff.flag = ifelse( jhu.state.avg > f*abs( cdc.state.avg ) |
                            cdc.state.avg > f*abs( jhu.state.avg ), 1, 0 ) ) %>%  # if the difference between the two is greater than f, we will flag it as underreporting
    ungroup() %>%
    distinct( state, jhu.state.avg, cdc.state.avg, diff.flag ) %>% data.frame() )

d.states <- data.frame( cbind( aggregate( mort.rate.jhu.adj~state, data = d, FUN = "mean"),
                               aggregate( age.adj.mort.rate~state, data = d, FUN = "mean")) ) %>%
  select( -state.1)

colnames( d.states ) <-c( "state", "jhu.state.avg", "cdc.state.avg")

d.states<- d.states %>%
  mutate( state.avg.diff = jhu.state.avg - cdc.state.avg,
          diff.flag = ifelse( jhu.state.avg > f*abs( cdc.state.avg ) |
                                cdc.state.avg > f*abs( jhu.state.avg ), 1, 0 ) ) %>% 
  distinct( state, jhu.state.avg, cdc.state.avg, diff.flag ) %>% data.frame() 

table( d.states$diff.flag) # 15 state levels removed

# check who and save for later sensitivity analyses outside this script
( d.states %>%
    filter( diff.flag == 1 ) %>%
    select( state ) )  %>% # remember not to consider HI and PR
  unlist() %>%
  saveRDS( "03-Data-Rodeo/03-state-names-vector-remove-sensitivity-analyses.rds" )

# create data set that removes those states from the analysis
d.cdc <- d %>%
  group_by( state ) %>%
  mutate( cdc.state.avg = mean( age.adj.mort.rate, na.rm = T ),
          jhu.state.avg = mean( mort.rate.jhu.adj, na.rm = T ),
          state.avg.diff = jhu.state.avg - cdc.state.avg,
          diff.flag = ifelse( jhu.state.avg > f*abs( cdc.state.avg ) |
                                cdc.state.avg > f*abs( jhu.state.avg ), 1, 0 ) ) %>%
  ungroup() %>%
  subset(., diff.flag != 1 ) 

# create final dataset that filters the underreporting states out
d.poly.cdc <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d.cdc, by = c( "fips" ) ) %>%
  filter( fips %in% d.cdc$fips ) %>%
  mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
          E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop / 100000 ), # expected death count for offset term (using CDC data)
          exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop / 100000 ), # expected case count for SIR computation (using JHU data)
          exp.cases.cdc = ceiling( mean( .$inc.prop.cdc, na.rm = TRUE )*pop / 100000 ), # expected case count for SIR computation (using CDC data)
          sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
          sir.cdc = cases.cdc / exp.cases.cdc, # SIR computation, county-level (CDC data)
          exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
          exp.cases.cdc.state = ceiling( mean( .$inc.prop.cdc.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using CDC data)
          sir.jhu.state = cases.jhu.state / exp.cases.jhu.state, # SIR computation, state-level (JHU data)
          sir.cdc.state = cases.cdc.state / exp.cases.cdc.state ) # SIR computation, state-level (JHU data)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (4.0) Fit the Model on National Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# states to omit from maps
states.out <- d.states[ which( d.states$diff.flag == 1 ), "state" ]
states.fips.out <- unique( c( d[ which( d$state %in% states.out), "state.code" ],
                              "69", "60", "66", "78", "15", "02", "72" ) )

## (4.1) Call INLA ##

res.cdc.1.sens <- car_inla_analysis( d.jhu = NULL, d.cdc = d.poly.cdc,
                                     formula.cdc = f.cdc.model.5,
                                     which.model = "cdc", E = "E_d.c",
                                     term = "fi", null.model.formula = f.cdc.model.null,
                                     model.2.formula = f.cdc.model.2, 
                                     model.3.formula = f.cdc.model.3,
                                     model.4.formula = f.cdc.model.4,
                                     state.codes.omit = states.fips.out,
                                     na.text = "Incomplete or Omitted Data",
                                     x.legend.pos = 0.78 )

 ## ---o--- ##


## (4.2) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/cdc-remove-underreporting-states" )

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/01-main-analysis" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/01-main-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/01-main-analysis/cdc-results" )

inla_results_save( res.cdc.1.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/01-main-analysis/" ),
                   tag = "underreport-cdc-main" )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (5.0) INLA Analysis: Northeast Census Region ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (5.1) Subset Data Based on Census Region ##

d.ne <- d.cdc %>%
  filter( state %in% c( "ME", "NH", "VT", "MA", "RI", "CT", "NY", "NJ", "PA" ) )

states.out.ne <- unlist( c( unique( d[ d$state %notin%  c( "ME", "NH", "VT", "MA", "RI", "CT", "NY", "NJ", "PA" ), "state.code" ] ),
                            "69", "60", "66", "78", "15", "02", "72" ) )

# wrangle/join shapefile data before plotting
( d.poly.cdc.ne <- t %>%
    mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
    select( -state.code ) %>%
    dplyr::rename( fips = GEOID ) %>%
    left_join( ., d.ne, by = c( "fips" ) ) %>%
    filter( fips %in% d.ne$fips ) %>%
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


## (5.2) Call INLA ##

# cdc
res.cdc.ne.sens <- car_inla_analysis( d.jhu = NULL, d.cdc = d.poly.cdc.ne,
                                      formula.cdc = f.cdc.model.5,
                                      which.model = "cdc", E = "E_d.c",
                                      term = "fi", null.model.formula = f.cdc.model.null,
                                      model.2.formula = f.cdc.model.2, 
                                      model.3.formula = f.cdc.model.3,
                                      model.4.formula = f.cdc.model.4,
                                      state.codes.omit = states.out.ne,
                                      na.text = "Incomplete or Omitted Data",
                                      x.legend.pos = 0.78 )

## ---o--- ##


## (5.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/02-northeast-analysis" )

inla_results_save( res.cdc.ne.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/02-northeast-analysis/" ),
                   tag = "underreport-cdc-ne", map.decomp.height = 2060  )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (6.0) INLA Analysis: South Census Region ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (6.1) Subset Data Based on Census Region ##

d.so <- d.cdc %>%
  filter( state %in% c( "MD", "DE", "DC", "WV", "VA", "NC", "KY", "TN", 
                        "SC","GA", "FL", "AL", "MS", "AR", "LA", "OK", "TX" ))

states.out.so <- unlist( c( unique( d[ d$state %notin%  c( "MD", "DE", "DC", "WV", "VA", "NC", "KY", "TN", 
                                                           "SC","GA", "FL", "AL", "MS", "AR", "LA", "OK", "TX" ), "state.code" ] ),
                            "69", "60", "66", "78", "15", "02", "72" ) )

# wrangle/join shapefile data before plotting
( d.poly.cdc.so <- t %>%
    mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
    select( -state.code ) %>%
    dplyr::rename( fips = GEOID ) %>%
    left_join( ., d.so, by = c( "fips" ) ) %>%
    filter( fips %in% d.so$fips ) %>%
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

# cdc
res.cdc.so.sens <- car_inla_analysis( d.jhu = NULL, d.cdc = d.poly.cdc.so,
                                      formula.cdc = f.cdc.model.5,
                                      which.model = "cdc", E = "E_d.c",
                                      term = "fi", null.model.formula = f.cdc.model.null,
                                      model.2.formula = f.cdc.model.2, 
                                      model.3.formula = f.cdc.model.3,
                                      model.4.formula = f.cdc.model.4,
                                      state.codes.omit = states.out.so,
                                      na.text = "Incomplete or Omitted Data",
                                      x.legend.pos = 0.78 )

## ---o--- ##


## (6.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/03-south-analysis" )

inla_results_save( res.cdc.so.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/03-south-analysis/" ),
                   tag = "underreport-cdc-so", map.decomp.height = 2060  )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (7.0) INLA Analysis: Midwest Census Region ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (7.1) Subset Data Based on Census Region ##

d.mw <- d.cdc %>%
  filter( state %in% c( "OH", "IN", "MI", "IL", "WI", "MN", "IA", "MO","ND", "SD","NE", "KS"))

states.out.mw <- unlist( c( unique( d[ d$state %notin%  c( "OH", "IN", "MI", "IL", "WI", "MN", "IA", "MO","ND", "SD","NE", "KS"), "state.code" ] ),
                            "69", "60", "66", "78", "15", "02", "72" ) )

# wrangle/join shapefile data before plotting
( d.poly.cdc.mw <- t %>%
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


## (7.2) Call INLA ##

# cdc
res.cdc.mw.sens <- car_inla_analysis( d.jhu = NULL, d.cdc = d.poly.cdc.mw,
                                      formula.cdc = f.cdc.model.5,
                                      which.model = "cdc", E = "E_d.c",
                                      term = "fi", null.model.formula = f.cdc.model.null,
                                      model.2.formula = f.cdc.model.2, 
                                      model.3.formula = f.cdc.model.3,
                                      model.4.formula = f.cdc.model.4,
                                      state.codes.omit = states.out.mw,
                                      na.text = "Incomplete or Omitted Data",
                                      x.legend.pos = 0.78 )

## ---o--- ##


## (7.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/04-midwest-analysis" )

inla_results_save( res.cdc.mw.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/04-midwest-analysis/" ),
                   tag = "underreport-cdc-mw", map.decomp.height = 2060  )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (8.0) INLA Analysis: West Census Region ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (8.1) Subset Data Based on Census Region ##

d.w <- d.cdc %>%
  filter( state %in% c( "NM", "CO", "WY", "MT", "ID", "UT", "AZ", "NV", "CA", "OR", "WA"))

states.out.w <- unlist( c( unique( d[ d$state %notin%  c( "NM", "CO", "WY", "MT", "ID", "UT", "AZ", "NV", "CA", "OR", "WA"), "state.code" ] ),
                           "69", "60", "66", "78", "15", "02", "72" ) )

# wrangle/join shapefile data before plotting
( d.poly.cdc.w <- t %>%
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


## (8.2) Call INLA ##

# cdc
res.cdc.w.sens <- car_inla_analysis( d.jhu = NULL, d.cdc = d.poly.cdc.w,
                                     formula.cdc = f.cdc.model.5,
                                     which.model = "cdc", E = "E_d.c",
                                     term = "fi", null.model.formula = f.cdc.model.null,
                                     model.2.formula = f.cdc.model.2, 
                                     model.3.formula = f.cdc.model.3,
                                     model.4.formula = f.cdc.model.4,
                                     state.codes.omit = states.out.w,
                                     na.text = "Incomplete or Omitted Data",
                                     x.legend.pos = 0.78 )

## (8.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/05-west-analysis" )

inla_results_save( res.cdc.w.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/05-west-analysis/" ),
                   tag = "underreport-cdc-w", map.decomp.height = 2060  )


## ---o--- ##


# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (9.0) INLA Analysis: Generate Forest Plot Plotting All Results From Models Above (CDC Underreporting Models) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# prepare lists with the INLA model objects
list.mods.cdc.sens <- list(  res.cdc.ne.sens$model.fit,
                             res.cdc.so.sens$model.fit,
                             res.cdc.mw.sens$model.fit,
                             res.cdc.w.sens$model.fit,
                             res.cdc.1.sens$model.fit )

# names for each of the regions to appear on plot
region.char <- c("Census Region: Northeast", "Census Region: South", 
                 "Census Region: Midwest", "Census Region: West", "National")


# forest plot cdc
( f.plot.cdc <- f_plot_inla( models = list.mods.cdc.sens,
                             term = "fi", 
                             x.label = unname( TeX( "Standardized Mortality Ratio (SMR)" ) ),
                             title = "Data Source: CDC",
                             x.names = region.char,
                             x.axis.limits = 1.2,
                             limits.direction = "right" ) +
    theme( plot.margin=unit(c( 0.1,0.3,0.1,0.1), "cm" ) ) +
    theme( plot.title = element_text( size = 15.2, color = "grey44", hjust = 1 ) ) )


# save to incorporate into a combination plot with the map decomposition of the JHU main model
ggsave( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/cdc/underreport-cdc-forest-plot-both-all-analyses.png",
        width = 6.5, height = 8.5)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------





### (9.0) Remove UT, NE, and FL from JHU Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

states.fips.out.jhu <- unique( c( d[ which( d$state %in% c( "UT", "NE", "FL")), "state.code" ],
                                  "69", "60", "66", "78", "15", "02", "72" ) )

# merge geometries to our dataset
d.poly.no.ne.fl.ut <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  dplyr::select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d, by = c( "fips" ) ) %>%
  filter( !state.code %in% c( "49", "31", "12" ) )

d.poly.no.ne.fl.ut <- d.poly.no.ne.fl.ut %>%
  mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
          E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using CDC data)
          exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
          exp.cases.cdc = ceiling( mean( .$inc.prop.cdc, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using CDC data)
          sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
          sir.cdc = cases.cdc / exp.cases.cdc, # SIR computation, county-level (CDC data)
          exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
          exp.cases.cdc.state = ceiling( mean( .$inc.prop.cdc.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using CDC data)
          sir.jhu.state = cases.jhu.state / exp.cases.jhu.state, # SIR computation, state-level (JHU data)
          sir.cdc.state= cases.cdc.state / exp.cases.cdc.state ) # SIR computation, state-level (JHU data)


## (9.1) Call INLA ##

res.jhu.1.sens <- car_inla_analysis( d.jhu = d.poly.no.ne.fl.ut, d.cdc = NULL,
                                     formula.jhu = f.jhu.model.5, E = "E_d.j",
                                     which.model = "jhu", null.model.formula = f.jhu.model.null,
                                     model.2.formula = f.jhu.model.2, 
                                     model.3.formula = f.jhu.model.3, 
                                     model.4.formula = f.jhu.model.4, 
                                     term = "fi", state.codes.omit = states.fips.out.jhu,
                                     na.text = "Incomplete or Omitted Data",
                                     x.legend.pos = 0.78 )

## ---o--- ##


## (9.2) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/jhu/01-main-analysis" )

inla_results_save( res.jhu.1.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/jhu/01-main-analysis/" ),
                   tag = "no-ne-fl-ut-jhu-main" )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (10.0) INLA Analysis: Northeast Census Region ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (10.1) Subset Data Based on Census Region ##

d.ne <- d.poly.no.ne.fl.ut %>%
  filter( state %in% c( "ME", "NH", "VT", "MA", "RI", "CT", "NY", "NJ", "PA" ) )

states.out.ne <- unlist( c( unique( d[ d$state %notin%  c( "ME", "NH", "VT", "MA", "RI", "CT", "NY", "NJ", "PA" ), "state.code" ] ),
                            "69", "60", "66", "78", "15", "02", "72" ) )

# wrangle/join shapefile data before plotting
( d.poly.jhu.ne <- t %>%
    mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
    select( -state.code ) %>%
    dplyr::rename( fips = GEOID ) %>%
    left_join( ., st_drop_geometry( d.ne ), by = c( "fips" ) ) %>%
    filter( fips %in% d.ne$fips ) %>%
    mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
            E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using jhu data)
            exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
            exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using jhu data)
            sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
            sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (jhu data)
            exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
            exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using jhu data)
            sir.jhu.state = cases.jhu.state / exp.cases.jhu.state, # SIR computation, state-level (JHU data)
            sir.jhu.state= cases.jhu.state / exp.cases.jhu.state ) ) # SIR computation, state-level (JHU data)
## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
## for each dataset we feed through `car_inla_analysis`


## ---o--- ##


## (10.2) Call INLA ##

# jhu
res.jhu.ne.sens <- car_inla_analysis( d.jhu = d.poly.jhu.ne, d.cdc = NULL,
                                      formula.jhu = f.jhu.model.5, E = "E_d.j",
                                      which.model = "jhu", null.model.formula = f.jhu.model.null,
                                      model.2.formula = f.jhu.model.2, 
                                      model.3.formula = f.jhu.model.3, 
                                      model.4.formula = f.jhu.model.4, 
                                      term = "fi",
                                      state.codes.omit = states.out.ne,
                                      na.text = "Incomplete or Omitted Data",
                                      x.legend.pos = 0.78 )

## ---o--- ##


## (10.3) Save Results ##
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/jhu/")
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/jhu/02-northeast-analysis" )

inla_results_save( res.jhu.ne.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/jhu/02-northeast-analysis/" ),
                   tag = "underreport-jhu-ne", map.decomp.height = 2060  )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (11.0) INLA Analysis: South Census Region ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (11.1) Subset Data Based on Census Region ##

d.so <- d.poly.no.ne.fl.ut %>%
  filter( state %in% c( "MD", "DE", "DC", "WV", "VA", "NC", "KY", "TN", 
                        "SC","GA", "FL", "AL", "MS", "AR", "LA", "OK", "TX" ))

states.out.so <- unlist( c( unique( d[ d$state %notin%  c( "MD", "DE", "DC", "WV", "VA", "NC", "KY", "TN", 
                                                           "SC","GA", "FL", "AL", "MS", "AR", "LA", "OK", "TX" ), "state.code" ] ),
                            "69", "60", "66", "78", "15", "02", "72" ) )

# wrangle/join shapefile data before plotting
( d.poly.jhu.so <- t %>%
    mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
    select( -state.code ) %>%
    dplyr::rename( fips = GEOID ) %>%
    left_join( ., st_drop_geometry( d.so ), by = c( "fips" ) ) %>%
    filter( fips %in% d.so$fips ) %>%
    mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
            E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using jhu data)
            exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
            exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using jhu data)
            sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
            sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (jhu data)
            exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
            exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using jhu data)
            sir.jhu.state = cases.jhu.state / exp.cases.jhu.state, # SIR computation, state-level (JHU data)
            sir.jhu.state= cases.jhu.state / exp.cases.jhu.state ) ) # SIR computation, state-level (JHU data)
## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
## for each dataset we feed through `car_inla_analysis`


## ---o--- ##


## (11.2) Call INLA ##

# jhu
res.jhu.so.sens <- car_inla_analysis( d.jhu = d.poly.jhu.so, d.cdc = NULL,
                                      formula.jhu = f.jhu.model.5, E = "E_d.j",
                                      which.model = "jhu", null.model.formula = f.jhu.model.null,
                                      model.2.formula = f.jhu.model.2, 
                                      model.3.formula = f.jhu.model.3, 
                                      model.4.formula = f.jhu.model.4, 
                                      term = "fi",
                                      state.codes.omit = states.out.so,
                                      na.text = "Incomplete or Omitted Data",
                                      x.legend.pos = 0.78 )

## ---o--- ##


## (11.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/jhu/03-south-analysis" )

inla_results_save( res.jhu.so.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/jhu/03-south-analysis/" ),
                   tag = "underreport-jhu-so", map.decomp.height = 2060  )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (12.0) INLA Analysis: Midwest Census Region ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (12.1) Subset Data Based on Census Region ##

d.mw <- d.poly.no.ne.fl.ut %>%
  filter( state %in% c( "OH", "IN", "MI", "IL", "WI", "MN", "IA", "MO","ND", "SD","NE", "KS"))

states.out.mw <- unlist( c( unique( d[ d$state %notin%  c( "OH", "IN", "MI", "IL", "WI", "MN", "IA", "MO","ND", "SD","NE", "KS"), "state.code" ] ),
                            "69", "60", "66", "78", "15", "02", "72" ) )

# wrangle/join shapefile data before plotting
( d.poly.jhu.mw <- t %>%
    mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
    select( -state.code ) %>%
    dplyr::rename( fips = GEOID ) %>%
    left_join( ., st_drop_geometry( d.mw ), by = c( "fips" ) ) %>%
    filter( fips %in% d.mw$fips ) %>%
    mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
            E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using jhu data)
            exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
            exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using jhu data)
            sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
            sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (jhu data)
            exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
            exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using jhu data)
            sir.jhu.state = cases.jhu.state / exp.cases.jhu.state, # SIR computation, state-level (JHU data)
            sir.jhu.state= cases.jhu.state / exp.cases.jhu.state ) ) # SIR computation, state-level (JHU data)
## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
## for each dataset we feed through `car_inla_analysis`

## ---o--- ##


## (12.2) Call INLA ##

# jhu
res.jhu.mw.sens <- car_inla_analysis( d.jhu = d.poly.jhu.mw, d.cdc = d.poly.jhu.mw,
                                      formula.jhu = f.jhu.model.5, E = "E_d.j",
                                      which.model = "jhu", null.model.formula = f.jhu.model.null,
                                      model.2.formula = f.jhu.model.2, 
                                      model.3.formula = f.jhu.model.3, 
                                      model.4.formula = f.jhu.model.4, 
                                      term = "fi",
                                      state.codes.omit = states.out.mw,
                                      na.text = "Incomplete or Omitted Data",
                                      x.legend.pos = 0.78 )

## ---o--- ##


## (12.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/jhu/04-midwest-analysis" )

inla_results_save( res.jhu.mw.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/jhu/04-midwest-analysis/" ),
                   tag = "underreport-jhu-mw", map.decomp.height = 2060  )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (13.0) INLA Analysis: West Census Region ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (13.1) Subset Data Based on Census Region ##

d.w <- d.poly.no.ne.fl.ut %>%
  filter( state %in% c( "NM", "CO", "WY", "MT", "ID", "UT", "AZ", "NV", "CA", "OR", "WA"))

states.out.w <- unlist( c( unique( d[ d$state %notin%  c( "NM", "CO", "WY", "MT", "ID", "UT", "AZ", "NV", "CA", "OR", "WA"), "state.code" ] ),
                           "69", "60", "66", "78", "15", "02", "72" ) )

# wrangle/join shapefile data before plotting
( d.poly.jhu.w <- t %>%
    mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
    select( -state.code ) %>%
    dplyr::rename( fips = GEOID ) %>%
    left_join( ., st_drop_geometry( d.w ), by = c( "fips" ) ) %>%
    filter( fips %in% d.w$fips ) %>%
    mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
            E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using jhu data)
            exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
            exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using jhu data)
            sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
            sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (jhu data)
            exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
            exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using jhu data)
            sir.jhu.state = cases.jhu.state / exp.cases.jhu.state, # SIR computation, state-level (JHU data)
            sir.jhu.state= cases.jhu.state / exp.cases.jhu.state ) ) # SIR computation, state-level (JHU data)
## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
## for each dataset we feed through `car_inla_analysis`


## ---o--- ##


## (13.2) Call INLA ##

# jhu
res.jhu.w.sens <- car_inla_analysis( d.jhu = d.poly.jhu.w, d.cdc = d.poly.jhu.w,
                                     formula.jhu = f.jhu.model.5, E = "E_d.j",
                                     which.model = "jhu", null.model.formula = f.jhu.model.null,
                                     model.2.formula = f.jhu.model.2, 
                                     model.3.formula = f.jhu.model.3, 
                                     model.4.formula = f.jhu.model.4, 
                                     term = "fi",
                                     state.codes.omit = states.out.w,
                                     na.text = "Incomplete or Omitted Data",
                                     x.legend.pos = 0.78 )

## (13.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/jhu/05-west-analysis" )

inla_results_save( res.jhu.w.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/jhu/05-west-analysis/" ),
                   tag = "underreport-jhu-w", map.decomp.height = 2060  )


## ---o--- ##


# ---------------------------------------------------------------------------------------------------------------------------------------------------------


### (9.0) INLA Analysis: Generate Forest Plot Plotting All Results From Models Above (jhu Underreporting Models) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# prepare lists with the INLA model objects
list.mods.jhu.sens <- list(  res.jhu.ne.sens$model.fit,
                             res.jhu.so.sens$model.fit,
                             res.jhu.mw.sens$model.fit,
                             res.jhu.w.sens$model.fit,
                             res.jhu.1.sens$model.fit )

# names for each of the regions to appear on plot
region.char <- c("Census Region: Northeast", "Census Region: South", 
                 "Census Region: Midwest", "Census Region: West", "National")


# forest plot jhu
( f.plot.jhu <- f_plot_inla( models = list.mods.jhu.sens,
                             term = "fi", 
                             x.label = unname( TeX( "Standardized Mortality Ratio (SMR)" ) ),
                             title = "Data Source: Johns Hopkins",
                             x.names = region.char,
                             x.axis.limits = 1.2,
                             limits.direction = "right" ) +
    theme( plot.margin=unit(c( 0.1,0.3,0.1,0.1), "cm" ) ) +
    theme( plot.title = element_text( size = 15.2, color = "grey44", hjust = 1 ) ) )


# save to incorporate into a combination plot with the map decomposition of the JHU main model
ggsave( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-states/jhu/underreport-jhu-forest-plot-both-all-analyses.png",
        width = 6.5, height = 8.5)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
