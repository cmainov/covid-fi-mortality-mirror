###-------------------------------------------------------------
###   SENSITIVITY-02: ADDITIONAL COVARIATES
###------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we fit the conditional autoregressive models in INLA (with the same specification as in
# "R/04-primary-analysis-inla.R") but we add county-level poverty rate, % Black, % Hispanic, and the unemployment rate
#  as covariates in the model that were not included in the primary analysis.
# 
# 
#
# INPUT DATA FILE: "03-Data-Rodeo/analytic-data.rds"
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
d <- readRDS( "03-Data-Rodeo/01-analytic-data.rds" ) 

# pull county-level shapefiles from `tigris`
t <- counties( cb = TRUE ) # USA Census tract shapefiles


## all data -- without removal of states from state sensitivity analysis ##
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


## Data with removal of UT, NE, and FL in JHU data ##

# merge geometries to our dataset
d.poly.jhu <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  dplyr::select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d, by = c( "fips" ) ) %>%
  filter( fips %in% d$fips ) %>%
  filter( !state.code %in% c( "49", "31", "12" ) ) %>% # remove potentially underreporting states (NE, FL, and UT)
  mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
          exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
          sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
          exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
          sir.jhu.state = cases.jhu.state / exp.cases.jhu.state ) # SIR computation, state-level (JHU data)
## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
## for each dataset we feed through `car_inla_analysis`

## cdc dataset with removal of states from state sensitivity analysis ##
these.out.cdc <- readRDS( "03-Data-Rodeo/03-state-names-vector-remove-sensitivity-analyses.rds" ) # potentially underreporting states identified in cdc data

d.poly.cdc <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  dplyr::select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d, by = c( "fips" ) ) %>%
  filter( fips %in% d$fips ) %>%
  filter( !state %in% these.out.cdc ) %>% # remove potentially underreporting states
  mutate(
    E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using CDC data)
    exp.cases.cdc = ceiling( mean( .$inc.prop.cdc, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using CDC data)
    sir.cdc = cases.cdc / exp.cases.cdc, # SIR computation, county-level (CDC data)
    exp.cases.cdc.state = ceiling( mean( .$inc.prop.cdc.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using CDC data)
    sir.cdc.state = cases.cdc.state / exp.cases.cdc.state ) # SIR computation, state-level (JHU data)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (2.0) Formula Specification ###
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

states.fips.out.jhu.sens <- unique( c( d[ which( d$state %in% c( "UT", "NE", "FL")), "state.code" ],
                                       "69", "60", "66", "78", "15", "02", "72" ) )

states.fips.out.cdc.sens <- unique( c( d[ which( d$state %in% these.out.cdc), "state.code" ],
                                       "69", "60", "66", "78", "15", "02", "72" ) )


# jhu
pov.jhu.1.sens.states <- car_inla_analysis( d.jhu = d.poly.jhu, d.cdc = NULL,
                                            formula.jhu = f.jhu.model.5.pov,
                                            which.model = "jhu", E = "E_d.j",
                                            term = "fi", null.model.formula = f.jhu.model.null,
                                            model.2.formula = f.jhu.model.2, 
                                            model.3.formula = f.jhu.model.3,
                                            model.4.formula = f.jhu.model.4,
                                            state.codes.omit = states.fips.out.jhu.sens,
                                            na.text = "Incomplete or Omitted Data",
                                            x.legend.pos = 0.78 )

pov.jhu.1 <- car_inla_analysis( d.jhu = d.poly, d.cdc = NULL,
                                formula.jhu = f.jhu.model.5.pov,
                                which.model = "jhu", E = "E_d.j",
                                term = "fi", null.model.formula = f.jhu.model.null,
                                model.2.formula = f.jhu.model.2, 
                                model.3.formula = f.jhu.model.3,
                                model.4.formula = f.jhu.model.4,
                                state.codes.omit = states.fips.out,
                                na.text = "Incomplete or Omitted Data",
                                x.legend.pos = 0.78 )


# cdc
pov.cdc.1.sens.states <- car_inla_analysis( d.jhu = NULL, d.cdc = d.poly.cdc,
                                   formula.cdc = f.cdc.model.5.pov,
                                   which.model = "cdc", E = "E_d.c",
                                   term = "fi", null.model.formula = f.cdc.model.null,
                                   model.2.formula = f.cdc.model.2, 
                                   model.3.formula = f.cdc.model.3,
                                   model.4.formula = f.cdc.model.4,
                                   state.codes.omit = states.fips.out.cdc.sens,
                                   na.text = "Incomplete or Omitted Data",
                                   x.legend.pos = 0.78 )

pov.cdc.1 <- car_inla_analysis( d.cdc = d.poly, 
                                formula.cdc = f.cdc.model.5.pov,
                                which.model = "cdc", E = "E_d.c",
                                term = "fi", null.model.formula = f.cdc.model.null,
                                model.2.formula = f.cdc.model.2, 
                                model.3.formula = f.cdc.model.3,
                                model.4.formula = f.cdc.model.4,
                                state.codes.omit = states.fips.out  )
## ---o--- ##


## (3.2) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/01-main-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/01-main-analysis/cdc-results" )

inla_results_save( pov.jhu.1.sens.states, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/01-main-analysis/jhu/" ),
                   tag = "add-vars-jhu-main", map.decomp.height = 2060  )

inla_results_save( pov.cdc.1.sens.states, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/01-main-analysis/cdc/" ),
                   tag = "add-vars-cdc-main", map.decomp.height = 2060  )

## ---o--- ##


## (3.3) Make Additional Supplementary Table with all Covariate Parameters ##

( jhu.all.cov <- pov.jhu.1.sens.states$model.summary$fixed %>%
    data.frame() %>%
    select( mean, `X0.025quant`, `X0.975quant` ) %>%
    filter( !row_number() %in% c(1) ) %>% # remove intercept row
    mutate( across( everything(), .fns = ~ formatC(signif(., digits=3), digits=3, format="fg", flag="#") ) ) ) %>%
  write.table( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/01-main-analysis/jhu/all-covariate-parameters-jhu.txt",
               sep = "," )

( cdc.all.cov <- pov.cdc.1.sens.states$model.summary$fixed %>%
    data.frame() %>%
    select( mean, `X0.025quant`, `X0.975quant` ) %>%
    filter( !row_number() %in% c(1) ) %>% # remove intercept row
    mutate( across( everything(), .fns = ~ formatC(signif(., digits=3), digits=3, format="fg", flag="#") ) ) ) %>%
  write.table( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/01-main-analysis/cdc/all-covariate-parameters-cdc.txt",
               sep = "," )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (4.0) INLA Analysis: Northeast Census Region ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

d.ne <- d %>%
  filter( state %in% c( "ME", "NH", "VT", "MA", "RI", "CT", "NY", "NJ", "PA" ) )

states.out.ne <- unlist( c( unique( d[ d$state %notin%  c( "ME", "NH", "VT", "MA", "RI", "CT", "NY", "NJ", "PA" ), "state.code" ] ),
                            "69", "60", "66", "78", "15", "02", "72" ) )

# wrangle/join shapefile data before plotting

# jhu dataset
d.poly.jhu.ne <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  dplyr::select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d, by = c( "fips" ) ) %>%
  filter( fips %in% d.ne$fips ) %>%
  filter( !state.code %in% c( "49", "31", "12" ) ) %>% # remove potentially underreporting states
  mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
          exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
          sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
          exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
          sir.jhu.state = cases.jhu.state / exp.cases.jhu.state ) # SIR computation, state-level (JHU data)
## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
## for each dataset we feed through `car_inla_analysis`

# cdc dataset
these.out.cdc <- readRDS( "03-Data-Rodeo/03-state-names-vector-remove-sensitivity-analyses.rds" )

d.poly.cdc.ne <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  dplyr::select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d, by = c( "fips" ) ) %>%
  filter( fips %in% d.ne$fips ) %>%
  filter( !state %in% these.out.cdc ) %>% # remove potentially underreporting states
  mutate(
    E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using CDC data)
    exp.cases.cdc = ceiling( mean( .$inc.prop.cdc, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using CDC data)
    sir.cdc = cases.cdc / exp.cases.cdc, # SIR computation, county-level (CDC data)
    exp.cases.cdc.state = ceiling( mean( .$inc.prop.cdc.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using CDC data)
    sir.cdc.state = cases.cdc.state / exp.cases.cdc.state ) # SIR computation, state-level (JHU data)


## ---o--- ##


## (4.2) Call INLA ##

# jhu
pov.jhu.ne <- car_inla_analysis( d.jhu = d.poly.jhu.ne, 
                                 formula.jhu = f.jhu.model.5.pov,
                                 which.model = "jhu", E = "E_d.c",
                                 term = "fi", null.model.formula = f.jhu.model.null,
                                 model.2.formula = f.jhu.model.2, 
                                 model.3.formula = f.jhu.model.3,
                                 model.4.formula = f.jhu.model.4,
                                 state.codes.omit = states.fips.out  )

# cdc
pov.cdc.ne <- car_inla_analysis( d.cdc = d.poly.cdc.ne, 
                                 formula.cdc = f.cdc.model.5.pov,
                                 which.model = "cdc", E = "E_d.c",
                                 term = "fi", null.model.formula = f.cdc.model.null,
                                 model.2.formula = f.cdc.model.2, 
                                 model.3.formula = f.cdc.model.3,
                                 model.4.formula = f.cdc.model.4,
                                 state.codes.omit = states.fips.out  )
## ---o--- ##


## (4.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/02-northeast-analysis" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/02-northeast-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/02-northeast-analysis/cdc-results" )

inla_results_save( pov.jhu.ne, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/02-northeast-analysis/jhu/" ),
                   tag = "add-vars-jhu-ne", map.decomp.width = (4124)*0.7, map.decomp.height = (3808)*0.7 )

inla_results_save( pov.cdc.ne, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/02-northeast-analysis/cdc/" ),
                   tag = "add-vars-cdc-ne", map.decomp.width = (4124)*0.7, map.decomp.height = (3808)*0.7 )

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
# jhu dataset
d.poly.jhu.so <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  dplyr::select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d, by = c( "fips" ) ) %>%
  filter( fips %in% d.so$fips ) %>%
  filter( !state.code %in% c( "49", "31", "12" ) ) %>% # remove potentially underreporting states
  mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
          exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
          sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
          exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
          sir.jhu.state = cases.jhu.state / exp.cases.jhu.state ) # SIR computation, state-level (JHU data)
## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
## for each dataset we feed through `car_inla_analysis`

d.poly.cdc.so <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  dplyr::select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d, by = c( "fips" ) ) %>%
  filter( fips %in% d.so$fips ) %>%
  filter( !state %in% these.out.cdc ) %>% # remove potentially underreporting states
  mutate(
    E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using CDC data)
    exp.cases.cdc = ceiling( mean( .$inc.prop.cdc, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using CDC data)
    sir.cdc = cases.cdc / exp.cases.cdc, # SIR computation, county-level (CDC data)
    exp.cases.cdc.state = ceiling( mean( .$inc.prop.cdc.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using CDC data)
    sir.cdc.state = cases.cdc.state / exp.cases.cdc.state ) # SIR computation, state-level (JHU data)


## ---o--- ##


## (5.2) Call INLA ##

# jhu
pov.jhu.so <- car_inla_analysis( d.jhu = d.poly.jhu.so, 
                                 formula.jhu = f.jhu.model.5.pov,
                                 which.model = "jhu", E = "E_d.c",
                                 term = "fi", null.model.formula = f.jhu.model.null,
                                 model.2.formula = f.jhu.model.2, 
                                 model.3.formula = f.jhu.model.3,
                                 model.4.formula = f.jhu.model.4,
                                 state.codes.omit = states.fips.out  )

# cdc
pov.cdc.so <- car_inla_analysis( d.cdc = d.poly.cdc.so, 
                                 formula.cdc = f.cdc.model.5.pov,
                                 which.model = "cdc", E = "E_d.c",
                                 term = "fi", null.model.formula = f.cdc.model.null,
                                 model.2.formula = f.cdc.model.2, 
                                 model.3.formula = f.cdc.model.3,
                                 model.4.formula = f.cdc.model.4,
                                 state.codes.omit = states.fips.out  )

## ---o--- ##


## (5.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/03-south-analysis" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/03-south-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/03-south-analysis/cdc-results" )

inla_results_save( pov.jhu.so, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/03-south-analysis/jhu/" ),
                   tag = "add-vars-jhu-so" )

inla_results_save( pov.cdc.so, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/03-south-analysis/cdc/" ),
                   tag = "add-vars-cdc-so" )

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
d.poly.jhu.mw <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  dplyr::select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d, by = c( "fips" ) ) %>%
  filter( fips %in% d.mw$fips ) %>%
  filter( !state.code %in% c( "49", "31", "12" ) ) %>% # remove potentially underreporting states
  mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
          exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
          sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
          exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
          sir.jhu.state = cases.jhu.state / exp.cases.jhu.state ) # SIR computation, state-level (JHU data)
## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
## for each dataset we feed through `car_inla_analysis`

d.poly.cdc.mw <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  dplyr::select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d, by = c( "fips" ) ) %>%
  filter( fips %in% d.mw$fips ) %>%
  filter( !state %in% these.out.cdc ) %>% # remove potentially underreporting states
  mutate(
    E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using CDC data)
    exp.cases.cdc = ceiling( mean( .$inc.prop.cdc, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using CDC data)
    sir.cdc = cases.cdc / exp.cases.cdc, # SIR computation, county-level (CDC data)
    exp.cases.cdc.state = ceiling( mean( .$inc.prop.cdc.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using CDC data)
    sir.cdc.state = cases.cdc.state / exp.cases.cdc.state ) # SIR computation, state-level (JHU data)

## ---o--- ##


## (6.2) Call INLA ##

# jhu
pov.jhu.mw <- car_inla_analysis( d.jhu = d.poly.jhu.mw, 
                                 formula.jhu = f.jhu.model.5.pov,
                                 which.model = "jhu", E = "E_d.c",
                                 term = "fi", null.model.formula = f.jhu.model.null,
                                 model.2.formula = f.jhu.model.2, 
                                 model.3.formula = f.jhu.model.3,
                                 model.4.formula = f.jhu.model.4,
                                 state.codes.omit = states.fips.out  )

# cdc
pov.cdc.mw <- car_inla_analysis( d.cdc = d.poly.cdc.mw, 
                                 formula.cdc = f.cdc.model.5.pov,
                                 which.model = "cdc", E = "E_d.c",
                                 term = "fi", null.model.formula = f.cdc.model.null,
                                 model.2.formula = f.cdc.model.2, 
                                 model.3.formula = f.cdc.model.3,
                                 model.4.formula = f.cdc.model.4,
                                 state.codes.omit = states.fips.out  )

## ---o--- ##


## (6.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/04-midwest-analysis" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/04-midwest-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/04-midwest-analysis/cdc-results" )

inla_results_save( pov.jhu.mw, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/04-midwest-analysis/jhu/" ),
                   tag = "add-vars-jhu-mw" )

inla_results_save( pov.cdc.mw, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/04-midwest-analysis/cdc/" ),
                   tag = "add-vars-cdc-mw" )

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
d.poly.jhu.w <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  dplyr::select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d, by = c( "fips" ) ) %>%
  filter( fips %in% d.w$fips ) %>%
  filter( !state.code %in% c( "49", "31", "12" ) ) %>% # remove potentially underreporting states
  mutate( E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
          exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
          sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
          exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
          sir.jhu.state = cases.jhu.state / exp.cases.jhu.state ) # SIR computation, state-level (JHU data)
## note: SIR and expected deaths counts are dataset-specific, which is why we compute them
## for each dataset we feed through `car_inla_analysis`

d.poly.cdc.w <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  dplyr::select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d, by = c( "fips" ) ) %>%
  filter( fips %in% d.w$fips ) %>%
  filter( !state %in% these.out.cdc ) %>% # remove potentially underreporting states
  mutate(
    E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using CDC data)
    exp.cases.cdc = ceiling( mean( .$inc.prop.cdc, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using CDC data)
    sir.cdc = cases.cdc / exp.cases.cdc, # SIR computation, county-level (CDC data)
    exp.cases.cdc.state = ceiling( mean( .$inc.prop.cdc.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using CDC data)
    sir.cdc.state = cases.cdc.state / exp.cases.cdc.state ) # SIR computation, state-level (JHU data)

## ---o--- ##


## (7.2) Call INLA ##

# jhu
pov.jhu.w <- car_inla_analysis( d.jhu = d.poly.jhu.w, 
                                 formula.jhu = f.jhu.model.5.pov,
                                 which.model = "jhu", E = "E_d.c",
                                 term = "fi", null.model.formula = f.jhu.model.null,
                                 model.2.formula = f.jhu.model.2, 
                                 model.3.formula = f.jhu.model.3,
                                 model.4.formula = f.jhu.model.4,
                                 state.codes.omit = states.fips.out  )

# cdc
pov.cdc.w <- car_inla_analysis( d.cdc = d.poly.cdc.w, 
                                 formula.cdc = f.cdc.model.5.pov,
                                 which.model = "cdc", E = "E_d.c",
                                 term = "fi", null.model.formula = f.cdc.model.null,
                                 model.2.formula = f.cdc.model.2, 
                                 model.3.formula = f.cdc.model.3,
                                 model.4.formula = f.cdc.model.4,
                                 state.codes.omit = states.fips.out  )

## ---o--- ##


## (6.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/05-west-analysis" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/05-west-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/05-west-analysis/cdc-results" )

inla_results_save( pov.jhu.w, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/05-west-analysis/jhu/" ),
                   tag = "add-vars-jhu-w" )

inla_results_save( pov.cdc.w, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/05-west-analysis/cdc/" ),
                   tag = "add-vars-cdc-w" )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (8.0) INLA Analysis: Generate Forest Plot Plotting All Results From Models Above ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# prepare lists with the INLA model objects
list.mods.jhu <- list(  pov.jhu.ne$model.fit,
                        pov.jhu.so$model.fit,
                        pov.jhu.mw$model.fit,
                        pov.jhu.w$model.fit,
                        pov.jhu.1.sens.states$model.fit )

list.mods.cdc <- list(  pov.cdc.ne$model.fit,
                        pov.cdc.so$model.fit,
                        pov.cdc.mw$model.fit,
                        pov.cdc.w$model.fit,
                        pov.cdc.1$model.fit )

# names for each of the regions to appear on plot
region.char <- c("Census Region: Northeast", "Census Region: South", 
                 "Census Region: Midwest", "Census Region: West", "National")


# forest plot jhu
( f.plot.jhu <- f_plot_inla( models = list.mods.jhu,
                             term = "fi", 
                             x.label = unname( TeX( "Standardized Mortality Ratio (SMR)" ) ),
                             title = "Data Source: Johns Hopkins",
                             x.names = region.char,
                             x.axis.limits = 1.55,
                             limits.direction = "right" ) +
    theme( plot.margin=unit(c( 0.1,0.3,0.5,0.1), "cm" ) ) + # unit(c(top, right, bottom, left), units)
    theme( plot.title = element_text( size = 12, color = "grey44", hjust = 1 ) ) ) # color title text and right-align it

# forest plot cdc
( f.plot.cdc <- f_plot_inla( models = list.mods.cdc,
                             term = "fi", 
                             x.label = unname( TeX( "Standardized Mortality Ratio (SMR)" ) ),
                             title = "Data Source: CDC",
                             x.names = region.char,
                             x.axis.limits = 1.35,
                             limits.direction = "right" ) +
    theme( plot.margin=unit(c( 0.5,0.3,0.1,0.1), "cm" ) ) +
    theme( plot.title = element_text( size = 12, color = "grey44", hjust = 1 ) ) )

# arrange the forest plots and save
ggarrange( f.plot.jhu + theme( axis.title.x = element_blank() ), f.plot.cdc, nrow = 2 )  

ggsave( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/add-vars-forest-plot-both-all-analyses.png",
        width = 6.5, height = 8 )

# save individual forest plots
ggsave( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/01-main-analysis/jhu/forest-plot-jhu-add-vars-analysis.png",
        width = 6.5, height = 8.5, plot = f.plot.jhu )

ggsave( "04-Tables-Figures/06-sensitivity-analyses/02-add-addtl-vars/01-main-analysis/cdc/forest-plot-cdc-add-vars-analysis.png",
        width = 6.5, height = 8.5, plot = f.plot.cdc )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
