###-------------------------------------------------------------
###   SENSITIVITY-04: NO COVARIATE SELECTION
###------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we fit the conditional autoregressive models in INLA (with the same specification as in
# "R/04-primary-analysis-inla.R") but we do not use any covariate/model selection alogirthm and information
# criterion. Instead, we fit the "full" models with the entire set of considered covariates.
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

# helper functions
source( "R/utils.R" )


### (1.0) Data Import and Preparation ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# read-in data and omit states not included in analysis
d <- readRDS( "03-Data-Rodeo/01-analytic-data.rds" ) 

# pull county-level shapefiles from `tigris`
t <- counties( cb = TRUE ) # USA Census tract shapefiles

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

# cdc dataset
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

# jhu data
f.jhu <- deaths.jhu.adj ~ I(fi.perc.20/4) + pct.emp.trans + unemp.rate + no.vehic + 
  disability + no.health.insur + perc.black + perc.female + perc.hisp + 
  perc.nh.white + perc.native + pop.density + ed.1less.than.hspct + 
  ed.5college.plus.pct + pct.emp.trade + median.age + perc.asian +
  poverty.rate + perc.vaccinated + gini.index + avg.hhsize + elec.2020.margin +
  ratio.pop.edp + health.index + sir.jhu + urb.cat.code + 
  sir.jhu.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +
  f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid") + f( fips, model = "iid" )

# cdc data
f.cdc <- deaths.cdc ~ I(fi.perc.20/4) + pct.emp.trans + unemp.rate + no.vehic + 
  disability + no.health.insur + perc.black + perc.female + perc.hisp + 
  perc.nh.white + perc.native + pop.density + ed.1less.than.hspct + 
  ed.5college.plus.pct + pct.emp.trade + median.age + perc.asian +
  poverty.rate + perc.vaccinated + gini.index + avg.hhsize + elec.2020.margin +
  ratio.pop.edp + health.index + sir.cdc + urb.cat.code + 
  sir.cdc.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +
  f( re.s, model = "besag", graph = g.cdc, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )
## ---o--- ##


## (2.2) Models with Random Effects and Offset Terms Only ("Null/baseline Model") ##

# jhu
f.jhu.model.null <- deaths.jhu.adj ~ f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

#cdc
f.cdc.model.null <- deaths.cdc ~ f( re.s, model = "besag", graph = g.cdc, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

## ---o--- ##


## (2.3) Models with Random Effects, Offset Terms, and Fixed Effects for Health Index and Median Age ("Basic Model") ##
f.jhu.model.2 <- deaths.jhu.adj ~ health.index + median.age + sir.jhu + f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

f.cdc.model.2 <- deaths.cdc ~ health.index + median.age + sir.cdc + f( re.s, model = "besag", graph = g.cdc, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )


## (2.4) Models with Random Effects, Offset Terms, and Fixed Effects for Health Index and Median Age
## and State-Level Variables ("Basic + State Model") ##

f.jhu.model.3 <- deaths.jhu.adj ~ health.index + median.age + sir.jhu +
  sir.jhu.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +
  f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

f.cdc.model.3 <- deaths.cdc ~ health.index + median.age + sir.cdc +
  sir.cdc.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +
  f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (3.0) INLA Analysis: National (Nearly All U.S. Counties) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (3.1) Call INLA ##

# states to remove from maps

states.fips.out.jhu <- unique( c( d[ which( d$state %in% c( "UT", "NE", "FL")), "state.code" ],
                                  "69", "60", "66", "78", "15", "02", "72" ) )

states.fips.out.cdc <- unique( c( d[ which( d$state %in% these.out.cdc), "state.code" ],
                              "69", "60", "66", "78", "15", "02", "72" ) )


# jhu
no.sel.jhu.1 <- car_inla_analysis( d.jhu = d.poly.jhu, d.cdc = d.poly.cdc,
                                    formula.jhu = f.jhu, formula.cdc = f.cdc,
                                    which.model = "jhu", E = "E_d.j",
                                    term = "fi", null.model.formula = f.jhu.model.null,
                                    model.2.formula = f.jhu.model.2, 
                                    model.3.formula = f.jhu.model.3,
                                    selection.criterion = NULL,
                                    state.codes.omit = states.fips.out.jhu,
                                   na.text = "Incomplete or Omitted Data",
                                   x.legend.pos = 0.78 )

# cdc
no.sel.cdc.1 <- car_inla_analysis( d.jhu = d.poly.jhu, d.cdc = d.poly.cdc,
                                    formula.jhu = f.jhu, formula.cdc = f.cdc,
                                    which.model = "cdc",  E = "E_d.c",
                                    term = "fi", null.model.formula = f.cdc.model.null,
                                    model.2.formula = f.cdc.model.2, 
                                    model.3.formula = f.cdc.model.3,
                                    selection.criterion = NULL,
                                    state.codes.omit = states.fips.out.cdc,
                                   na.text = "Incomplete or Omitted Data",
                                   x.legend.pos = 0.78 )

## ---o--- ##


## (3.2) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/01-main-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/01-main-analysis/cdc-results" )

inla_results_save( no.sel.jhu.1, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/01-main-analysis/jhu-results/" ),
                   tag = "no-sel-jhu-main", map.decomp.height = 2060  )

inla_results_save( no.sel.cdc.1, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/01-main-analysis/cdc-results/" ),
                   tag = "no-sel-cdc-main", map.decomp.height = 2060  )

## ---o--- ##


## (3.3) Make Additional Supplementary Table with all Covariate Parameters ##

( jhu.all.cov <- no.sel.jhu.1$model.summary$fixed %>%
    data.frame() %>%
    select( mean, `X0.025quant`, `X0.975quant` ) %>%
    filter( !row_number() %in% c(1) ) %>% # remove intercept row
    mutate( across( everything(), .fns = ~ formatC(signif(., digits=3), digits=3, format="fg", flag="#") ) ) ) %>%
  write.table( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/01-main-analysis/jhu-results/all-covariate-parameters-jhu.txt",
               sep = "," )

( cdc.all.cov <- no.sel.cdc.1$model.summary$fixed %>%
    data.frame() %>%
    select( mean, `X0.025quant`, `X0.975quant` ) %>%
    filter( !row_number() %in% c(1) ) %>% # remove intercept row
    mutate( across( everything(), .fns = ~ formatC(signif(., digits=3), digits=3, format="fg", flag="#") ) ) ) %>%
  write.table( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/01-main-analysis/cdc-results/all-covariate-parameters-cdc.txt",
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
no.sel.jhu.ne <- car_inla_analysis( d.jhu = d.poly.jhu.ne, d.cdc = d.poly.cdc.ne,
                                     formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                     term = "fi", which.model = "jhu", E = "E_d.j", 
                                     null.model.formula = f.jhu.model.null,
                                     model.2.formula = f.jhu.model.2, 
                                     model.3.formula = f.jhu.model.3,
                                     selection.criterion = NULL,
                                     state.codes.omit = states.out.ne,
                                    na.text = "Incomplete or Omitted Data",
                                    x.legend.pos = 0.78 )

# cdc
no.sel.cdc.ne <- car_inla_analysis( d.jhu = d.poly.jhu.ne, d.cdc = d.poly.cdc.ne,
                                     formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                     term = "fi", which.model = "cdc", E = "E_d.c", 
                                     null.model.formula = f.cdc.model.null,
                                     model.2.formula = f.cdc.model.2, 
                                     model.3.formula = f.cdc.model.3,
                                     selection.criterion = NULL,
                                     state.codes.omit = states.out.ne,
                                    na.text = "Incomplete or Omitted Data",
                                    x.legend.pos = 0.78 )
## ---o--- ##


## (4.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/02-northeast-analysis" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/02-northeast-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/02-northeast-analysis/cdc-results" )

inla_results_save( no.sel.jhu.ne, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/02-northeast-analysis/jhu-results/" ),
                   tag = "no-sel-jhu-ne", map.decomp.width = (4124)*0.7, map.decomp.height = (3808)*0.7 )

inla_results_save( no.sel.cdc.ne, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/02-northeast-analysis/cdc-results/" ),
                   tag = "no-sel-cdc-ne", map.decomp.width = (4124)*0.7, map.decomp.height = (3808)*0.7 )

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
no.sel.jhu.so <- car_inla_analysis( d.jhu = d.poly.jhu.so, d.cdc = d.poly.cdc.so,
                                     formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                     term = "fi", which.model = "jhu", E = "E_d.j", 
                                     null.model.formula = f.jhu.model.null,
                                     model.2.formula = f.jhu.model.2, 
                                     model.3.formula = f.jhu.model.3,
                                     selection.criterion = NULL,
                                     state.codes.omit = states.out.so,
                                    na.text = "Incomplete or Omitted Data",
                                    x.legend.pos = 0.78 )
# cdc
no.sel.cdc.so <- car_inla_analysis( d.jhu = d.poly.jhu.so, d.cdc = d.poly.cdc.so,
                                     formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                     term = "fi", which.model = "cdc", E = "E_d.c", 
                                     null.model.formula = f.cdc.model.null,
                                     model.2.formula = f.cdc.model.2, 
                                     model.3.formula = f.cdc.model.3,
                                     selection.criterion = NULL,
                                     state.codes.omit = states.out.so,
                                    na.text = "Incomplete or Omitted Data",
                                    x.legend.pos = 0.78 )

## ---o--- ##


## (5.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/03-south-analysis" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/03-south-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/03-south-analysis/cdc-results" )

inla_results_save( no.sel.jhu.so, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/03-south-analysis/jhu-results/" ),
                   tag = "no-sel-jhu-so" )

inla_results_save( no.sel.cdc.so, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/03-south-analysis/cdc-results/" ),
                   tag = "no-sel-cdc-so" )

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
no.sel.jhu.mw <- car_inla_analysis( d.jhu = d.poly.jhu.mw, d.cdc = d.poly.cdc.mw,
                                     formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                     term = "fi", which.model = "jhu", E = "E_d.j", 
                                     null.model.formula = f.jhu.model.null,
                                     model.2.formula = f.jhu.model.2, 
                                     model.3.formula = f.jhu.model.3,
                                     selection.criterion = NULL,
                                     state.codes.omit = states.out.mw,
                                    na.text = "Incomplete or Omitted Data",
                                    x.legend.pos = 0.78 )

no.sel.cdc.mw <- car_inla_analysis( d.jhu = d.poly.jhu.mw, d.cdc = d.poly.cdc.mw,
                                     formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                     term = "fi", which.model = "cdc", E = "E_d.c", 
                                     null.model.formula = f.cdc.model.null,
                                     model.2.formula = f.cdc.model.2, 
                                     model.3.formula = f.cdc.model.3,
                                     selection.criterion = NULL,
                                     state.codes.omit = states.out.mw,
                                    na.text = "Incomplete or Omitted Data",
                                    x.legend.pos = 0.78 )

## ---o--- ##


## (6.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/04-midwest-analysis" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/04-midwest-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/04-midwest-analysis/cdc-results" )

inla_results_save( no.sel.jhu.mw, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/04-midwest-analysis/jhu-results/" ),
                   tag = "no-sel-jhu-mw" )

inla_results_save( no.sel.cdc.mw, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/04-midwest-analysis/cdc-results/" ),
                   tag = "no-sel-cdc-mw" )

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
no.sel.jhu.w <- car_inla_analysis( d.jhu = d.poly.jhu.w, d.cdc = d.poly.cdc.w,
                                    formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                    term = "fi", which.model = "jhu", E = "E_d.j", 
                                    null.model.formula = f.jhu.model.null,
                                    model.2.formula = f.jhu.model.2, 
                                    model.3.formula = f.jhu.model.3,
                                    selection.criterion = NULL,
                                    state.codes.omit = states.out.w,
                                   na.text = "Incomplete or Omitted Data",
                                   x.legend.pos = 0.78 )
# cdc
no.sel.cdc.w <- car_inla_analysis( d.jhu = d.poly.jhu.w, d.cdc = d.poly.cdc.w,
                                    formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                    term = "fi", which.model = "cdc", E = "E_d.c", 
                                    null.model.formula = f.cdc.model.null,
                                    model.2.formula = f.cdc.model.2, 
                                    model.3.formula = f.cdc.model.3,
                                    selection.criterion = NULL,
                                    state.codes.omit = states.out.w,
                                   na.text = "Incomplete or Omitted Data",
                                   x.legend.pos = 0.78 )

## ---o--- ##


## (6.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/05-west-analysis" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/05-west-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/05-west-analysis/cdc-results" )

inla_results_save( no.sel.jhu.w, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/05-west-analysis/jhu-results/" ),
                   tag = "no-sel-jhu-w" )

inla_results_save( no.sel.cdc.w, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/05-west-analysis/cdc-results/" ),
                   tag = "no-sel-cdc-w" )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (8.0) INLA Analysis: Generate Forest Plot Plotting All Results From Models Above ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# prepare lists with the INLA model objects
list.mods.jhu <- list(  no.sel.jhu.ne$model.fit,
                        no.sel.jhu.so$model.fit,
                        no.sel.jhu.mw$model.fit,
                        no.sel.jhu.w$model.fit,
                        no.sel.jhu.1$model.fit )

list.mods.cdc <- list(  no.sel.cdc.ne$model.fit,
                        no.sel.cdc.so$model.fit,
                        no.sel.cdc.mw$model.fit,
                        no.sel.cdc.w$model.fit,
                        no.sel.cdc.1$model.fit )

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

ggsave( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/no-sel-forest-plot-both-all-analyses.png",
        width = 6.5, height = 8 )

# save individual forest plots
ggsave( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/01-main-analysis/jhu-results/forest-plot-jhu-no-sel-analysis.png",
        width = 6.5, height = 8.5, plot = f.plot.jhu )

ggsave( "04-Tables-Figures/06-sensitivity-analyses/05-no-covariate-selection/01-main-analysis/cdc-results/forest-plot-cdc-no-sel-analysis.png",
        width = 6.5, height = 8.5, plot = f.plot.cdc )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------