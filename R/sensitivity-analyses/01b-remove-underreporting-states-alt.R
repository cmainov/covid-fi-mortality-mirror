##----------------------------------------------------------------------------------------------------
###   SENSITIVITY-01B: REMOVE UNDERREPORTING STATES IN CDC AND JHU DATA--NO COVARIATE SELECTION
###---------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we fit the conditional autoregressive models in INLA (with the same specification as in
# "R/04-primary-analysis-inla.R") but we remove underreporting states. To identify underreporting states,
# we implement an algorithm whereby the states with COVID-19 mortality rates (based on CDC data) that are
# more than 2x smaller than the COVID-19 mortality rates based on the JHU data are identified and removed 
#
# In addition, we also remove Utah, Nebraska, and Florida counties from the analysis 
# using JHU data as they appeared suspicious in the chloropleths (see R/03-descriptives).
# 
# This analysis is the same as "01-remove-underreporting-states.R"  but instead of implementing covariate selection with the 
# Watanabe–Akaike information criterion (WAIC), it extracts the predictors from the main analysis 
# for the final model in each of the national stratified analyses and uses those predictors in the model 
# specification. Thus, covariate selection is turned off. NOTE: we used the first file, 
# "01-remove-underreporting-states.R", for this sensitivity analysis and wrote the code for this alternate
# file in case we need it in the future. Be mindful that running the code in this file overwrites all the
# figures and tables that are generated in "01-remove-underreporting-states.R".
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



### (2.0) Formula Specification ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (2.1) Models with Full Set of Covariates ##

# jhu data
f.jhu <- deaths.jhu.adj ~ I(fi.perc.20/4) + pct.emp.trans + no.vehic + 
  disability + no.health.insur + perc.female + 
  perc.nh.white + perc.native + pop.density + ed.1less.than.hspct + 
  ed.5college.plus.pct + pct.emp.trade + median.age + perc.asian +
  perc.vaccinated + gini.index + avg.hhsize + elec.2020.margin +
  ratio.pop.edp + health.index + sir.jhu + urb.cat.code + 
  sir.jhu.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +
  f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid") + f( fips, model = "iid" )

# cdc data
f.cdc <- deaths.cdc ~ I(fi.perc.20/4) + pct.emp.trans + no.vehic + 
  disability + no.health.insur + perc.female + 
  perc.nh.white + perc.native + pop.density + ed.1less.than.hspct + 
  ed.5college.plus.pct + pct.emp.trade + median.age + perc.asian +
  perc.vaccinated + gini.index + avg.hhsize + elec.2020.margin +
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
  f( re.s, model = "besag", graph = g.cdc, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

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
            cdc.state.avg = mean( age.adj.mort.rate, na.rm = T ), # compute mean state mortality rate using CDC data
            jhu.state.avg = mean( mort.rate.jhu.adj, na.rm = T ),     # compute mean state mortality rate using JHU data
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
          E_d.c = ceiling( mean( .$age.adj.mort.rate, na.rm = TRUE )*pop /100000 ), # expected death count for offset term (using CDC data)
          exp.cases.jhu = ceiling( mean( .$inc.prop.jhu, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using JHU data)
          exp.cases.cdc = ceiling( mean( .$inc.prop.cdc, na.rm = TRUE )*pop /100000 ), # expected case count for SIR computation (using CDC data)
          sir.jhu = cases.jhu / exp.cases.jhu, # SIR computation, county-level (JHU data)
          sir.cdc = cases.cdc / exp.cases.cdc, # SIR computation, county-level (CDC data)
          exp.cases.jhu.state = ceiling( mean( .$inc.prop.jhu.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using JHU data)
          exp.cases.cdc.state = ceiling( mean( .$inc.prop.cdc.state, na.rm = TRUE )*state.pop /100000 ), # expected case count for state-level SIR computation (using CDC data)
          sir.jhu.state = cases.jhu.state / exp.cases.jhu.state, # SIR computation, state-level (JHU data)
          sir.cdc.state= cases.cdc.state / exp.cases.cdc.state ) # SIR computation, state-level (JHU data)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (4.0) Fit the Model on National Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# states to omit from maps
states.out <- d.states[ which( d.states$diff.flag == 1 ), "state" ]
states.fips.out <- unique( c( d[ which( d$state %in% states.out), "state.code" ],
                              "69", "60", "66", "78", "15", "02", "72" ) )

## (4.1) Call INLA ##

# specify formula (we specify the same formula with the same predictors that were selected for in the main analysis)
this.file.cdc.1 <- "/Users/mainovieytesca/Documents/GitHub/COVID-FI-Mortality/04-Tables-Figures/01-main-analysis/cdc-results/model-selection-log.txt" # this obtains the file with the list of predictors that were included in the final model

# clunky way of reading in the predictors .txt file and preparing it for creating the formula
f.1.preds.cdc <- readLines( this.file.cdc.1 ) %>%
  noquote( . ) %>%
  str_extract(., "(?<=R:\\s).*") %>% noquote( . ) %>%
  str_replace(., '\\"', '') %>%
  .[ !is.na(.) ]

f.cdc.a <- as.formula( paste0( "deaths.cdc ~ I(fi.perc.20/4) +", paste0( f.1.preds.cdc, collapse = '+' ),
                   "+ f( re.s, model = 'besag', graph = g.cdc, scale.model = T ) + f( state, model = 'iid') + f( fips, model = 'iid' )" ) )

res.cdc.1.sens <- car_inla_analysis( d.jhu = d.poly, d.cdc = d.poly.cdc,
                                     formula.jhu = f.jhu, formula.cdc = f.cdc,
                                     which.model = "cdc", null.model.formula = f.cdc.model.null,
                                     selection.criterion = "waic",
                                     criterion.threshold = 5,
                                     model.2.formula = f.cdc.model.2, 
                                     model.3.formula = f.cdc.model.3, 
                                     term = "fi", E = "E_d.c", 
                                     state.codes.omit = states.fips.out)

## ---o--- ##


## (4.2) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/cdc-remove-underreporting-states" )

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/01-main-analysis" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/01-main-analysis/jhu-results" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/01-main-analysis/cdc-results" )

inla_results_save( res.cdc.1.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/01-main-analysis/" ),
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

# specify formula (we specify the same formula with the same predictors that were selected for in the main analysis)
this.file.cdc.2 <- "/Users/mainovieytesca/Documents/GitHub/COVID-FI-Mortality/04-Tables-Figures/02-northeast-analysis/cdc-results/model-selection-log.txt" # this obtains the file with the list of predictors that were included in the final model

# clunky way of reading in the predictors .txt file and preparing it for creating the formula
f.2.preds.cdc <- readLines( this.file.cdc.2 ) %>%
  noquote( . ) %>%
  str_extract(., "(?<=R:\\s).*") %>% noquote( . ) %>%
  str_replace(., '\\"', '') %>%
  .[ !is.na(.) ]

f.cdc.b <- as.formula( paste0( "deaths.cdc ~ I(fi.perc.20/4) +", paste0( f.2.preds.cdc, collapse = '+' ),
                               "+ f( re.s, model = 'besag', graph = g.cdc, scale.model = T ) + f( state, model = 'iid') + f( fips, model = 'iid' )" ) )

# cdc
res.cdc.ne.sens <- car_inla_analysis( d.jhu = d.poly.cdc.ne, d.cdc = d.poly.cdc.ne,
                                 formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                 term = "fi", which.model = "cdc", E = "E_d.c", 
                                 selection.criterion = "waic",
                                 criterion.threshold = 5,
                                 null.model.formula = f.cdc.model.null,
                                 model.2.formula = f.cdc.model.2, 
                                 model.3.formula = f.cdc.model.3,
                                 state.codes.omit = states.out.ne )

## ---o--- ##


## (5.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/02-northeast-analysis" )

inla_results_save( res.cdc.ne.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/02-northeast-analysis/" ),
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

# specify formula (we specify the same formula with the same predictors that were selected for in the main analysis)
this.file.cdc.3 <- "/Users/mainovieytesca/Documents/GitHub/COVID-FI-Mortality/04-Tables-Figures/03-south-analysis/cdc-results/model-selection-log.txt" # this obtains the file with the list of predictors that were included in the final model

# clunky way of reading in the predictors .txt file and preparing it for creating the formula
f.3.preds.cdc <- readLines( this.file.cdc.3 ) %>%
  noquote( . ) %>%
  str_extract(., "(?<=R:\\s).*") %>% noquote( . ) %>%
  str_replace(., '\\"', '') %>%
  .[ !is.na(.) ]

f.cdc.c <- as.formula( paste0( "deaths.cdc ~ I(fi.perc.20/4) +", paste0( f.3.preds.cdc, collapse = '+' ),
                               "+ f( re.s, model = 'besag', graph = g.cdc, scale.model = T ) + f( state, model = 'iid') + f( fips, model = 'iid' )" ) )

# cdc
res.cdc.so.sens <- car_inla_analysis( d.jhu = d.poly.cdc.so, d.cdc = d.poly.cdc.so,
                                 formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                 term = "fi", which.model = "cdc", E = "E_d.c", 
                                 selection.criterion = "waic",
                                 criterion.threshold = 5,
                                 null.model.formula = f.cdc.model.null,
                                 model.2.formula = f.cdc.model.2, 
                                 model.3.formula = f.cdc.model.3,
                                 state.codes.omit = states.out.so )

## ---o--- ##


## (6.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/03-south-analysis" )

inla_results_save( res.cdc.so.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/03-south-analysis/" ),
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

# specify formula (we specify the same formula with the same predictors that were selected for in the main analysis)
this.file.cdc.4 <- "/Users/mainovieytesca/Documents/GitHub/COVID-FI-Mortality/04-Tables-Figures/04-midwest-analysis/cdc-results/model-selection-log.txt" # this obtains the file with the list of predictors that were included in the final model

# clunky way of reading in the predictors .txt file and preparing it for creating the formula
f.4.preds.cdc <- readLines( this.file.cdc.4 ) %>%
  noquote( . ) %>%
  str_extract(., "(?<=R:\\s).*") %>% noquote( . ) %>%
  str_replace(., '\\"', '') %>%
  .[ !is.na(.) ]

f.cdc.d <- as.formula( paste0( "deaths.cdc ~ I(fi.perc.20/4) +", paste0( f.4.preds.cdc, collapse = '+' ),
                               "+ f( re.s, model = 'besag', graph = g.cdc, scale.model = T ) + f( state, model = 'iid') + f( fips, model = 'iid' )" ) )

# cdc
res.cdc.mw.sens <- car_inla_analysis( d.jhu = d.poly.cdc.mw, d.cdc = d.poly.cdc.mw,
                                 formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                 term = "fi", which.model = "cdc", E = "E_d.c", 
                                 selection.criterion = "waic",
                                 criterion.threshold = 5,
                                 null.model.formula = f.cdc.model.null,
                                 model.2.formula = f.cdc.model.2, 
                                 model.3.formula = f.cdc.model.3,
                                 state.codes.omit = states.out.mw )

## ---o--- ##


## (7.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/04-midwest-analysis" )

inla_results_save( res.cdc.mw.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/04-midwest-analysis/" ),
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

# specify formula (we specify the same formula with the same predictors that were selected for in the main analysis)
this.file.cdc.5 <- "/Users/mainovieytesca/Documents/GitHub/COVID-FI-Mortality/04-Tables-Figures/05-west-analysis/cdc-results/model-selection-log.txt" # this obtains the file with the list of predictors that were included in the final model

# clunky way of reading in the predictors .txt file and preparing it for creating the formula
f.5.preds.cdc <- readLines( this.file.cdc.5 ) %>%
  noquote( . ) %>%
  str_extract(., "(?<=R:\\s).*") %>% noquote( . ) %>%
  str_replace(., '\\"', '') %>%
  .[ !is.na(.) ]

f.cdc.e <- as.formula( paste0( "deaths.cdc ~ I(fi.perc.20/4) +", paste0( f.5.preds.cdc, collapse = '+' ),
                               "+ f( re.s, model = 'besag', graph = g.cdc, scale.model = T ) + f( state, model = 'iid') + f( fips, model = 'iid' )" ) )


# cdc
res.cdc.w.sens <- car_inla_analysis( d.jhu = d.poly.cdc.w, d.cdc = d.poly.cdc.w,
                                formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                term = "fi", which.model = "cdc", E = "E_d.c", 
                                selection.criterion = "waic",
                                criterion.threshold = 5,
                                null.model.formula = f.cdc.model.null,
                                model.2.formula = f.cdc.model.2, 
                                model.3.formula = f.cdc.model.3,
                                state.codes.omit = states.out.w )

## (8.3) Save Results ##

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/05-west-analysis" )

inla_results_save( res.cdc.w.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/05-west-analysis/" ),
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
ggsave( "04-Tables-Figures/06-sensitivity-analyses/01-remove-underreporting-cdc/underreport-cdc-forest-plot-both-all-analyses.png",
        width = 6.5, height = 8.5 )

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

# specify formula (we specify the same formula with the same predictors that were selected for in the main analysis)
this.file.jhu.1 <- "/Users/mainovieytesca/Documents/GitHub/COVID-FI-Mortality/04-Tables-Figures/01-main-analysis/jhu-results/model-selection-log.txt" # this obtains the file with the list of predictors that were included in the final model

# clunky way of reading in the predictors .txt file and preparing it for creating the formula
f.1.preds.jhu <- readLines( this.file.jhu.1 ) %>%
  noquote( . ) %>%
  str_extract(., "(?<=R:\\s).*") %>% noquote( . ) %>%
  str_replace(., '\\"', '') %>%
  .[ !is.na(.) ]

f.jhu.a <- as.formula( paste0( "deaths.jhu.adj ~ I(fi.perc.20/4) +", paste0( f.1.preds.jhu, collapse = '+' ),
                   "+ f( re.s, model = 'besag', graph = g, scale.model = T ) + f( state, model = 'iid') + f( fips, model = 'iid' )" ) )

res.jhu.1.sens <- car_inla_analysis( d.jhu = d.poly.no.ne.fl.ut, d.cdc = d.poly.no.ne.fl.ut,
                                     formula.jhu = f.jhu, formula.cdc = f.cdc,
                                     selection.criterion = "waic",
                                     criterion.threshold = 5,
                                     which.model = "jhu", null.model.formula = f.jhu.model.null,
                                     model.2.formula = f.jhu.model.2, 
                                     model.3.formula = f.jhu.model.3, 
                                     term = "fi", E = "E_d.j",
                                     state.codes.omit = states.fips.out.jhu )

## ---o--- ##


## (9.2) Save Results ##

 # dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu" )
 # dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu/01-main-analysis" )

inla_results_save( res.jhu.1.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu/01-main-analysis/" ),
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

# specify formula (we specify the same formula with the same predictors that were selected for in the main analysis)
this.file.jhu.2 <- "/Users/mainovieytesca/Documents/GitHub/COVID-FI-Mortality/04-Tables-Figures/02-northeast-analysis/jhu-results/model-selection-log.txt" # this obtains the file with the list of predictors that were included in the final model

# clunky way of reading in the predictors .txt file and preparing it for creating the formula
f.2.preds.jhu <- readLines( this.file.jhu.2 ) %>%
  noquote( . ) %>%
  str_extract(., "(?<=R:\\s).*") %>% noquote( . ) %>%
  str_replace(., '\\"', '') %>%
  .[ !is.na(.) ]

f.jhu.b <- as.formula( paste0( "deaths.jhu.adj ~ I(fi.perc.20/4) +", paste0( f.2.preds.jhu, collapse = '+' ),
                   "+ f( re.s, model = 'besag', graph = g, scale.model = T ) + f( state, model = 'iid') + f( fips, model = 'iid' )" ) )

# jhu
res.jhu.ne.sens <- car_inla_analysis( d.jhu = d.poly.jhu.ne, d.cdc = d.poly.jhu.ne,
                                      formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                      selection.criterion = "waic",
                                      criterion.threshold = 5,
                                      term = "fi", which.model = "jhu", E = "E_d.c", 
                                      null.model.formula = f.jhu.model.null,
                                      model.2.formula = f.jhu.model.2, 
                                      model.3.formula = f.jhu.model.3,
                                      state.codes.omit = states.out.ne )

## ---o--- ##


## (10.3) Save Results ##
 # dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu/")
 # dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu/02-northeast-analysis" )

inla_results_save( res.jhu.ne.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu/02-northeast-analysis/" ),
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

# specify formula (we specify the same formula with the same predictors that were selected for in the main analysis)
this.file.jhu.3 <- "/Users/mainovieytesca/Documents/GitHub/COVID-FI-Mortality/04-Tables-Figures/03-south-analysis/jhu-results/model-selection-log.txt" # this obtains the file with the list of predictors that were included in the final model

# clunky way of reading in the predictors .txt file and preparing it for creating the formula
f.3.preds.jhu <- readLines( this.file.jhu.3 ) %>%
  noquote( . ) %>%
  str_extract(., "(?<=R:\\s).*") %>% noquote( . ) %>%
  str_replace(., '\\"', '') %>%
  .[ !is.na(.) ]

f.jhu.c <- as.formula( paste0( "deaths.jhu.adj ~ I(fi.perc.20/4) +", paste0( f.3.preds.jhu, collapse = '+' ),
                   "+ f( re.s, model = 'besag', graph = g, scale.model = T ) + f( state, model = 'iid') + f( fips, model = 'iid' )" ) )

# jhu
res.jhu.so.sens <- car_inla_analysis( d.jhu = d.poly.jhu.so, d.cdc = d.poly.jhu.so,
                                      formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                      selection.criterion = "waic",
                                      criterion.threshold = 5,
                                      term = "fi", which.model = "jhu", E = "E_d.c", 
                                      null.model.formula = f.jhu.model.null,
                                      model.2.formula = f.jhu.model.2, 
                                      model.3.formula = f.jhu.model.3,
                                      state.codes.omit = states.out.so )

## ---o--- ##


## (11.3) Save Results ##

 # dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu/03-south-analysis" )

inla_results_save( res.jhu.so.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu/03-south-analysis/" ),
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

# specify formula (we specify the same formula with the same predictors that were selected for in the main analysis)
this.file.jhu.4 <- "/Users/mainovieytesca/Documents/GitHub/COVID-FI-Mortality/04-Tables-Figures/04-midwest-analysis/jhu-results/model-selection-log.txt" # this obtains the file with the list of predictors that were included in the final model

# clunky way of reading in the predictors .txt file and preparing it for creating the formula
f.4.preds.jhu <- readLines( this.file.jhu.4 ) %>%
  noquote( . ) %>%
  str_extract(., "(?<=R:\\s).*") %>% noquote( . ) %>%
  str_replace(., '\\"', '') %>%
  .[ !is.na(.) ]

f.jhu.d <- as.formula( paste0( "deaths.jhu.adj ~ I(fi.perc.20/4) +", paste0( f.4.preds.jhu, collapse = '+' ),
                   "+ f( re.s, model = 'besag', graph = g, scale.model = T ) + f( state, model = 'iid') + f( fips, model = 'iid' )" ) )

# jhu
res.jhu.mw.sens <- car_inla_analysis( d.jhu = d.poly.jhu.mw, d.cdc = d.poly.jhu.mw,
                                      formula.jhu = f.jhu, formula.cdc = f.cdc, 
                                      term = "fi", which.model = "jhu", E = "E_d.c", 
                                      selection.criterion = "waic",
                                      criterion.threshold = 5,
                                      null.model.formula = f.jhu.model.null,
                                      model.2.formula = f.jhu.model.2, 
                                      model.3.formula = f.jhu.model.3,
                                      state.codes.omit = states.out.mw )

## ---o--- ##


## (12.3) Save Results ##

 # dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu/04-midwest-analysis" )

inla_results_save( res.jhu.mw.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu/04-midwest-analysis/" ),
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

# specify formula (we specify the same formula with the same predictors that were selected for in the main analysis)
this.file.jhu.5 <- "/Users/mainovieytesca/Documents/GitHub/COVID-FI-Mortality/04-Tables-Figures/05-west-analysis/jhu-results/model-selection-log.txt" # this obtains the file with the list of predictors that were included in the final model

# clunky way of reading in the predictors .txt file and preparing it for creating the formula
f.5.preds.jhu <- readLines( this.file.jhu.5 ) %>%
  noquote( . ) %>%
  str_extract(., "(?<=R:\\s).*") %>% noquote( . ) %>%
  str_replace(., '\\"', '') %>%
  .[ !is.na(.) ]

f.jhu.e <- as.formula( paste0( "deaths.jhu.adj ~ I(fi.perc.20/4) +", paste0( f.5.preds.jhu, collapse = '+' ),
                   "+ f( re.s, model = 'besag', graph = g, scale.model = T ) + f( state, model = 'iid') + f( fips, model = 'iid' )" ) )

# jhu
res.jhu.w.sens <- car_inla_analysis( d.jhu = d.poly.jhu.w, d.cdc = d.poly.jhu.w,
                                     formula.jhu = f.jhu, formula.cdc = f.jhu, 
                                     term = "fi", which.model = "jhu", E = "E_d.c", 
                                     selection.criterion = "waic",
                                     criterion.threshold = 5,
                                     null.model.formula = f.jhu.model.null,
                                     model.2.formula = f.jhu.model.2, 
                                     model.3.formula = f.jhu.model.3,
                                     state.codes.omit = states.out.w )

## (13.3) Save Results ##

 # dir.create( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu/05-west-analysis" )

inla_results_save( res.jhu.w.sens, path = paste0( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu/05-west-analysis/" ),
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
ggsave( "04-Tables-Figures/06-sensitivity-analyses/02-remove-underreporting-jhu/underreport-jhu-forest-plot-both-all-analyses.png",
        width = 6.5, height = 8.5)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

