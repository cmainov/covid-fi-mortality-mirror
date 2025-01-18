###-------------------------------------------------------------
###   07-RESIDUALS ANALYSIS
###-------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we fit the conditional autoregressive models in INLA and keep
# the same fixed effects but modify the number of random effects in each model. In
# total, we fit 3 models:
# i. one with the fixed effects and only a county-level residual
# ii. one that further adds the state-level random effect (RE)
# iii. one that adds the iCAR RE
#
# We then plot the residuals to examine the behavior of the residuals and if the
# assumption of those residuals can be assumed across the three models.
# We keep the fixed effects constant, including only a food insecurity fixed effect
# and a population density fixed effect, which is the specification that arose 
# from the covariate selection process in the main analysis.
#
# INPUT DATA FILE: "03-Data-Rodeo/analytic-data.rds"
#
#
# Resources: 
# i. Pub with more detailed description of CAR model, its implementation in INLA, and examples: https://doi.org/10.1016/j.sste.2013.07.003
# ii. the R-INLA Project: https://www.r-inla.org/
#
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

library( INLA )
library( tidyverse )
library( spdep ) # for creating neighborhood matrix
library( tigris ) # shapefiles

# helper functions
source( "R/utils.R" )



### (1.0) Data Import and Preparation ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# read-in data and omit states not included in analysis
d <- readRDS( "03-Data-Rodeo/01-analytic-data.rds" ) 

# pull county-level shapefiles from `tigris`
t <- counties( cb = TRUE ) # USA Census tract shapefiles

# merge geometries to our dataset
d.jhu <- t %>%
  mutate( state.code = str_extract( GEOID, "^.{2}") ) %>%
  dplyr::select( -state.code ) %>%
  dplyr::rename( fips = GEOID ) %>%
  left_join( ., d, by = c( "fips" ) ) %>%
  filter( fips %in% d$fips ) %>%
  mutate( 
          E_d.j = ceiling( mean( .$amr, na.rm = TRUE )*pop ), # expected death count for offset term (using JHU data)
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



### (2.0) INLA Components ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# create neighborhood list
nb <- poly2nb( d.jhu )

td <- tempdir() # create temporary directory

# INLA-formatted neighbors list (adjacency matrix)
am <- nb2INLA( paste( td, "inla-mat.adj", sep="/" ), nb )

# create graph object ( I don't like using super assignment operators within function environments but there is no other way to get the g vector to work)
g <<- inla.read.graph( filename = paste( td, "inla-mat.adj", sep="/" ) )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (3.0) County-level Residual Only ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (3.1) Formula ##


# specify formula (we specify the same formula with the same predictors that were selected for in the main analysis)
this.file.jhu.1 <- "/Users/mainovieytesca/Documents/GitHub/COVID-FI-Mortality/04-Tables-Figures/01-main-analysis/jhu-results/model-selection-log.txt" # this obtains the file with the list of predictors that were included in the final model

# clunky way of reading in the predictors .txt file and preparing it for creating the formula
f.1.preds.jhu <- readLines( this.file.jhu.1 ) %>%
  noquote( . ) %>%
  str_extract(., "(?<=R:\\s).*") %>% noquote( . ) %>%
  str_replace(., '\\"', '') %>%
  .[ !is.na(.) ]

f.jhu.1 <- as.formula( paste0( "deaths.jhu.adj ~ I(fi.perc.20/4) +", paste0( f.1.preds.jhu, collapse = '+' ),
                               "+ median.age + sir.jhu + health.index", # basic model variables
                               "+ sir.jhu.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +", # basic + state model variables
                               "+ f( fips, model = 'iid' )" ) ) # random effects


## ---o--- ##


## (3.2) Call INLA ##
this.fit <- NULL

while( is.null( this.fit ) ){ # condition for continuing to loop (iterations until model converges)
  
  this.fit <- 
    tryCatch(
      inla( f.jhu.1, family = "nbinomial", data = data.frame( d.jhu ),
      control.predictor = list( compute = TRUE ), E = E_d.j, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
      control.compute = list( return.marginals.predictor = TRUE, waic = TRUE ), # compute marginals
      control.fixed = list( mean = 0.0, prec = 0.001 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects

error = function( e )
  NULL ) # assign model1 as NULL if `inla` function breaks
}

## ---o--- ##


## (3.3) Extract Residuals ##
res.sp <- this.fit$marginals.random$fips[ 1:nrow( d.jhu ) ]
zeta <- sapply( res.sp, function(x) inla.emarginal( fun = function(x) x, x ) )

## ---o--- ##


## (3.4) Generate Plots ##

# arrange with `layout` and save
lo <- layout( matrix( c(1,1,2,3), nrow = 2,
                      byrow = TRUE ), # matrix of plot arrangements to feed to `layout`
              widths=c(1,1), 
              heights=c(2,2)
)

dev.copy( png,'04-Tables-Figures/08-other-supplementary-files/04a-residuals-analysis-cty-resids.png' )

hist( zeta,
      main = "Fixed Effects + County-Level Residual",
      xlab = "County Residual",
      breaks = 100 )
addfiglab( "A" ) # function for adding letter label to base R plot (see utils.R; this is not MY function)
qqnorm( zeta )

# add LOESS smoother to residuals plot
x <- 1:length( zeta )
line.c <- loess( zeta ~ x )
plot( zeta, ylab = "County Residual",
      main = "Residuals Plot")
lines( predict( line.c ), col = 'red', lwd = 2 )

dev.off()

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (4.0) Add State-level Random Effect ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (4.1) Formula ##
f.jhu.2 <- deaths.jhu.adj ~ I(fi.perc.20/4) +
  f( fips, model = "iid" ) +
  f( state, model = "iid" )

f.jhu.2 <- as.formula( paste0( "deaths.jhu.adj ~ I(fi.perc.20/4) +", paste0( f.1.preds.jhu, collapse = '+' ),
                               "+ median.age + sir.jhu + health.index", # basic model variables
                               "+ sir.jhu.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +", # basic + state model variables
                               "+ f( fips, model = 'iid' ) + f( state, model = 'iid' )" ) )# random effects
## ---o--- ##


## (4.2) Call INLA ##
this.fit.2 <- NULL

while( is.null( this.fit.2 ) ){ # condition for continuing to loop (iterations until model converges)
  
  this.fit.2 <- 
    tryCatch(
      inla( f.jhu.2, family = "nbinomial", data = data.frame( d.jhu ),
            control.predictor = list( compute = TRUE ), E = E_d.j, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
            control.compute = list( return.marginals.predictor = TRUE, waic = TRUE ), # compute marginals
            control.fixed = list( mean = 0.0, prec = 0.001 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
      
      error = function( e )
        NULL ) # assign model1 as NULL if `inla` function breaks
}

## ---o--- ##


## (4.3) Extract Residuals ##
res.sp.s <- this.fit.2$marginals.random$fips[1:nrow( d.jhu )]
zeta.s <- sapply( res.sp.s, function(x) inla.emarginal( fun = function(x) x, x ) )

## ---o--- ##


## (4.4) Generate Plots ##

# arrange with `layout` and save
lo <- layout( matrix( c(1,1,2,3), nrow = 2,
                      byrow = TRUE ), # matrix of plot arrangements to feed to `layout`
              widths=c(1,1), 
              heights=c(2,2)
)

dev.copy( png,'04-Tables-Figures/08-other-supplementary-files/04b-residuals-analysis-cty-st-resids.png' )

hist( zeta.s,
      main = "Fixed Effects + State + County-Level Residuals",
      xlab = "County Residual",
      breaks = 100 )
addfiglab( "B" ) # function for adding letter label to base R plot (see utils.R; this is not MY function)
qqnorm( zeta.s )

# add LOESS smoother to residuals plot
x <- 1:length( zeta.s )
line.c <- loess( zeta.s ~ x )
plot( zeta.s, ylab = "County Residual",
      main = "Residuals Plot")
lines( predict( line.c ), col = 'red', lwd = 2 )

dev.off()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (5.0) Add iCAR Random Effect ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (5.1) Formula ##

# index for county-level iCAR RE
d.jhu$re.s <- 1:nrow( d.jhu )

f.jhu.3 <- deaths.jhu.adj ~ I(fi.perc.20/4) + 
  f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( fips, model = "iid" ) +
  f( state, model = "iid" )

f.jhu.3 <- as.formula( paste0( "deaths.jhu.adj ~ I(fi.perc.20/4) +", paste0( f.1.preds.jhu, collapse = '+' ),
                               "+ median.age + sir.jhu + health.index", # basic model variables
                               "+ sir.jhu.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +", # basic + state model variables
                               "+ f( fips, model = 'iid' ) + f( state, model = 'iid' ) + f( re.s, model = 'besag', graph = g, scale.model = T )" ) )# random effects

## ---o--- ##


## (5.2) Call INLA ##
this.fit.3 <- NULL

while( is.null( this.fit.3 ) ){ # condition for continuing to loop (iterations until model converges)
  
  this.fit.3 <- 
    tryCatch(
      inla( f.jhu.3, family = "nbinomial", data = data.frame( d.jhu ),
            control.predictor = list( compute = TRUE ), E = E_d.j, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
            control.compute = list( return.marginals.predictor = TRUE, waic = TRUE ), # compute marginals
            control.fixed = list( mean = 0.0, prec = 0.001 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
      
      error = function( e )
        NULL ) # assign model1 as NULL if `inla` function breaks
}

## ---o--- ##


## (5.3) Extract Residuals ##
res.sp.sg <- this.fit.3$marginals.random$fips[1:nrow( d.jhu )]
zeta.sg <- sapply( res.sp.sg, function(x) inla.emarginal( fun = function(x) x, x ) )

## ---o--- ##


## (5.4) Generate Plots ##

# arrange with `layout` and save
lo <- layout( matrix( c(1,1,2,3), nrow = 2,
                      byrow = TRUE ), # matrix of plot arrangements to feed to `layout`
              widths=c(1,1), 
              heights=c(2,2)
)

dev.copy( png,'04-Tables-Figures/08-other-supplementary-files/04c-residuals-analysis-cty-st-icar-resids.png' )

hist( zeta.sg,
      main = "Fixed Effects + iCAR + State\n+ County-Level Residuals",
      xlab = "County Residual",
      breaks = 100 )
addfiglab( "C" ) # function for adding letter label to base R plot (see utils.R; this is not MY function)
qqnorm( zeta.sg )

# add LOESS smoother to residuals plot
x <- 1:length( zeta.sg )
line.c <- loess( zeta.sg ~ x )
plot( zeta.sg, ylab = "County Residual",
      main = "Residuals Plot")
lines( predict( line.c ), col = 'red', lwd = 2 )

dev.off()

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (6.0) Evaluate Changes to Fixed Effects Coefficients ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (6.1) County-level Residual Model ##

this.fit.4 <- NULL

while( is.null( this.fit.4 ) ){ # condition for continuing to loop (iterations until model converges)
  
  this.fit.4 <- 
    tryCatch(
      inla( f.jhu.1, family = "nbinomial", data = data.frame( d.jhu ),
            control.predictor = list( compute = TRUE ), E = E_d.j, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
            control.compute = list( return.marginals.predictor = TRUE, waic = TRUE ), # compute marginals
            control.fixed = list( mean = 0.0, prec = 0.001 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
      
      error = function( e )
        NULL ) # assign model1 as NULL if `inla` function breaks
}

fe.4 <- this.fit.4$summary.fixed

# generate table of results
res.4 <- data.frame( 
  coef = exp( fe.4[ str_detect( rownames(fe.4), "fi.perc"), "mean"] ), # food insecurity coefficient,
  low = exp( fe.4[ str_detect( rownames(fe.4), "fi.perc"), "0.025quant"] ), # 95% CrI lower bound
  upper = exp( fe.4[ str_detect( rownames(fe.4), "fi.perc"), "0.975quant"] ), # 95% CrI upper bound
  waic = this.fit.4$waic$waic # waic
) %>% round( 3 )

## ---o--- ##


## (6.2) Add State-level Residual ##

this.fit.5 <- NULL

while( is.null( this.fit.5 ) ){ # condition for continuing to loop (iterations until model converges)
  
  this.fit.5 <- 
    tryCatch(
      inla( f.jhu.2, family = "nbinomial", data = data.frame( d.jhu ),
            control.predictor = list( compute = TRUE ), E = E_d.j, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
            control.compute = list( return.marginals.predictor = TRUE, waic = TRUE ), # compute marginals
            control.fixed = list( mean = 0.0, prec = 0.001 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
      
      error = function( e )
        NULL ) # assign model1 as NULL if `inla` function breaks
}

fe.5 <- this.fit.5$summary.fixed

# generate table of results
res.5 <- data.frame( 
  coef = exp( fe.5[ str_detect( rownames(fe.5), "fi.perc"), "mean"] ), # food insecurity coefficient,
  low = exp( fe.5[ str_detect( rownames(fe.5), "fi.perc"), "0.025quant"] ), # 95% CrI lower bound
  upper = exp( fe.5[ str_detect( rownames(fe.5), "fi.perc"), "0.975quant"] ), # 95% CrI upper bound
  waic = this.fit.5$waic$waic # waic
) %>% round( 3 )

## ---o--- ##


## (6.3) Add iCAR Term ##

this.fit.6 <- NULL

while( is.null( this.fit.6 ) ){ # condition for continuing to loop (iterations until model converges)
  
  this.fit.6 <- 
    tryCatch(
      inla( f.jhu.3, family = "nbinomial", data = data.frame( d.jhu ),
            control.predictor = list( compute = TRUE ), E = E_d.j, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
            control.compute = list( return.marginals.predictor = TRUE, waic = TRUE ), # compute marginals
            control.fixed = list( mean = 0.0, prec = 0.001 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
      
      error = function( e )
        NULL ) # assign model1 as NULL if `inla` function breaks
}

fe.6 <- this.fit.6$summary.fixed

# generate table of results
res.6 <- data.frame( 
  coef = exp( fe.6[ str_detect( rownames(fe.6), "fi.perc"), "mean"] ), # food insecurity coefficient,
  low = exp( fe.6[ str_detect( rownames(fe.6), "fi.perc"), "0.025quant"] ), # 95% CrI lower bound
  upper = exp( fe.6[ str_detect( rownames(fe.6), "fi.perc"), "0.975quant"] ), # 95% CrI upper bound
  waic = this.fit.6$waic$waic # waic
) %>% round( 3 )

## ---o--- ##


## (6.4) Bind All Results Together and Save ##

rbind( res.4, res.5, res.6 ) %>%
  write.table( "04-Tables-Figures/08-other-supplementary-files/05-rand-eff-model-analysis.txt",
               sep = "," )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------


