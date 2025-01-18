###-------------------------------------------------------------
###   SENSITIVITY-03: DIFFERENT CHOICES OF PRIORS
###------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we fit the conditional autoregressive models in INLA (with the same specification as in
# "R/04-primary-analysis-inla.R") but we implement different specifications of the priors for the hyper-
# parameters of the model (i.e., the log precisions of the random effects). We only perform this analysis
# using the JHU data to get a sense of how different priors are changing the the results. If the results 
# change drastically, this indicates that the priors are having a greater influence on the posteriors
# than the likelihood. If the results do not change drastically, it indicates the "data are speaking" and
# the posteriors are largely being influenced by the likelihood as opposed to priors.
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

# specify formula (we specify the same formula with the same predictors that were selected for in the main analysis)
this.file.jhu.1 <- "/Users/mainovieytesca/Documents/GitHub/COVID-FI-Mortality/04-Tables-Figures/01-main-analysis/jhu-results/model-selection-log.txt" # this obtains the file with the list of predictors that were included in the final model

# clunky way of reading in the predictors .txt file and preparing it for creating the formula
f.1.preds.jhu <- readLines( this.file.jhu.1 ) %>%
  noquote( . ) %>%
  str_extract(., "(?<=R:\\s).*") %>% noquote( . ) %>%
  str_replace(., '\\"', '') %>%
  .[ !is.na(.) ]

f.jhu.a <- as.formula( paste0( "deaths.jhu ~ I(fi.perc.20/4) +", paste0( f.1.preds.jhu, collapse = '+' ),
                               "+ median.age + sir.jhu + health.index", # basic model variables
                               "+ sir.jhu.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +", # basic + state model variables
                               "+ f( re.s, model = 'besag', graph = g, scale.model = T ) + f( state, model = 'iid') + f( fips, model = 'iid' )" ) ) # random effects

## ---o--- ##


## (2.2) Models with Random Effects and Offset Terms Only ("Null/baseline Model") ##

f.jhu.model.null <- deaths.jhu.adj ~ f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

## ---o--- ##

## (2.3) Models with Random Effects, Offset Terms, and Fixed Effects for Health Index and Median Age ("Basic Model") ##
f.jhu.model.2 <- deaths.jhu.adj ~ health.index + median.age + sir.jhu + f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

## (2.4) Models with Random Effects, Offset Terms, and Fixed Effects for Health Index and Median Age
## and State-Level Variables ("Basic + State Model") ##
f.jhu.model.3 <- deaths.jhu.adj ~ health.index + median.age + sir.jhu +
  sir.jhu.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +
  f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (3.0) Set the Parameters of the Priors to Loop On ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# 4 log-gamma specifications
priors.list <- list( 
  lg.prior.1 = list( prec = list( prior = "loggamma", param = c(1, 0.005), initial = 4 ) ),
  lg.prior.2 = list( prec = list( prior = "loggamma", param = c(0.5, 0.01), initial = 4 ) ),
  lg.prior.3 = list( prec = list( prior = "loggamma", param = c(0.25, 1), initial = 4 ) ),
  lg.prior.4 = list( prec = list( prior = "loggamma", param = c(0.5, 1), initial = 4 ) )
)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (4.0) INLA Analysis: National (Nearly All U.S. Counties) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (4.1) Call INLA ##

# states to remove from maps

states.fips.out <- unique( c( "69", "60", "66", "78", "15", "02", "72" ) )

# initialize list to store results
res.priors.jhu <- list()

# loop through each specification of the priors 
for( i in seq_along( priors.list ) ){
  
  res.priors.jhu[[i]] <- car_inla_analysis( d.jhu = d.poly, d.cdc = d.poly,
                                formula.jhu = f.jhu.a, 
                                which.model = "jhu", E = "E_d.j", car.prior = priors.list[[i]],
                                ur.prior = priors.list[[i]],
                                term = "fi", null.model.formula = f.jhu.model.null,
                                model.2.formula = f.jhu.model.2,
                                model.3.formula = f.jhu.model.3,
                                selection.criterion = NULL,
                                state.codes.omit = states.fips.out )

}

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (5.0) Generate Combined Plot of Fractional Variance Explained by Random Effects ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# save spatial fractional variance plots arranged
ggarrange(  res.priors.jhu[[1]]$rand.eff.var.plot +
             ggtitle( TeX( "log-$\\gamma(1, 0.005)$" ) ) +
             theme( plot.title = element_text( color = "grey44", hjust = 1 ) ),
           res.priors.jhu[[2]]$rand.eff.var.plot +
             ggtitle( TeX( "log-$\\gamma(0.5, 0.01)$" ) ) +
             theme( plot.title = element_text( color = "grey44", hjust = 1 ) ),
           res.priors.jhu[[3]]$rand.eff.var.plot +
             ggtitle( TeX( "log-$\\gamma(0.25, 1)$" ) ) +
             theme( plot.title = element_text( color = "grey44", hjust = 1 ) ),
           res.priors.jhu[[4]]$rand.eff.var.plot +
             ggtitle( TeX( "log-$\\gamma(0.5, 1)$" ) ) +
             theme( plot.title = element_text( color = "grey44", hjust = 1 ) ),
           nrow = 3, ncol = 2,
           common.legend = TRUE, legend = "bottom" )

# dir.create( "04-Tables-Figures/06-sensitivity-analyses/04-choice-of-priors/" )
# dir.create( "04-Tables-Figures/06-sensitivity-analyses/04-choice-of-priors/jhu-results" )

ggsave( "04-Tables-Figures/06-sensitivity-analyses/04-choice-of-priors/jhu-results/priors-fractional-spatial-variance-plot.png",
        width = 6, height = 8 )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (6.0) Generate Table of Fixed Effects to Show How they Change with Different Priors###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


# fixed effects
bind_cols( data.frame( prior.choice = c( "log-gamma(1, 0.005)", "log-gamma(0.5, 0.01)",
                                         "log-gamma(0.25, 1)","log-gamma(0.5, 1)") ),
bind_rows( res.priors.jhu[[1]]$fixed.term,
           res.priors.jhu[[2]]$fixed.term,
           res.priors.jhu[[3]]$fixed.term,
           res.priors.jhu[[4]]$fixed.term)
) %>%
  write.table( ., "04-Tables-Figures/06-sensitivity-analyses/04-choice-of-priors/jhu-results/priors-fi-fixed-effects-estimates.txt",
               sep = "," )
           
           
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

