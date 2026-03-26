###------------------------------------------------------------
###   02-COVID-19 MORTALITY DATA IMPORT, WRANGLING, AND MERGE
###------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we will import mortality data from two sources:
# i. the restricted access files from the Centers from Disease Control and Prevention (CDC) and
# ii. COVID-19 Data Repository by the Center for Systems Science and Engineering (CSSE) 
# at Johns Hopkins University (JHU).
#
# For the CDC data, we will aggregate deaths counts up to the county level and compute
# age-standardized mortality counts and rates using the direct standardization method. For the
# JHU data, we will obtain the raw mortality counts and compute crude mortality rates for 
# each county. We also import vaccination data and state testing data
# 
#
# INPUT DATA FILES: 
# i. CDC Case Surveillance Restricted Access Detailed Data
# ii. JHU COVID-19 time series data containing cases and mortality counts
# iii. CDC COVID-19 vaccination data
# iv. JHU COVID-19 state testing rate data
# v. "02-Data-Wrangled/01-covariate-merge.rds"
#
#
# OUTPUT DATA FILE: "03-Data-Rodeo/01-analytic-data.rds"
#
# A Special Note: The restricted access files from the CDC are not publicly available. 
# They can be requested using the link provided below. 
#
#
# Resources (Accessed 10 November 2022): 
# i. CDC Data Request Link: https://data.cdc.gov/Case-Surveillance/COVID-19-Case-Surveillance-Restricted-Access-Detai/mbd7-r32t
# ii. JHU Data Repository: https://github.com/CSSEGISandData/COVID-19
# iii. Method for age-standardization of mortality rates: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3406211/
# iv. MMWR week calendar: (https://health.maryland.gov/phpa/OIDEOR/CIDSOR/NEDSS/MMWR_Calendar.pdf) 
# v. Pub on lessening health disparities as time wore on in pandemic: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9098236/
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

library( tidyverse )
library( lubridate )  # for date/time wrangling
library( tidycensus ) # for Census API queries
library( RCurl )      # for URL fetching

# read in helper functions
source( "R/utils.R")

cov.data <- readRDS( "02-Data-Wrangled/01-covariate-merge.rds")


### (1.0) Read-in CDC Restricted Access COVID Data and Wrangle ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# **NOTE: These data ARE NOT publicly available and are being read in from a personal hard drive.

## Use for-loop to read in .csv files ##

covid.data.list <- list( )
for( i in 1:10 ){
  
  start.time <- Sys.time( ) # track start time
  
  # restricted access CDC data is on my personal hard drive and not in the repo
  covid.data.list[[i]] <- read.csv( paste0( '/Volumes/T7 Shield/Arthur-Lab/CDC-COVID-FI/CDC Restricted Access Data/04_04_22/COVID_Cases_Restricted_Detailed_04042022_Part_', i, '.csv' ) )[ , c( 'county_fips_code', 'death_yn', 'cdc_report_dt','cdc_case_earliest_dt', 'age_group' ) ]
  
  end.time <- Sys.time( )  # track end time
  print( end.time - start.time )  # print time elapsed
  print( paste0( "Iteration ", i, "/10 complete." ) )  # track loop iteration
}

# bind into single dataset
d.raw <- do.call( "rbind", covid.data.list )


# Assume deaths come, on average, 8 days after onset since we are not given date of expiration
# source: https://pubmed.ncbi.nlm.nih.gov/33709641/
d.raw$death_dt <- as.Date( d.raw$cdc_case_earliest_dt ) + 8
d.raw[ d.raw$death_yn != 'Yes', 'death_dt' ] <- NA # only those with death = 'yes' keep death dates

# fix fips codes so they are standard 5 digits
d.raw <- d.raw %>%
  mutate( county_fips_code = ifelse( as.numeric( county_fips_code ) < 10000, paste0( "0", county_fips_code ), county_fips_code ) )


# Maximum date in this dataset
max( d.raw$death_dt, na.rm = T ) # Max is 12-26-2021

min( d.raw$death_dt, na.rm = T )  # Min is 01/22/2020 ( cases )

# prep data for entire study period (1st year of pandemic)
( d.2 <- death_data_prep( df = d.raw, mort.date.start = "2020-03-25",
                        mort.date.end = "2021-12-25",
                        case.date.start = "2020-03-17",
                        case.date.end = "2021-12-17" ) ) %>%
    saveRDS( "02-Data-Wrangled/02-cdc-raw-mortality-aggregated.rds" )

# read in working dataset so that above code does not need to be re-run in case script is closed
d.2 <- readRDS( "02-Data-Wrangled/02-cdc-raw-mortality-aggregated.rds" )


# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (2.0) Compute Crude Mortality Rates Across Age Strata Within Counties ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# bring in county data on age structure from `tidycensus` to compute
# COVID-19 mortality rates within each strata with each county
cnty.ref <- get_acs( year = 2020,
                     variables = c( AGE04.TOT = "S0101_C01_002",
                                    AGE59.TOT = "S0101_C01_003",
                                    AGE1014.TOT = "S0101_C01_004",
                                    AGE1519.TOT = "S0101_C01_005",
                                    AGE2024.TOT = "S0101_C01_006",
                                    AGE2529.TOT = "S0101_C01_007",
                                    AGE3034.TOT = "S0101_C01_008",
                                    AGE3539.TOT = "S0101_C01_009",
                                    AGE4044.TOT = "S0101_C01_010",
                                    AGE4549.TOT = "S0101_C01_011",
                                    AGE5054.TOT = "S0101_C01_012",
                                    AGE5559.TOT = "S0101_C01_013",
                                    AGE6064.TOT = "S0101_C01_014",
                                    AGE6569.TOT = "S0101_C01_015",
                                    AGE7074.TOT = "S0101_C01_016",
                                    AGE7579.TOT = "S0101_C01_017",
                                    AGE8084.TOT = "S0101_C01_018",
                                    AGE85PLUS.TOT = "S0101_C01_019" ), 
                     geography = "county",
                     cache = T ) %>%
  
  # pivot to wide format
  pivot_wider( id_cols = GEOID, names_from = variable, values_from = estimate ) %>%
  rename( fips = GEOID ) %>% # rename identifier column
  data.frame()               # convert to df


# reclassify age columns
d.3 <- left_join( d.2, cnty.ref, by = 'fips' ) %>%
  mutate( `0-9 Years` = AGE04.TOT + AGE59.TOT,
          `10-19 Years` = AGE1014.TOT + AGE1519.TOT,
          `20-29 Years` = AGE2024.TOT + AGE2529.TOT,
          `30-39 Years` = AGE3034.TOT + AGE3539.TOT,
          `40-49 Years` = AGE4044.TOT + AGE4549.TOT,
          `50-59 Years` = AGE5054.TOT + AGE5559.TOT,
          `60-69 Years` = AGE6064.TOT + AGE6569.TOT,
          `70-79 Years` = AGE7074.TOT + AGE7579.TOT,
          `80+ Years` = AGE8084.TOT + AGE85PLUS.TOT ) 


levs.agegrp <- levels( factor( d.3$age_group ) )


# compute age-group-specific crude mortality rates and also keep the raw county mortality counts as well
d.raw.cts <- d.3 %>%
  group_by( fips ) %>%
  mutate( raw.deaths.cdc = sum( cum.deaths.cdc ) ) %>%
  ungroup() %>%
  distinct( fips, raw.deaths.cdc )

d.4 <- d.3 %>%
  mutate( grp.age.pop = ifelse( age_group == levs.agegrp[ 1 ], `0-9 Years`,
                                ifelse( age_group == levs.agegrp[ 2 ], `10-19 Years`,
                                        ifelse( age_group == levs.agegrp[ 3 ], `20-29 Years`,
                                                ifelse( age_group == levs.agegrp[ 4 ], `30-39 Years`,
                                                        ifelse( age_group == levs.agegrp[ 5 ], `40-49 Years`,
                                                                ifelse( age_group == levs.agegrp[ 6 ], `50-59 Years`,
                                                                        ifelse( age_group == levs.agegrp[ 7 ], `60-69 Years`,
                                                                                ifelse( age_group == levs.agegrp[ 8 ], `70-79 Years`,
                                                                                        ifelse( age_group == levs.agegrp[ 9 ], `80+ Years`, NA ) ) ) ) ) ) ) ) ) ) %>%
  mutate( crude.mort.rate = ( cum.deaths.cdc / grp.age.pop ) ) %>%
  left_join( ., d.raw.cts, by = "fips" )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (3.0) Direct Standardization: Age-Standardize the Crude Mortality Rates ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# method for direct standardization is described here: 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3406211/

# bring in reference population data (U.S. age structure)
ref.age <- get_acs( year = 2020,
                    table = "S0101", 
                    geography = "us",
                    cache = T ) %>% 
  pivot_wider( id_cols = GEOID, names_from = variable, values_from = estimate ) %>%
  rename( AGE04.TOT = S0101_C01_002,
   AGE59.TOT = S0101_C01_003,
   AGE1014.TOT = S0101_C01_004,
   AGE1519.TOT = S0101_C01_005,
   AGE2024.TOT = S0101_C01_006,
   AGE2529.TOT = S0101_C01_007,
   AGE3034.TOT = S0101_C01_008,
   AGE3539.TOT = S0101_C01_009,
   AGE4044.TOT = S0101_C01_010,
   AGE4549.TOT = S0101_C01_011,
   AGE5054.TOT = S0101_C01_012,
   AGE5559.TOT = S0101_C01_013,
   AGE6064.TOT = S0101_C01_014,
   AGE6569.TOT = S0101_C01_015,
   AGE7074.TOT = S0101_C01_016,
   AGE7579.TOT = S0101_C01_017,
   AGE8084.TOT = S0101_C01_018,
   AGE85PLUS.TOT = S0101_C01_019 ) %>%
  mutate( `0-9 Years` = AGE04.TOT + AGE59.TOT,
   `10-19 Years` = AGE1014.TOT + AGE1519.TOT,
   `20-29 Years` = AGE2024.TOT + AGE2529.TOT,
   `30-39 Years` = AGE3034.TOT + AGE3539.TOT,
   `40-49 Years` = AGE4044.TOT + AGE4549.TOT,
   `50-59 Years` = AGE5054.TOT + AGE5559.TOT,
   `60-69 Years` = AGE6064.TOT + AGE6569.TOT,
   `70-79 Years` = AGE7074.TOT + AGE7579.TOT,
   `80+ Years` = AGE8084.TOT + AGE85PLUS.TOT ) %>%
  select( `0-9 Years`, `10-19 Years`, `20-29 Years`, `30-39 Years`,
          `40-49 Years`, `50-59 Years`, `60-69 Years`, `70-79 Years`,
          `80+ Years` ) %>%
  pivot_longer( cols = contains( "Years" ),
                names_to = "age_group", 
                values_to = "ref.pop" ) # put into long format for easy merge and to facilitate computations 




## Multiply county crude mortality rates by age group population in reference population ##

# fix levels of age_group variable to ensure smooth join
d.4 %>%
  mutate( age_group = as.factor( age_group ) )

# recode levels
d.4$age_group <- recode_factor( d.4$age_group, "0 - 9 Years" = "0-9 Years",
                                "10 - 19 Years" = "10-19 Years",
                                "20 - 29 Years" = "20-29 Years",
                                "30 - 39 Years" = "30-39 Years",
                                "40 - 49 Years" = "40-49 Years",
                                "50 - 59 Years" = "50-59 Years",
                                "60 - 69 Years" = "60-69 Years",
                                "70 - 79 Years" = "70-79 Years" )


# join reference population data to working data and conduct computations 
d.5 <- left_join( d.4, ref.age, by = "age_group" ) %>%
  
  # compute expected deaths in reference population
  mutate( exp.deaths.cdc = ceiling( crude.mort.rate*ref.pop ) )%>% 
  group_by( fips ) %>%
  
  # sum expected deaths within county fips and divide by reference population size
  # Note: mort rates are multiplied by 100,000 to represent deaths per 100,000 in the population
  mutate( age.adj.mort.rate = 100000*sum( exp.deaths.cdc ) / sum( ref.age$ref.pop ) ) %>%
  ungroup() %>%
  distinct( fips, age.adj.mort.rate, cases.cdc, raw.deaths.cdc ) %>% # keep distinct rows with adj mortality rates, count, and fips
  arrange(fips)
  



## Join mortality rates back to working data ##

d.6 <- left_join( cov.data, d.5, by = "fips") %>%
  mutate( deaths.cdc = ceiling( ( age.adj.mort.rate / 100000 )*pop ) ) # compute age-adjusted death count

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (4.0) Johns Hopkins University COVID-19 Mortality Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## JHU data download (accessed 8 November 2022)
# https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data

# observed data from JHU CRC 
jhu.mort <- read.csv( text = getURL( "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv" ) ) %>%
  rename( fips = FIPS ) %>%
  mutate( deaths.jhu = as.numeric( X12.25.21 - X3.25.20 ),
          fips = ifelse( as.numeric( fips ) < 10000, paste0( "0", fips ), fips ) ) %>%
  select( fips, deaths.jhu )

# get national reference COVID mortality data
mort.age.us <- read.csv( text = getURL('https://data.cdc.gov/api/views/vsak-wrfu/rows.csv?accessType=DOWNLOAD') ) %>%
  mutate( End.Week = format( End.Week ,format = '%m/%d/%Y' ) ) %>%
  filter( End.Week >= format( '03/25/2020', format = '%m/%d/%Y' ) &
            End.Week <= format( '12/25/2021', format = '%m/%d/%Y' ) & 
            Age.Group != "All Ages" & Sex == "All Sex") %>%
  group_by( Age.Group ) %>%
  mutate( covid.deaths.cum = sum( COVID.19.Deaths ) ) %>%
  ungroup() %>%
  distinct( state = State, age.group = Age.Group, covid.deaths.cum )

# combine ages under 1 years and 1-4 years into single group since covid data combine these groups
covid.mort.age <- rbind( data.frame( state = "United States",
                                     age.group = "Under 5 Years",
                                     covid.deaths.cum = sum( mort.age.us[ which( mort.age.us$age.group %in% c( "Under 1 year", "1-4 Years" ) ), "covid.deaths.cum"]
                                     ) ),
                         mort.age.us[ which( mort.age.us$age.group %notin% c( "Under 1 year", "1-4 Years" ) ), ]
)


us.ref.21 <- get_acs( year = 2021,
                      variables = c( AGE04.TOT = "S0101_C01_002",
                                     AGE0509.TOT = "S0101_C01_003",
                                     AGE1014.TOT = "S0101_C01_004",
                                     AGE1519.TOT = "S0101_C01_005",
                                     AGE2024.TOT = "S0101_C01_006",
                                     AGE2529.TOT = "S0101_C01_007",
                                     AGE3034.TOT = "S0101_C01_008",
                                     AGE3539.TOT = "S0101_C01_009",
                                     AGE4044.TOT = "S0101_C01_010",
                                     AGE4549.TOT = "S0101_C01_011",
                                     AGE5054.TOT = "S0101_C01_012",
                                     AGE5559.TOT = "S0101_C01_013",
                                     AGE6064.TOT = "S0101_C01_014",
                                     AGE6569.TOT = "S0101_C01_015",
                                     AGE7074.TOT = "S0101_C01_016",
                                     AGE7579.TOT = "S0101_C01_017",
                                     AGE8084.TOT = "S0101_C01_018",
                                     AGE85PLUS.TOT = "S0101_C01_019" ), 
                      geography = "us",
                      cache = T ) %>%
  
  # pivot to wide format
  pivot_wider( id_cols = GEOID, names_from = variable, values_from = estimate ) %>%
  rename( fips = GEOID ) %>% # rename identifier column
  
  # re-categorize ages
  mutate( 
    AGE0514.TOT = AGE0509.TOT + AGE1014.TOT,
    AGE1524.TOT = AGE1519.TOT + AGE2024.TOT,
    AGE2534.TOT = AGE2529.TOT + AGE3034.TOT,
    AGE3544.TOT = AGE3539.TOT + AGE4044.TOT,
    AGE4554.TOT = AGE4549.TOT + AGE5054.TOT,
    AGE5564.TOT = AGE5559.TOT + AGE6064.TOT,
    AGE6574.TOT = AGE6569.TOT + AGE7074.TOT,
    AGE7584.TOT = AGE7579.TOT + AGE8084.TOT )  %>%
  select( AGE04.TOT, AGE0514.TOT, AGE1524.TOT, AGE2534.TOT, AGE3544.TOT,
          AGE4554.TOT, AGE5564.TOT, AGE6574.TOT, AGE7584.TOT, AGE85PLUS.TOT ) %>%
  unlist() # convert to vector

covid.standard.ref <- covid.mort.age %>%
  mutate( pop.size = us.ref.21, # add population size of age group to dataset
          covid.mort.rate = covid.deaths.cum / pop.size ) # compute mortality rate in reference population

wide.covid.ref <- covid.standard.ref %>%
  select( age.group, covid.mort.rate ) %>%
  pivot_wider( names_from = age.group, 
               values_from = covid.mort.rate,
               names_sep = "." )

## now get county level data to compute expected deaths and age adjusted death count

d.7 <- get_acs( year = 2021,
                     variables = c( AGE04.TOT = "S0101_C01_002",
                                    AGE0509.TOT = "S0101_C01_003",
                                    AGE1014.TOT = "S0101_C01_004",
                                    AGE1519.TOT = "S0101_C01_005",
                                    AGE2024.TOT = "S0101_C01_006",
                                    AGE2529.TOT = "S0101_C01_007",
                                    AGE3034.TOT = "S0101_C01_008",
                                    AGE3539.TOT = "S0101_C01_009",
                                    AGE4044.TOT = "S0101_C01_010",
                                    AGE4549.TOT = "S0101_C01_011",
                                    AGE5054.TOT = "S0101_C01_012",
                                    AGE5559.TOT = "S0101_C01_013",
                                    AGE6064.TOT = "S0101_C01_014",
                                    AGE6569.TOT = "S0101_C01_015",
                                    AGE7074.TOT = "S0101_C01_016",
                                    AGE7579.TOT = "S0101_C01_017",
                                    AGE8084.TOT = "S0101_C01_018",
                                    AGE85PLUS.TOT = "S0101_C01_019" ), 
                     geography = "county",
                     cache = T ) %>%
  
  # pivot to wide format
  pivot_wider( id_cols = GEOID, names_from = variable, values_from = estimate ) %>%
  rename( fips = GEOID ) %>% # rename identifier column
  # re-categorize ages
  mutate( 
    AGE0514.TOT = AGE0509.TOT + AGE1014.TOT,
    AGE1524.TOT = AGE1519.TOT + AGE2024.TOT,
    AGE2534.TOT = AGE2529.TOT + AGE3034.TOT,
    AGE3544.TOT = AGE3539.TOT + AGE4044.TOT,
    AGE4554.TOT = AGE4549.TOT + AGE5054.TOT,
    AGE5564.TOT = AGE5559.TOT + AGE6064.TOT,
    AGE6574.TOT = AGE6569.TOT + AGE7074.TOT,
    AGE7584.TOT = AGE7579.TOT + AGE8084.TOT )  %>%
  mutate( exp.deaths.05under =  AGE04.TOT*wide.covid.ref$`Under 5 Years` ,
          exp.deaths.0514 =  AGE0514.TOT*wide.covid.ref$`5-14 Years` ,
          exp.deaths.1524 =  AGE1524.TOT*wide.covid.ref$`15-24 Years` ,
          exp.deaths.2534 =  AGE2534.TOT*wide.covid.ref$`25-34 Years` ,
          exp.deaths.3544 =  AGE3544.TOT*wide.covid.ref$`35-44 Years` ,
          exp.deaths.4554 =  AGE4554.TOT*wide.covid.ref$`45-54 Years` ,
          exp.deaths.5564 =  AGE5564.TOT*wide.covid.ref$`55-64 Years` ,
          exp.deaths.6574 =  AGE6574.TOT*wide.covid.ref$`65-74 Years` ,
          exp.deaths.7584 =  AGE7584.TOT*wide.covid.ref$`75-84 Years` ,
          exp.deaths.over85 =  AGE85PLUS.TOT*wide.covid.ref$`85 Years and Over` ,
          
          # compute total expected deaths for a county
          exp.deaths = exp.deaths.05under + exp.deaths.0514 + exp.deaths.1524 +
            exp.deaths.2534 + exp.deaths.3544 + exp.deaths.4554 + exp.deaths.5564 +
            exp.deaths.6574 + exp.deaths.7584 + exp.deaths.over85 ) %>%
  
  select( fips, exp.deaths ) %>%
  left_join( jhu.mort, ., by = "fips" ) %>%
  left_join( d.6, ., by = "fips" ) %>%
  
  # final computations: smr, amr, and age standardized death count
  mutate( exp.deaths = ifelse( is.na( exp.deaths ), 0, exp.deaths ),
          smr = deaths.jhu / exp.deaths,
          amr = smr*( deaths.jhu / pop ),
          deaths.jhu.adj = ceiling( pop*amr) )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (5.0) JHU County-Level Case Time Series Data for Computing Incidence ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# download data (Accessed 10 November 2022)
jhu.case <- read.csv( text = getURL( "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv" ) )

## Average Lag between case and death due to COVID 19: 8.053 days (will round to 8) ##
# source: https://pubmed.ncbi.nlm.nih.gov/33709641/

d.8 <- jhu.case %>% 
  rename( fips = FIPS ) %>%     # select dates based on 8-day lag
  mutate( cases.jhu = as.numeric( X12.17.21 - X3.17.20 ) ) %>%
  select( fips, cases.jhu ) %>%
  # make fips column character and standardized to 5-digits
  mutate( fips = ifelse( as.numeric( fips ) < 10000, paste0( "0", fips ), fips ) ) %>%
  
  left_join( d.7, ., by = "fips" ) %>%# join to working data
  
  # compute county COVID-19 cumulative incidence proportions using both CDC and JHU data during study period.
  # units are cases / 100k in the population (also compute state means and countsas a state-level variable)
  mutate( 
    inc.prop.jhu = 100000*( cases.jhu / pop ),
    
    inc.prop.cdc = 100000*( cases.cdc / pop ) ) %>%
  group_by( state ) %>%
  mutate( inc.prop.jhu.state = weighted.mean( x = inc.prop.jhu,
                                              w = pop ),
          inc.prop.cdc.state = weighted.mean( x = inc.prop.cdc,
                                              w = pop ),
          cases.jhu.state = sum( cases.jhu, na.rm = TRUE ),
          cases.cdc.state = sum( cases.cdc, na.rm = TRUE ),
          state.pop = sum( pop, na.rm = TRUE ) ) %>% # also add a state-level population variable
  ungroup()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (6.0) CDC County Vaccination Rate Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# download data (Accessed 11 Nov 2022)
jhu.vax <- read.csv( text = getURL( "https://data.cdc.gov/api/views/8xkx-amqh/rows.csv?accessType=DOWNLOAD" ) )

# wrangle and merge
( d.9 <- left_join( d.8, ( jhu.vax %>%
                        mutate( Date = mdy( Date ) ) %>% # format date column
                        
                        # we will consider those vaccinated as of the end of the analytic period (03-25-2021)
                        filter( Date == as_date( "2021-12-25" ) ) %>%  
                        select( fips = FIPS, perc.vaccinated = Series_Complete_Pop_Pct ) %>% 
  filter( ! fips == "UNK" ) ) ) %>%
    group_by( state ) %>% # compute state mean vaccination rate as a state-level variable
    mutate( perc.vaccinated.state = weighted.mean( x = perc.vaccinated,
                                                w = pop ) ) %>%
  ungroup() )
    
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (7.0) Exclusions ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (7.1) Exclude Counties with Missing CDC Mortality Data ##

ex.1 <- nrow( d.9 )

( ( step.1 <- d.9 %>%
  
  # create column for state fips code
  mutate( state.code = str_extract( fips, "^.{2}"),
          all = 1 ) %>%
  
  # some exclusions
  filter( !is.na( deaths.cdc ) & !is.na( deaths.jhu ) ) ) %>% 
    summarise( counties.remaining = n( ) , counties.excluded = ex.1 - n( ) ) -> step.1.rem ) 

# counties.remaining counties.excluded
#               3219                 2

# see who got dropped
d.9[ which( d.9$fips %notin% step.1$fips ), c( "fips", "state", "county" ) ]

## ---o--- ##


## (7.2) Exclude Counties from PR ##

( ( step.2 <- step.1 %>%
      # Puerto Rico omitted
  filter( !state.code %in% c( "72" )  ) ) %>%
    summarise( counties.remaining = n( ) , counties.excluded = step.1.rem$counties.remaining - n( ) ) -> step.2.rem ) 
    
# counties.remaining counties.excluded
#               3141               78

## ---o--- ##


## (7.2) Exclude Counties with Incomplete Covariate Data ##

these.covs <- c( "pop", "pop.density", "unemp.rate", "fi.perc.20", "mort.rate.jhu", "age.adj.mort.rate", "median.age","perc.female", 
                    "perc.native","perc.hisp","perc.black","perc.asian", "urb.cat.code", "no.health.insur", "disability",
                    "perc.nh.white", "perc.fb", "pct.emp.trade","pct.emp.trans","ed.1less.than.hspct","no.vehic",
                    "ed.5college.plus.pct","poverty.rate","cases.jhu", "avg.hhsize", "gini.index", "elec.2020.margin", 
                    "inc.prop.jhu","perc.vaccinated","ratio.pop.edp","health.index" )

( step.3 <- step.2 %>%
      drop_na( any_of( these.covs ) ) ) %>%
    summarise( counties.remaining = n( ) , 
               counties.excluded = step.2.rem$counties.remaining - n( ) ) 

# counties.remaining counties.excluded
#               3101                 5

# see who got dropped
( cov.drop <- step.2[ which( step.2$fips %notin% step.3$fips ), c( "fips", "state", "county" ) ] ) # NOTE: all of Alaska gets omitted give missingness in some data

# ensure there are no county duplicates
sum( duplicated( step.3$fips ) )
# no dups

## ---o--- ##

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### Save Final Dataset ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS( step.3, "03-Data-Rodeo/01-analytic-data.rds")
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

