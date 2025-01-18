###------------------------------------------------------------
###   05-LONG-FORM DATA WRANGLE FOR WEEKLY TIME SERIES
###------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we create a long-form dataset by pivoting on `fips` and
# `week`. The purpose of this is to generate a dataset that will allow us to 
# evaluate weekly changes in mortality counts for both the JHU and CDC datasets.
# This dataset is subsequently used in "06-time-series-plot-deaths.R" to generate
# the plot.
# 
#
# INPUT DATA FILE: "03-Data-Rodeo/01-analytic-data.rds", restricted access
# CDC data files (on my hard drive), public datasets from the Johns Hopkins CRC
# for the crude mortality and case counts.
#
# OUTPUT DATA FILE: "03-Data-Rodeo/02-long-form-analytic-data.rds"
#
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


library( tidyverse )
library( rvest )      # for webscraping html tables
library( RCurl )      # url
library( lubridate )  # date/time wrangling
library( tidycensus )

# import helper functions
source( "R/utils.R" )


### (0.0) Web Scrape for MMWR HTML Tables ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

url <-  "https://ibis.doh.nm.gov/resource/MMWRWeekCalendar.html"

df.list <- url %>% 
  read_html() %>% 
  html_nodes("table") %>% 
  html_table( fill = T ) 

mmwr <- do.call("cbind", df.list )
colnames( mmwr ) <- mmwr[ 1, ]

# intervals 
int.2020 <- interval(ymd("2020-01-01"), ymd("2020-12-31"))
int.2021 <- interval(ymd("2021-01-01"), ymd("2021-12-31"))
int.total <- interval(ymd("2020-03-25"), ymd("2021-03-25"))

# put MMWR table together
mmwr.df <- data.frame( mmwr[ -1, ] ) %>% # remove additional header row
  select( `MMWR.Week`, `X2021`, `X2020` ) %>% # select years of interest
  pivot_longer( values_to = "end.date", cols = c( `X2021`, `X2020` ) ) %>% # long format
  select( mmwr.week = MMWR.Week, end.date ) %>%
  
  # first, some regex to clean up date variables
  mutate( end.date = str_replace_all( end.date, "^(\\d\\/)", "0\\1"),
          end.date = str_replace_all( end.date, "(?<=\\/)(\\d\\/)", "0\\1"),
          
          # convert to date format
          end.date = as.Date( end.date, format = "%m/%d/%Y" ) ) %>%
  arrange( end.date ) %>%
  group_by( end.date ) %>%
  mutate( week = cur_group_id( ) ) %>%
  ungroup() %>% 
  data.frame() %>% 
  na.omit()

# now add start date as the day after the end date of the last period

for( i in 1: nrow( mmwr.df ) ) {
  
  mmwr.df[ i, "start.date" ] <- mmwr.df[ i-1, "end.date" ] + 1

  }

# remove NAs at tail ends of rnge
mmwr.df <- mmwr.df %>% na.omit()


# Create a frame with date and MMWR Week no. to join to data ##

dat.mmwr <- data.frame()
for( i in 1:nrow( mmwr.df ) ){
  
  # store dates and week information
  date.s<-  mmwr.df[ i, "start.date" ] 
  date.e <- mmwr.df[ i, "end.date" ] 
  
  mmwr.wk <- mmwr.df[ i, "mmwr.week" ] 
  week.no <- mmwr.df[ i, "week" ] 
  
  
  # make date sequence for range between start and end date of each MMWR week
  date.range <- seq.Date(from = as.Date( date.s ) , 
                        to = as.Date( date.e ), 
                        by = 'days') 
  
  # make a frame
  d.int <- data.frame( date = date.range,
                       mmwr = mmwr.wk,
                       week = week.no,
                       week.start = date.s,
                       week.end = date.e )
  
  # iteratively update out frame
  dat.mmwr <- rbind( dat.mmwr, d.int )
}
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (1.0) CDC Restricted Access Data Import and Wrangling ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


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


## (1.2) Put into long format ( long step, may take several hours) ##

cov.data <- readRDS( "03-Data-Rodeo/01-analytic-data.rds")
week.prep <- distinct( dat.mmwr, week, week.start, week.end )

# initialize frame to store data from each iteration
long.raw.cdc <- data.frame()
start.wk <- 13 # start MMWR week containing 3/25/2020 (official start is 3/22/2020)
max.wk <- length( unique( week.prep$week.end ) ) # last week to use for study
times.it <- vector() # vector to hold times of each iteration

for ( i in start.wk:max.wk ){
 
  start.time <- Sys.time( ) # track start time
  
  # use `death_data_prep` fct to prep data for each MMWR month and rowbind to previous iteration of data
  d.month <- death_data_prep( df = d.raw, mort.date.start = week.prep[ week.prep$week == start.wk  , "week.start"],
                                       mort.date.end = week.prep[ week.prep$week == i , "week.end"],
                                       case.date.start = week.prep[ week.prep$week == start.wk , "week.start"],
                                       case.date.end = week.prep[ week.prep$week == i , "week.end"] ) %>%
    mutate( week = week.prep[ week.prep$week == i , "week"],
            week.start = week.prep[ week.prep$week == i , "week.start"],
            week.end = week.prep[ week.prep$week == i , "week.end"] ) %>%
    arrange( fips, week ) 
  
  long.raw.cdc <- rbind( long.raw.cdc, 
                          d.month )
  
  
  
  end.time <- Sys.time( )  # track end time
  print( end.time - start.time )  # print time elapsed
  times.it <- c( times.it, end.time - start.time )
  est.time.left <- mean( times.it )*((nrow( week.prep ) - start.wk) - (i-start.wk+1 ) ) / 60
  print( paste0( "Iteration ", ((i-start.wk)+1), "/", length( seq_along(start.wk:max.wk) ), " complete." ) )  # track loop iteration
  print( paste0( "Estimated time remaining: ", est.time.left," minutes" ) )
  }

# save to local hard drive to not have to run above loop again
saveRDS( long.raw.cdc, "/Volumes/T7 Shield/Arthur-Lab/CDC-COVID-FI/CDC Restricted Access Data/long-form-week-case-mortality-data.rds")

# load in saved file
long.raw.cdc <- readRDS( "/Volumes/T7 Shield/Arthur-Lab/CDC-COVID-FI/CDC Restricted Access Data/long-form-week-case-mortality-data.rds")


## (1.3) Age-Standardize the Death Counts ##


# method for direct standardization is described here: 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3406211/

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

# convert columns to numeric 
for( i in which( str_detect( colnames( cnty.ref ), '.TOT' ) ) ){
  cnty.ref[ , i ] <- as.numeric( cnty.ref[ , i ] )
}

# reclassify age columns
long.cdc.3 <- left_join( long.raw.cdc, cnty.ref, by = 'fips' ) %>%
  mutate( `0-9 Years` = AGE04.TOT + AGE59.TOT,
          `10-19 Years` = AGE1014.TOT + AGE1519.TOT,
          `20-29 Years` = AGE2024.TOT + AGE2529.TOT,
          `30-39 Years` = AGE3034.TOT + AGE3539.TOT,
          `40-49 Years` = AGE4044.TOT + AGE4549.TOT,
          `50-59 Years` = AGE5054.TOT + AGE5559.TOT,
          `60-69 Years` = AGE6064.TOT + AGE6569.TOT,
          `70-79 Years` = AGE7074.TOT + AGE7579.TOT,
          `80+ Years` = AGE8084.TOT + AGE85PLUS.TOT ) 

# compute age-group-specific crude mortality rates

levs.agegrp <- levels( factor( long.cdc.3$age_group ) )

long.cdc.4 <- long.cdc.3 %>%
  mutate( grp.age.pop = ifelse( age_group == levs.agegrp[ 1 ], `0-9 Years`,
                                ifelse( age_group == levs.agegrp[ 2 ], `10-19 Years`,
                                        ifelse( age_group == levs.agegrp[ 3 ], `20-29 Years`,
                                                ifelse( age_group == levs.agegrp[ 4 ], `30-39 Years`,
                                                        ifelse( age_group == levs.agegrp[ 5 ], `40-49 Years`,
                                                                ifelse( age_group == levs.agegrp[ 6 ], `50-59 Years`,
                                                                        ifelse( age_group == levs.agegrp[ 7 ], `60-69 Years`,
                                                                                ifelse( age_group == levs.agegrp[ 8 ], `70-79 Years`,
                                                                                        ifelse( age_group == levs.agegrp[ 9 ], `80+ Years`, NA ) ) ) ) ) ) ) ) ) ) %>%
  mutate( crude.mort.rate = ( cum.deaths.cdc / grp.age.pop ) )


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


## (1.4) Multiply county crude mortality rates by age group population in reference population ##

# fix levels of age_group variable to ensure smooth join
long.cdc.4 %>%
  mutate( age_group = as.factor( age_group ) )

# recode levels
long.cdc.4$age_group <- recode_factor( long.cdc.4$age_group, "0 - 9 Years" = "0-9 Years",
                                "10 - 19 Years" = "10-19 Years",
                                "20 - 29 Years" = "20-29 Years",
                                "30 - 39 Years" = "30-39 Years",
                                "40 - 49 Years" = "40-49 Years",
                                "50 - 59 Years" = "50-59 Years",
                                "60 - 69 Years" = "60-69 Years",
                                "70 - 79 Years" = "70-79 Years" )

long.cdc.5 <- left_join( long.cdc.4, ref.age, by = "age_group" ) %>%
  
  # compute expected deaths in reference population
  mutate( exp.deaths.cdc = crude.mort.rate*ref.pop ) %>% 
  group_by( fips, week ) %>%
  
  # sum expected deaths within county fips and divide by reference population size
  # Note: mort rates are multiplied by 100,000 to represent deaths per 100,000 in the population
  mutate( age.adj.mort.rate = 100000*sum( exp.deaths.cdc ) / sum( ref.age$ref.pop ),
          crude.deaths.cdc = sum( cum.deaths.cdc)) %>%
  ungroup() %>%
  distinct( fips, week, week.start, week.end, age.adj.mort.rate, cases.cdc, crude.deaths.cdc )

## (1.5) Join mortality rates back to working data ##

long.cdc.6 <- left_join( long.cdc.5, cov.data[,1:47], by = "fips") %>%
  mutate( deaths.cdc = ceiling( ( age.adj.mort.rate / 100000 )*pop ) ) %>%
  filter( fips %in% cov.data$fips )   # filter out some rows summarizing state-level data

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (2.0) JHU Mortality and Case Data Import and Wrangling ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (2.1) Raw mortality data import ##
jhu.mort <- read.csv( text = getURL( "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv" ) )


# modify column names in raw mortality data for easy selection
colnames( jhu.mort ) <- str_remove( colnames( jhu.mort ), "X" ) %>%
  str_replace_all(., "(^\\d\\.)", "0\\1") %>%              # append "0" to months less than "10"
  str_replace_all(., "(?<=\\d\\d\\.)(\\d)\\.", "0\\1\\." ) # append "0" to days less than "10"

# subset date column names and convert to dates for easy subset
cols.as.dates <- as.Date( colnames( jhu.mort )[ 13:length( colnames( jhu.mort ) ) ], "%m.%d.%y" ) 

## (2.2) Columns to extract (based on study period--in this case the last MMWR week extends into 2022, so will need to get those data as well)

# start and end dates ( see section 1.2 above to get `start.week` and `week.prep`)
date.start <- week.prep[ week.prep$week == start.wk  , "week.start"]
date.end <- week.prep[ nrow( week.prep ) , "week.end"]

date.cols.to.extract <- cols.as.dates[ cols.as.dates <= date.end & cols.as.dates >= date.start ] %>%
  format( ., "%m.%d.%y" )
  
# now subset those date columns and other relevant columns
# and put data into long format ( each row is a county-week combination )
jhu.mort.long <- jhu.mort %>% # change `date` column to class date
  select( fips = FIPS, any_of( date.cols.to.extract ) ) %>% 
  pivot_longer( cols = c( contains( ".20" ),
                          contains( ".21" ),
                          contains( ".22" )),
                values_to = "deaths.jhu", 
                names_to = "date" ) %>%
  mutate( date = mdy( date ) ) %>%
  
  # join MMWR week columns
  left_join( ., dat.mmwr, by = "date" ) %>%
  group_by( fips, week ) %>%
  filter( deaths.jhu == max( deaths.jhu ) ) %>%      # this accounts for adjustments in death counts since counts are cumulative
  ungroup() %>%
  # keep distinct combinations of the following columns
  distinct( fips, week, mmwr, week.start, week.end, deaths.jhu ) %>%
  arrange( fips, week ) %>%
  mutate( fips = ifelse( fips < 10000, paste0( "0", fips ), fips ) )  # convert fips to character and 5-digits before merge
  



## (2.3) Raw case data import ##

jhu.case <- read.csv( text = getURL( "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv" ) )

# modify column names in raw mortality data for easy selection
colnames( jhu.case ) <- str_remove( colnames( jhu.case ), "X" ) %>%
  str_replace_all(., "(^\\d\\.)", "0\\1") %>%              # append "0" to months less than "10"
  str_replace_all(., "(?<=\\d\\d\\.)(\\d)\\.", "0\\1\\." ) # append "0" to days less than "10"

# pivot case data into long form ( each row is a county-week combination )
jhu.case.long <- jhu.case %>% # change `date` column to class date
  select( fips = FIPS, any_of( date.cols.to.extract ) ) %>% 
  pivot_longer( cols = c( contains( ".20" ), 
                          contains( ".21" ),
                          contains( ".21" )),
                values_to = "cases.jhu", 
                names_to = "date" ) %>%
  mutate( date = mdy( date ) ) %>%
  
  # join MMWR week columns
  left_join( ., dat.mmwr, by = "date" ) %>%
  group_by( fips, week ) %>%
  filter( cases.jhu == max( cases.jhu ) ) %>%      # this accounts for adjustments in death counts since counts are cumulative
  ungroup() %>%
  # keep distinct combinations of the following columns
  distinct( fips, week, mmwr, week.start, week.end, cases.jhu ) %>%
  arrange( fips, week ) %>%
  mutate( fips = ifelse( fips < 10000, paste0( "0", fips ), fips ) )  # convert fips to character and 5-digits before merge


# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (3.0) Merge All Data into Final Frame ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

( d.out <- left_join( jhu.mort.long, jhu.case.long ) %>%
  left_join( long.cdc.6, . )  %>%
  mutate( week = week - min( .$week ) + 1 ) )  %>%  # fix week number to start at "1" (there are two week columns in the dataset, `week` which pertains to week in the study period, and `mmwr` which is the mmwr week)
  
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (4.0) Save Long-Form Dataset ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
  
saveRDS( "03-Data-Rodeo/02-long-form-analytic-data.rds" )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------


