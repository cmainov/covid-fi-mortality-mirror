###---------------------------------------------------
###   01-COVARIATE DATA IMPORT, WRANGLING, AND MERGE
###---------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we will import covariate data from a number of public sources. 
# Certain variables will be kept and the data will be wrangled. The common identifier
# across datasets will be FIPS codes.
# 
# INPUT DATA FILES: 
# Several from publicly available sources. Links are provided below and will be used
# for download and accessing the data.
#
# OUTPUT DATA FILE: "02-Data-Wrangled/01-covariate-merge.rds"
#
# Resources: 
# i. Computation of county health index scores: https://pubmed.ncbi.nlm.nih.gov/34100936/ 
# ii. County Health Rankings Variable Descriptions: https://www.countyhealthrankings.org/explore-health-rankings/county-health-rankings-measures
# iii. Using `tidycensus`: https://walker-data.com/tidycensus/articles/basic-usage.html
# iv. Gini Index and COVID-19 mortality, correlational analysis: https://link.springer.com/article/10.1007/s11606-020-05971-3
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

library( tidycensus )  # For Census API Queries
library( tidyverse )
library( tigris )      # for TIGER shapefiles
library( downloader )  # for downloading external files
library( haven )       # for reading in SAS and other foreign data files
library( readxl )

options( tigris_class = "sf" )
options( tigris_use_cache = TRUE )



### American Community Survey (U.S. Census Bureau) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# Census variable codes
dp <- load_variables(2021, "acs5", cache = TRUE) %>%
  filter( str_detect( label, "high school" ) )

## Query TIGER Shapefiles API for shapefiles/land area ##
cy <- counties()


## Query ACS/Census API to obtain key variables, wrangle, and compute percentages ##

acs <- get_acs( geography = "county", year = 2020,
         variables = c( pop = "B01003_001", pov = "B17001_002",
                        female = "B01001_026", male = "B01001_002",
                        black = "B02001_003", asian = "B02001_005",
                        hispanic = "B03001_003", am.indian = "B02001_004",
                        nat.hawaii = "B02001_006", nh.white = "B03002_003",
                        median.age = "B01002_001", unemp = "B27011_008",
                        foreign.born = "B05002_013", gini.index = "B19083_001",
                        health.ins = "B992701_003", disable = "B10052_002",
                        no.veh = "B08014_002" ) ) %>%
  
  # arrange data in wide format
  pivot_wider( id_cols = GEOID, names_from = variable, values_from = estimate ) %>%
  left_join( ., cy %>% select( GEOID, ALAND ) ) %>% # join shapefile data
  
  # compute percentages
  mutate( native = am.indian + nat.hawaii,
          poverty.rate = ( pov / pop )*100,
          perc.female = ( female / pop )*100,
          perc.male = ( male / pop )*100,
          perc.hisp = ( hispanic / pop )*100,
          perc.native = ( native / pop )*100,
          perc.nh.white = ( nh.white / pop )*100,
          perc.black = ( black / pop )*100,
          perc.asian = ( asian / pop )*100,
          ALAND = ALAND / 2589988.1103,          # convert land area from sq. meters to sq. miles
          pop.density = pop / ALAND,
          unemp.rate = (unemp / pop )*100,
          perc.fb = ( foreign.born / pop )*100,
          no.vehic = (no.veh / pop )*100,
          no.health.insur = ( 1 - (health.ins / pop ) )*100,
          disability = (disable / pop )*100 ) %>%
  select( -c( pov, female, male, hispanic, am.indian, native, nat.hawaii, nh.white,
              ALAND, unemp, foreign.born, health.ins, disable, no.veh ) ) %>% # remove some variables no longer needed
  rename( fips = GEOID ) %>% # rename fips code column
  select( -geometry ) %>%  # drop geometry column
  data.frame() 

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### Atlas of Rural and Small Town America ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# file url (accessed 04 November 2022)
url <- "https://www.ers.usda.gov/webdocs/DataFiles/82701/RuralAtlasData23.xlsx?v=4875.1"

# download file and save in directory
download.file( url, 
              method = "curl", 
              destfile = "01-Data-Raw/rsa.xlsx" )

# merge 3 sheets from raw .xlsx file and join to ACS data
d.1 <- read_xlsx( "01-Data-Raw/rsa.xlsx", sheet = "Income") %>%
  left_join( ., read_xlsx( "01-Data-Raw/rsa.xlsx", sheet = "People") ) %>%
               left_join( ., read_xlsx( "01-Data-Raw/rsa.xlsx", sheet = "Jobs" ) ) %>%
  select( fips = FIPS, State, County, FIPS, AvgHHSize, Ed5CollegePlusPct, Ed1LessThanHSPct, 
          PerCapitaInc, MedHHInc, PctEmpTrans, PctEmpTrade ) %>%  # keep columns
  left_join( acs, . )        # join to ACS data

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### Area Health Resource File/American Medical Association ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# file url (accessed 04 November 2022)
url.2 <- "https://data.hrsa.gov//DataDownload/AHRF/AHRF_2020-2021_SAS.zip"

# download file into repository
download.file( url.2, destfile = "01-Data-Raw/ahrf.zip")

# unzip file in repository
unzip( "01-Data-Raw/ahrf.zip", exdir = "01-Data-Raw" )

# read in .sas7bdat file from repository and subset columns of interest
aha <- read_sas( "01-Data-Raw/AHRF_2020-2021_SAS/AHRF2021.sas7bdat" ) %>%

# total emergency physician variables (by age in 019) are f10759-19,f10760-19,f10761-19,f10762-19,f12090-19,12091-19,
# f14740-19 (DO osteopathic emergency doctors)
# f12592-19 (# of hospitals with on-campus emergency department)
# f12565-19 (# of hospitals with medical/surgical intensive care)
# f12566-19 (# of hospitals with cardiac intensive care)

  select(  f00002, f1075919, f1076019, f1076119, f1076219,
           f1209019, f1209119, f1474019, f1256519, f1256619,
           f1259219 ) %>%
  mutate( total.ed.phys = rowSums( .[-1] ) ) %>% # sum total emergency physicians (do not include fips column)
  select( fips = f00002, total.ed.phys, no.hosp.ed = f1259219, no.hosp.ms.icu = f1256519,
          no.hosp.card.icu = f1256619 ) # final keep variables

# join to working data set and compute ED physician:population ratio (interpreted as number of ED physicians per 10000 county inhabitants)
d.2 <- left_join( d.1, aha ) %>%
  mutate( ratio.pop.edp = ( total.ed.phys / pop )*10000 )

# remove data file given size >25 MB (for the purposes of uploading to GitHub)
file.remove( "01-Data-Raw/AHRF_2020-2021_SAS/AHRF2021.sas7bdat")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### UW Institute for Health Metrics and Evaluation Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## Download data ##

# url for respiratory mortality (accessed 4 November 2022)
url.3 <- "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_USA_COUNTY_RESP_DISEASE_MORTALITY_1980_2014_NATIONAL_XLSX.zip"

# url for cardiovascular mortality (accessed 4 November 2022)
url.4 <- "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_USA_COUNTY_CVD_MORTALITY_RATES_1980_2014_NATIONAL_XLSX.zip"

# url for cancer mortality (accessed 4 November 2022)
url.5 <- "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_USA_COUNTY_CANCER_MORTALITY_RATES_1980_2014_NATIONAL_XLSX.zip"

# download zip files
download.file( url.3, destfile = "01-Data-Raw/uw.resp.zip")
download.file( url.4, destfile = "01-Data-Raw/uw.cvd.zip")
download.file( url.5, destfile = "01-Data-Raw/uw.neoplasm.zip")

# unzip files into directory
unzip( "01-Data-Raw/uw.resp.zip", exdir="01-Data-Raw" )
unzip( "01-Data-Raw/uw.cvd.zip", exdir="01-Data-Raw" )
unzip( "01-Data-Raw/uw.neoplasm.zip", exdir="01-Data-Raw" )

# filenames for import (read from current directory)
filenm.resp <- dir( "01-Data-Raw")[ str_detect( dir( "01-Data-Raw"), "RESP" ) ]
filenm.cvd <- dir( "01-Data-Raw")[ str_detect( dir( "01-Data-Raw"), "CVD" ) ]
filenm.neoplasm <- dir( "01-Data-Raw")[ str_detect( dir( "01-Data-Raw"), "CANCER" ) ]


## Read in data ##

# chronic respiratory dz mortality
resp <- read_xlsx( paste0( "01-Data-Raw/", filenm.resp ),
           sheet = "Chronic respiratory diseases", # only the CRD sheet
           skip = 1, col_names = T ) %>%
  select( fips = FIPS, resp.mort = contains( "Mortality Rate, 2014" ) ) %>% # keep only 2014 data
  filter( fips > 1000 ) %>% # keep only county data (remove national and state data)
  mutate( resp.mort = as.numeric( str_extract( resp.mort, "^[^\\s\\()]*") ) ) # extract everything before space and open parenthesis (to remove confidence interval) and the convert to numeric
  
# cardiovascular dz mortality
cvd <- read_xlsx( paste0( "01-Data-Raw/", filenm.cvd ),
                       sheet = "Cardiovascular diseases", # only the CRD sheet
                       skip = 1, col_names = T ) %>%
  select( fips = FIPS, cvd.mort = contains( "Mortality Rate, 2014" ) ) %>% # keep only 2014 data
  filter( fips > 1000 ) %>% # keep only county data (remove national and state data)
  mutate( cvd.mort = as.numeric( str_extract( cvd.mort, "^[^\\s\\()]*") ) ) # extract everything before space and open parenthesis (to remove confidence interval) and the convert to numeric

# cancer mortality
neoplasm <- read_xlsx( paste0( "01-Data-Raw/", filenm.neoplasm ),
                  sheet = "Neoplasms", # only the CRD sheet
                  skip = 1, col_names = T ) %>%
  select( fips = FIPS, ca.mort = contains( "Mortality Rate, 2014" ) ) %>% # keep only 2014 data
  filter( fips > 1000 ) %>% # keep only county data (remove national and state data)
  mutate( ca.mort = as.numeric( str_extract( ca.mort, "^[^\\s\\()]*") ) ) # extract everything before space and open parenthesis (to remove confidence interval) and the convert to numeric


## Join and convert fips to character for merge with working data ##

d.3 <- left_join( resp, cvd ) %>%
  left_join( ., neoplasm ) %>%
  mutate( fips = ifelse( fips < 10000, paste0( "0", fips ), fips ) ) %>%
  
  # join to working data 
  left_join( d.2, . )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### USDA/ERS Food Environment Atlas ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## Download and unzip data ##

# url for FEA data (accessed 5 November 2022)
url.6 <- "https://www.ers.usda.gov/webdocs/DataFiles/80526/FoodEnvironmentAtlas.zip"

# download file and save in repository
download.file( url.6, 
               destfile = "01-Data-Raw/fea.zip" )
# unzip file into directory
unzip( zipfile = "01-Data-Raw/fea.zip", 
       files = "StateAndCountyData.csv",
       exdir = "01-Data-Raw")



## Merge 2 sheets from raw .xls file and join to working data ##

d.4 <- read.csv( "01-Data-Raw/StateAndCountyData.csv" ) %>%
  
  # arrange data into wide format
  pivot_wider( id_cols = FIPS, names_from = Variable_Code, values_from = Value ) %>%
  select( fips = FIPS, PCT.LACCESS.POP15 = PCT_LACCESS_POP15,
          PCT.SNAP17 = PCT_SNAP17,
          PCT.OBESE.ADULTS17 = PCT_OBESE_ADULTS17, 
          PCT.DIABETES.ADULTS13 = PCT_DIABETES_ADULTS13 ) %>% # keep variables
  
  # convert fips to character and 5-digits before merge
  mutate( fips = ifelse( fips < 10000, paste0( "0", fips ), fips ) ) %>%
  left_join( d.3, . ) # join to working data
  
# PCT.LACCESS.POP15 = Population, low access to store (%), 2015
# PCT.OBESE.ADULTS17 = Adult obesity rate, 2017*
# PCT.DIABETES.ADULTS13 = Adult diabetes rate, 2013
# PCT.SNAP17 = % population on SNAP benefits, 2017

# remove data file given size >25 MB (for the purposes of uploading to GitHub)
file.remove( "01-Data-Raw/StateAndCountyData.csv" )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### County Health Rankings Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## Download data ##

#  url for CHR data (accessed 20 June 2022)
url.7 <- "https://www.countyhealthrankings.org/sites/default/files/media/document/analytic_data2023_0.csv"

download.file( url.7, destfile = "01-Data-Raw/chr.csv" )


## Wrangle and join ##

# read in data and skip national summary data 
d.5 <- read.csv( "01-Data-Raw/chr.csv", 
                 header = T )[ -c( 1, 2 ), ]  %>%
  filter( ! Name %in% state.name ) %>% # remove state summary data
  
  # keep variables
  select( fips = X5.digit.FIPS.Code, 
          fi.perc.20 = Food.Insecurity.raw.value, 
          lim.acc.food.19 = Limited.Access.to.Healthy.Foods.raw.value,
          perc.adult.smoke.19 = Adult.Smoking.raw.value ) %>%
  
  # convert to percentage and numeric type
  mutate( fi.perc.20 = as.numeric( fi.perc.20 )*100,
          lim.acc.food.19 = as.numeric( lim.acc.food.19 )*100,
          perc.adult.smoke.19 = as.numeric( perc.adult.smoke.19 )*100 ) %>%
  left_join( d.4, . ) # join to working data

# fi.perc.20 = Percentage of population who lack adequate access to food 2020 (Map the Meal Gap)
# lim.acc.food.19 = Percentage of population who are low-income and do not live close to a grocery store.
# perc.adult.smoke.19 = Smoking prevalence (% Adults smoking), 2019 (BRFSS)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### NCHS Urban-Rural Classification Scheme ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## Download data ##

# file url (accessed 05 November 2022)
url.8 <- "https://www.cdc.gov/nchs/data/data_acces_files/NCHSURCodes2013.xlsx"

# download file and save in directory
download.file( url.8, 
               method = "curl", 
               destfile = "01-Data-Raw/nchs.xlsx" )

## Wrangle and join ##

d.6 <- read_xlsx( "01-Data-Raw/nchs.xlsx" ) %>%
  rename( nchs.2013.code = `2013 code`) %>%  # rename 2013 code column
  # classify urban and rural based on NCHS codes (see :https://www.cdc.gov/nchs/data/series/sr_02/sr02_166.pdf)
  mutate( urb.cat.code = as.factor( ifelse( nchs.2013.code %in% c( 1:4 ),"Metropolitan",
                  ifelse( nchs.2013.code %in% c( 5:6 ), 'Non-metropolitan', NA ) ) ),
          # convert original column with six levels to factor
          nchs.2013.code = as.factor( nchs.2013.code ) ) %>%
  select( fips = `FIPS code`, urb.cat.code, nchs.2013.code ) %>%  # keep binary column and column with all levels
  
  # convert fips to character and 5-digits before merge
  mutate( fips = ifelse( fips < 10000, paste0( "0", fips ), fips ) ) %>%
  left_join( d.5, . ) # join to working data
  
# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### 2020 Presidential Election Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

url.9 <- "https://raw.githubusercontent.com/tonmcg/US_County_Level_Election_Results_08-20/master/2020_US_County_Level_Presidential_Results.csv"

# download file and save in directory
download.file( url.9, destfile = "01-Data-Raw/elec-2020.csv" )

d.7 <- read.csv( "01-Data-Raw/elec-2020.csv" ) %>%
  # convert fips to character and 5-digits before merge 
  rename( fips = county_fips ) %>%
  mutate( fips = ifelse( fips < 10000, paste0( "0", fips ), fips ),
          elec.2020.margin = per_point_diff*100 ) %>% # select vote differential column ( positive number indicates a county won by republicans and a negative number indicates a democratic county)
  select( fips, elec.2020.margin ) %>%
  left_join( d.6, . ) %>%
  group_by( State ) %>% # group by state for computing state 2020 election margins
  mutate( elec.2020.margin.state = weighted.mean( x = elec.2020.margin,
                                                  w = pop ) ) %>% # compute state 2020 election margin as weighted mean of county 2020 election margins weighted by population size
  ungroup()

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### Clean-up Column Formatting ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# add periods to column names to separate key words
colnames( d.7 ) <- str_replace_all( colnames( d.7 ), 
                                    "([:lower:])(?=[:upper:])", "\\1\\.")
# regex: find lower case followed by an upper case; keep lower case and append a period after it

# add periods to column names to separate numbers from characters
colnames( d.7 ) <- str_replace_all( colnames( d.7 ), 
                                    "([:alpha:])(?=\\d)", "\\1\\.")
# regex: find letter followed by a digit; keep letter and append a period after it

# relocate state and county variables and make all columns same case
d.8 <- d.7 %>%
  relocate( State, .after = fips ) %>%
  relocate( County, .after = State ) %>%
  rename_all( tolower )   # make all column names lowercase

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### Create Health Index Composite Variable via PCA (1st PC) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# source: https://pubmed.ncbi.nlm.nih.gov/34100936/

## PCA ##

prhealth <- princomp( ~ ca.mort + cvd.mort + resp.mort + perc.adult.smoke.19 +
                     pct.diabetes.adults.13 + pct.obese.adults.17, data = d.8,
                   cor = T, scores = T )


## Examine loadings ##

prhealth$loadings # all factors load positively on 1st principal component


## Generate scores ##

scores <- data.frame( predict( prhealth, newdata = d.8 ) )%>%  # generate scores
  mutate( rowid = rownames( . ) ) %>%  # add row identifier for merge
  select( rowid, health.index = Comp.1 )  # keep first PC

# merge PC scores to working data 
( d.9 <- d.8 %>%
  mutate( rowid = rownames( . ) ) %>% # add row identifier for merge
  left_join( ., scores ) %>%          # join scores data to working data
  select( -rowid ) %>%                # remove rownames column
  relocate( health.index, .before = perc.adult.smoke.19 ) ) %>% # relocate new column
  group_by( state ) %>%
  mutate( health.index.state = weighted.mean( x = health.index,
                                              w = pop ) ) %>% # compute state mean health index by taking the mean of the counties health index score weighted by population
  ungroup() %>%
  relocate( health.index,state, .after = health.index.state ) %>% # relocate new column
  
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### Save ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS( "02-Data-Wrangled/01-covariate-merge.rds")
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

