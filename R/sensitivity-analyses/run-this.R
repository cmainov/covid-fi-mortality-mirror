# this script runs all the sensitivity analyses scripts. all results files are stored in their respective
# folders under 04-Tables-Figures/06-sensitivity-analyses

these.files <- dir( "R/sensitivity-analyses" )

sapply( paste0( "R/sensitivity-analyses/", these.files[1:3] ),  function( x ) source( x ) )

