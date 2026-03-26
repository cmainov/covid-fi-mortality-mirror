###---------------------------------------------------
###   RESULT-GENERATING AND OTHER HELPER FUNCTIONS
###---------------------------------------------------

#####################################
########## %notin% operator #########
#####################################
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
`%notin%` <- Negate( `%in%` )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

#####################################
###### list to rbind function #######
#####################################
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
list_it <- function( list ) {
  do.call( "rbind" , list )
}
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


#####################################
######### Nearest Function ##########
#####################################
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

nearest <- function( x, your.number ){
  # x = a vector of values (numeric)
  # your.number = the number you want to find the value closest to
  # returns the index for the entry nearest to `your.number` 
  
  which.min( abs( x - your.number ) )
  
}
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


#################################################
######### Equally Spaced Bins Function ##########
#################################################
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

eq_spaced_bins <- function( x = NULL, minimum = NULL, maximum = NULL, bins ){
  # x = vector of values you want to cateogrize into equally spaced bins
  # bins = no. of bins desired
  # returns a data frame with each bin, the bin start value and end value
  # NOTE: the start value for the subsequent row is the same as the end value for the 
  # previous row, and so the user needs to specify inclusive or exclusive notation
  # when using this data frame to assign values to bins 
  
  # checks 
  if( !is.null( x ) & ( !is.null( minimum ) | !is.null( maximum ) ) ){
    stop( "Either `x` or `minimum` and `maximum` must be supplied but not all three" )
  }
  
  if( is.null( x ) & ( ( is.null( minimum ) & !is.null( maximum ) ) |
                       ( !is.null( minimum ) & is.null( maximum ) ) ) ){
    stop( "One of `minimum` of `maximum` was supplied but the other value needed (`minimum` of `maximum`) was not supplied. Ensure both values are supplied or a vector, `x`, is supplied." )
  }
  
  # if a vector (`x`) is supplied do this 
  if ( !is.null( x ) ){
    
    bin.size <- ( max( x ) - min( x ) ) / ( bins )

    start <- min( x )
  }
  
  # if minimum and maximum values are supplied do this 
  if ( !is.null( minimum ) &  !is.null( maximum ) ){
    
    bin.size <- ( maximum - minimum ) / ( bins )
    
    start <- minimum
  }
  
  # initiate data frame to hold table
  hold.frame <- data.frame()
  for( i in 1:bins ){
    
    end <- start + bin.size
    
    hold.frame <- rbind( hold.frame,
                         data.frame( bin = i, 
                                     bin.start = start,
                                     bin.end = end ) )
    
    start <- end 
  }
  
  return( hold.frame )
}

# examples: 
# using `x`:
# eq_spaced_bins( x = qs, bins = 20 )
#
# specifying `minimum` and `maximum`:
# eq_spaced_bins( minimum = min( qs ), 
#                 maximum = max( qs ),
#                 bins = 20 )
# 
# specifying both sets of parameters which prompts error:
# eq_spaced_bins( x = qs,
#                 minimum = min( qs ), 
#                 maximum = max( qs ),
#                 bins = 20 )
#
# not supplying both `minimum` or `maximum` which prompts error:
# eq_spaced_bins(
#                 maximum = max( qs ),
#                 bins = 20 )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



####################################################################################################
#################################### Quantile Cutting Function #####################################
####################################################################################################

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

quant_cut<-function(var,x,df){
  
  xvec<-vector() # initialize null vector to store
  
  for (i in 1:x){
    xvec[i]<-i/x
  }
  
  qs<-c(min(df[[var]],na.rm=T), quantile(df[[var]],xvec,na.rm=T))
  
  df[['new']]=x+1 # initialize variable
  
  for (i in 1:(x)){
    df[['new']]<-ifelse(df[[var]]<qs[i+1] & df[[var]]>=qs[i],
                        c(1:length(qs))[i],
                        ifelse(df[[var]]==qs[qs==max(qs)],x,df[['new']]))
  }
  
  return(df[['new']])
}
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



####################################################################################################
#################################### Mortality Data Prep Function ##################################
####################################################################################################

# ---------------------------------------------------------------------------------------------------------------------------------------------------------


death_data_prep <- function( df, mort.date.start, mort.date.end,
                             case.date.start, case.date.end ){
  
  ## Aggregate cases at the level of county FIPS code ##
  
  d.1b <- df %>%
    group_by( county_fips_code ) %>%
    filter( cdc_case_earliest_dt >= format( case.date.start, format = '%Y-%m-%d' ) & 
              cdc_case_earliest_dt <= format( case.date.end, format = '%Y-%m-%d' ) ) %>%
    mutate( cases.cdc =  n() ) %>% # cumulative county COVID-19 cases between specified dates above
    ungroup() %>%
    select( fips = county_fips_code, cases.cdc ) %>%
    distinct() %>%
    mutate( fips = as.numeric( fips ),
            fips = ifelse( fips < 10000, paste0( "0", fips ), fips ) ) %>%
    filter( !str_detect( fips, "000$") ) %>%
    arrange( fips )
  
  ## Aggregate deaths at the level of county FIPS code and age group ##
  
  d.1 <- df %>%
    mutate( death_yn = as.numeric( recode( death_yn, `Yes` = 1, `Unknown` = 0, `No` = 0, `Missing` = 0 ) ) ) %>%
    filter( death_dt >= format( mort.date.start, format = '%Y-%m-%d' ) & 
              death_dt <= format( mort.date.end, format = '%Y-%m-%d' ) ) %>% # set max and min date for analysis
    group_by( county_fips_code, age_group ) %>%  # aggregate
    mutate( cum.deaths.cdc = sum( death_yn ) ) %>%
    select( fips = county_fips_code, age_group, cum.deaths.cdc ) %>%
    ungroup() %>%
    mutate( fips = as.numeric( fips ),
            fips = ifelse( fips < 10000, paste0( "0", fips ), fips ) ) %>%
    distinct() %>%
    left_join( ., d.1b, by = "fips" ) %>%     # join to cases data
    arrange( fips, age_group ) 
  
  
  ## Fill in "0's" for age groups with no deaths within a county ##
  
  # get levels for loop
  levs.age.grp <- levels( factor( d.1$age_group ) )
  levs.age.grp <- levs.age.grp[ levs.age.grp != "Missing" ]
  fips.codes <- levels( factor( cov.data$fips ) )
  
  # initialize data frame
  age.frame.out <- data.frame( age_group = NA,
                               fips = NA,
                               cum.deaths.cdc = NA,
                               cases.cdc = NA )
  
  # loop for replacement
  for( i in 1:length( fips.codes ) ){
    
    start.time <- Sys.time()
    
    cases.in.fips <- if (nrow( d.1b %>%filter( fips == fips.codes[i] ) %>% select( cases.cdc ) ) == 0 ) 0 else c( na.omit( d.1b %>%filter( fips == fips.codes[i] ) %>% select( cases.cdc ) ) )
    
    if ( fips.codes[i] %notin% levels( as.factor( filter( d.1, fips == fips.codes[i] )$fips ) ) ){
      
      age.frame.out <- na.omit( rbind( age.frame.out, data.frame( age_group = levs.age.grp,
                                                                  fips = rep( fips.codes[i], length( levs.age.grp ) ),
                                                                  cum.deaths.cdc = rep( 0, length( levs.age.grp ) ) ,
                                                                  cases.cdc = cases.in.fips ) ) )
    }
    
    else if ( fips.codes[i] %in% levels( as.factor( filter( d.1, fips == fips.codes[i] )$fips ) ) ){
      
      for( j in seq_along( levs.age.grp ) ){
        
        
        if ( levs.age.grp[j] %notin% levels( as.factor( filter( d.1, fips == fips.codes[i] )$age_group ) ) ){
          
          age.frame.out <- na.omit( rbind( age.frame.out, data.frame( age_group = levs.age.grp[j],
                                                                      fips = fips.codes[i],
                                                                      cum.deaths.cdc = 0,
                                                                      cases.cdc = cases.in.fips ) ) )
          
        }
      }
      
    }
    
    end.time <- Sys.time( )  # track end time
    print( end.time - start.time )  # print time elapsed
    print( paste0( "Iteration ", i, "/",length( fips.codes ), " complete." ) )  # track loop iteration
  }
  
  
  
  
  ## Row bind the original and new dataset and save intermediate dataset ##
  
  d.2 <- rbind( d.1, age.frame.out ) %>%
    arrange( fips, age_group ) %>%
    filter( age_group != 'Missing' ) %>%
    group_by( fips ) %>%
    mutate( cases.cdc = max( cases.cdc, na.rm = TRUE ),
            cases.cdc = ifelse( is.infinite( cases.cdc ),NA, cases.cdc ) )
  
  
  return( d.2 )
}
# ---------------------------------------------------------------------------------------------------------------------------------------------------------




####################################################################################################
################################## Table 1 (Continuous Variables) ##################################
####################################################################################################

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

tab1.var.mean<-function(var.name,df,table.var.name,strata.var=NULL,strata.level=NULL){
  
  if(is.null(strata.var)==T){
    df<- data.frame( df )
  }
  else {
    df<-data.frame( df[df[[strata.var]]==strata.level,] )
  }
  
  
  rowvar.name<-c(table.var.name)
  rowvar.mean<-c(paste0(round(mean(df[[var.name]],na.rm=T),digits=1),' (',round(sd(df[[var.name]]),digits=1),')'))
  
  partial.table<-data.frame(cbind(rowvar.name,rowvar.mean))
  colnames(partial.table)<-c('Characteristic','Frequency (%) or Mean (SD)')
  return(partial.table)
  
}
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



####################################################################################################
################################# Table 1 (Categorical Variables) ##################################
####################################################################################################

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

tab1.var.freq<-function(var.name,df,table.var.name,strata.var=NULL,strata.level=NULL){ #var.name is quoted string of how
  #variable is stored in dfset, df=is the dfset stored in R environment
  #and table.var.name is a character string of how that section of table 1 should be titled
  #strata.var is the variable, quoted, to stratify on, and strata.level is a quoted string
  #containing the level of strata.var that is to be examined
  
  
  df <- data.frame( df )
  
  
  
  if (is.null(strata.var)==T){
    
    df2<-data.frame( df )
    
    df2[[ var.name ]] <- factor( df2[[var.name]] )
  }
  
  if (is.null(strata.var)==F) {
    df2<-data.frame( df[df[[strata.var]] %in% strata.level,] )
    
    df2[[ var.name ]] <- factor(df2[[ var.name ]], 
                                levels =c( levels( factor( df2[[var.name]] ) ), levels( factor( df[[var.name]] ) )[which( levels( factor( df[[var.name]] ) ) %notin% levels( factor( df2[[var.name]] ) ) ) ] ))
    
  }
  
  rowvar.name<-vector()
  levelvec<-levels(df2[[var.name]])[order(levels(df2[[var.name]]))]
  
  for (i in 1:length(levelvec)){
    rowvar.name[i]<-paste0(levelvec[i])
  }
  
  # add to variable name header
  rowvar.name <- c(table.var.name,rowvar.name)
  
  rowvar.freq<-vector()
  for (i in 1:length(levelvec)){
    rowvar.freq[i]<-paste0(table(df2[[var.name]])[levelvec[i]],' (',round(100*table(df2[[var.name]])[levelvec[i]]/sum(table(df2[[var.name]]), na.rm = T),digits=1),')')
  }
  
  rowvar.freq<-c('',rowvar.freq)
  rowvar.freq<-ifelse(rowvar.freq=='NA (NA)',paste0('0 (0.0)'),rowvar.freq)
  
  partial.table<-data.frame(cbind(rowvar.name,rowvar.freq))
  colnames(partial.table)<-c('Characteristic','Frequency (%) or Mean (SD)')
  return(partial.table)
  
  
}
# ---------------------------------------------------------------------------------------------------------------------------------------------------------




####################################################################################################
##################################### Weighted Mean (Table 1) ######################################
####################################################################################################
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


tab1.var.mean<-function(var.name,data,table.var.name,strata.var=NULL,strata.level=NULL,weight=NULL){
  
  if(is.null(strata.var)==T){
    data<-data
  }
  else {
    data<-data[data[[strata.var]]==strata.level,]
  }
  
  invisible(require(radiant.data))
  
  if(is.null(weight)==F){
    rowvar.name<-c(table.var.name)
    rowvar.w.mean<-c(paste0(round(weighted.mean(data[[var.name]],w=data[[weight]],na.rm=T),digits=2),' (',round(weighted.sd(data[[var.name]],wt=data[[weight]],na.rm=T),digits=2),')'))
    partial.table<-data.frame(cbind(rowvar.name,rowvar.w.mean))
    colnames(partial.table)<-c('Characteristic','(Weighted) Mean (SD)')
  }
  
  if(is.null(weight)==T){
    rowvar.name<-c(table.var.name)
    rowvar.w.mean<-c(paste0(round(mean(data[[var.name]],na.rm=T),digits=2),' (',round(sd(data[[var.name]],na.rm=T),digits=2),')'))
    partial.table<-data.frame(cbind(rowvar.name,rowvar.w.mean))
    colnames(partial.table)<-c('Characteristic','(Weighted) Mean (SD)')
  }
  
  return(partial.table)
  
}

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




####################################################################################################
##################################### INLA Model Forest Plot #######################################
####################################################################################################

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

f_plot_inla <- function( models, term, x.names, title = NULL, x.label, x.axis.limits = NULL, limits.direction = NULL ){
  
  # models: a list with the model objects (INLA models only)
  # term: a character with the name of the variable to plot
  # x.names: a character vector with the names to put on the x axis (which is flipped in this case)
  # title: a character title for the plot
  # controls size of x axis text (which is on the y axis when)
  # x.label: a character, the x axis title
  # x.axis.limits: a numeric value that will multiply the x axis limits on the plot. must be used in conjunction with `limits.direction`. the default is the ggplot-default axis limits
  # limits.direction: a string, one of "left", "right", "both". it selects the axis limit to amplify.
  
  # dependencies: dplyr, ggplot2, scales
  
  plot.these <- data.frame() # initialize frame to prepare data for plotting
  for( i in 1:length( models ) ){
    
    d.plot <- models[[i]]$summary.fixed 
    
    
    
    # wrangle table and plot forest plot with ggplot2
    plot.these <- bind_rows( plot.these, d.plot %>%
                               mutate( row = rownames(.) ) %>%
                               filter( str_detect( row, term ) ) %>%
                               mutate( sig = as.factor( as.numeric( sign( `0.025quant`) == sign( `0.975quant`) ) ),
                                       x = x.names[i],
                                       y.short = round( exp( mean ), 2 ),
                                       y.min.short = round( exp( `0.025quant` ), 2 ),
                                       y.max.short = round( exp( `0.975quant` ), 2 ),
                                       wo = paste0( sprintf( "%.2f", y.short ), " (", sprintf( "%.2f", y.min.short ), " - ", sprintf( "%.2f", y.max.short ), ")")) ) # `sprintf` controls the number of decimal places to paste the number at (here 2 decimal places)
    
  }
  
  f.plot.inla <- ggplot( data = plot.these, 
                         mapping = aes( x = x, 
                                        y = exp( mean ), 
                                        ymin = exp( `0.025quant` ), 
                                        ymax = exp( `0.975quant` ), 
                                        label = wo,
                                        shape = sig ) ) +
    geom_pointrange( size = 1 ) +
    geom_errorbar( aes( ymin = exp( `0.025quant` ), ymax = exp( `0.975quant` ), width = 0 ) ) + # add error bars
    coord_flip() + 
    geom_text( mapping = aes( color = sig), hjust = -0.17, vjust = -0.7, family = "Avenir", size = 5.8 ) +
    scale_color_manual( values = c( "grey50", "black"), breaks = c( 0, 1 )) +
    scale_shape_manual(values = c( 1, 16 ), breaks = c( 1, 0 )) + # different shapes based on stat. significance
    geom_hline( yintercept = 1, lty = "dashed" ) +
    theme_classic() +
    theme( text = element_text( family = "Avenir" ),
           axis.title.y = element_blank(),
           plot.title = element_text( size = 15),
           axis.text.y =  element_text( size = 16.7 ),
           axis.text.x =  element_text( size = 16.7 ),
           axis.title.x =  element_text( size = 17 ),
           legend.position = "none") +
    ylab( x.label )  +
    ggtitle ( title ) +
    scale_y_continuous( 
      labels = scales::number_format(accuracy = 0.01, # the scales package allows us to control the digits following a decimal (here to two decimal places) and the decimal mark (here a ".")
                                     decimal.mark = '.' ) )
  
  if( !is.null( x.axis.limits ) ){
    
    if( is.null( limits.direction ) ) stop( '`limits.direction` must be one of "right", "left" or "both" if `x.axis.limits` is specified and not NULL.')
    
    scale.y <- layer_scales( f.plot.inla )$y$get_limits()
    
    if( limits.direction == "right" ) scale.y.low <- scale.y[1]
    if( limits.direction == "left" ) scale.y.up <- scale.y[2]
    if( limits.direction == "left" | limits.direction == "both" ) scale.y.low <- scale.y[1]*( 1 / x.axis.limits )
    if( limits.direction == "right" | limits.direction == "both" ) scale.y.up <- scale.y[2]*x.axis.limits
    
    f.plot.inla <- f.plot.inla + 
      scale_y_continuous( limits = c( scale.y.low, scale.y.up ),
                          labels = scales::number_format(accuracy = 0.01, # the scales package allows us to control the digits following a decimal (here to two decimal places) and the decimal mark (here a ".")
                                                         decimal.mark = '.' ) )
    
  }
  
  return( f.plot.inla )
  
}

# Example:
# f_plot_inla( models = list( res.jhu.so, res.cdc.so ),
#   term = "fi.perc.19", 
#   x.names = c( "Data Source: Johns Hopkins", "Data Source: CDC" ) )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



###########################################################################################
################################## Generate NA Color/Text for Legend ############################
###########################################################################################
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

generate_NA_legend <- function( na.text, color, text.size=17 ){
  # na.text = a string that will show as the value for the NA's in the plot
  # color = a ggplot2 color (a string)
  # dependencies: ggpubr, dplyr, ggplot2
  
  # generate a plot and extract the legend 
  legend.plot <- data.frame( main = c( rnorm(20), NA ) ) %>%
    mutate( main.two = ifelse( !is.na( main ), NA, na.text ),
            y = c( rnorm(20), NA ) ) %>%
    filter( is.na( main ) ) %>%
    
    ggplot(mapping = aes( x = main, y = y, color = main.two ) ) +
    geom_point( shape = 15, size = 5 ) +
    scale_color_manual( values = c( color ) ) +
    theme_classic() +
    theme(
      legend.text = element_text( family = "Avenir", size = text.size),
      legend.title = element_blank()
    ) 
  
  # extract the legend
  leg.p <- get_legend( legend.plot ) %>%
    ggpubr::as_ggplot(.) +
    theme(
      plot.background = element_rect(fill = 'black'))
  
  return( leg.p )
}

# example:
# generate_NA_legend( "Lost data", "blue")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



####################################################################################################
##################################### INLA Model Rel. Risk Plot ####################################
####################################################################################################

# ---------------------------------------------------------------------------------------------------------------------------------------------------------


plot_rr <- function( df, model, legend.scale = 0.3, share.legend, plot.this, 
                     which.fixed = NULL, E = FALSE, 
                     leg.position = "top", state.codes.omit,
                     na.text = "Incomplete Data" ){
  
  # df: data frame w/ geometries  used in INLA model
  # model: INLA model object
  # legend.scale: scales the size of the legend. default is 0.3
  # plot.this: a character string--one of "overall.risk", "fixed.effect", "spatial.effect", "unstructured.random.state", "unstructured.random.cty"
  # if the "fixed.effect" option is supplied to the plot.this argument, then this is a character string with the variable name of the fixed effect to plot
  # state.codes.omit = a character vector with the 2 digit states fips codes to omit from the shapefile
  # na.text = a character string. text for denoting the missing or omitted data for the missing data legend. defaults to "Incomplete Data",
  # E = a logical. was the E arguement specified in the INLA call? default is FALSE.
  # dependencies: viridis, ggpubr, cowplot, dplyr, INLA, sf
  
  if( plot.this == "unstructured.random.state" ) res.sp <- model$marginals.random$state[1:nrow( model$summary.random$state )] # extract posterior distributions for state unstructured random effects and remove NULL elements of list
  if( plot.this == "spatial.effect" ) res.sp <- model$marginals.random$re.s[1:nrow( df )] # extract posterior distributions for spatial effects
  if( plot.this == "unstructured.random.cty" ) res.sp <- model$marginals.random$fips[1:nrow( df )] # extract posterior distributions for spatial effects
  if( plot.this == "overall.risk" ) res.sp <- model$marginals.fitted.values[1:nrow( df )] # extract posterior distributions for fitted values
  if( plot.this == "fixed.effect" ){ res.sp <- model$marginals.fixed$`I(fi.perc.20/4)`
  
  
  # extract coefficient of term to plot
  this.coef <- model$summary.fixed %>%
    mutate( row = rownames(.) ) %>%
    filter( str_detect( row, which.fixed ) ) %>%
    select( mean ) %>% unlist()
  
  # compute beta*x_i
  zeta <- this.coef %*% as.vector( df[[which.fixed]] ) %>% as.vector() 
  
  df$rr <- zeta
  
  } 
  
  # compute expectation of the posterior marginals. see:  https://www.paulamoraga.com/book-geospatial/sec-inla.html?
  if( plot.this %in% c( "spatial.effect", 
                        "unstructured.random.cty")  ){
    
    zeta <- sapply( res.sp, function(x) inla.emarginal( fun = function(x) x, x ) ) # select mean (we will keep in log scale)
    
    df$rr <- zeta
    
  }
  
  if( plot.this == "unstructured.random.state" ){
    
    zeta <- sapply( res.sp, function(x) inla.emarginal( fun = function(x) x, x ) ) # select mean (we will keep in log scale)
    
    df <- data.frame( zeta ) %>%
      mutate( state = model$summary.random$state$ID ) %>%
      left_join( df, ., by = "state" ) %>% # add rownames for state random effects
      rename( rr = zeta )
    
  }
  
  if( plot.this == "overall.risk" ){
    
    zeta <- sapply( res.sp, function(x) inla.emarginal( fun = function(x) x, x ) )
    
    if( E ) zeta <- log( zeta )
    
    df$rr <- zeta
  }
  
  
  
  ## plot ##
  
  if( share.legend ){
    # extract coefficient of term to plot
    this.coef <- model$summary.fixed %>%
      mutate( row = rownames(.) ) %>%
      filter( str_detect( row, which.fixed ) ) %>%
      select( mean ) %>% unlist()
    
    random.1 <- model$marginals.random$state[1:nrow( model$summary.random$state )] # extract posterior distributions for state unstructured random effects and remove NULL elements of list
    random.2 <- model$marginals.random$re.s[1:nrow( df )] # extract posterior distributions for spatial effects
    random.3 <- model$marginals.random$fips[1:nrow( df )]  # extract posterior distributions for unstructured county residual
    fitted.vals <- model$marginals.fitted.values[1:nrow( df )]
    
    # commbine all values into single vector to get maximum and minimum among all values for the plot
    these.all <- c( sapply( random.1, function(x) inla.emarginal( fun = function(x) x, x ) ),
                    sapply( random.2, function(x) inla.emarginal( fun = function(x) x, x ) ),
                    sapply( random.3, function(x) inla.emarginal( fun = function(x) x, x ) ),
                    if( E )  log(  sapply( fitted.vals, function(x) inla.emarginal( fun = function(x) x, x ) ) ) else if( !E )  sapply( fitted.vals, function(x) inla.emarginal( fun = function(x) x, x ) ) , # extract posterior distributions for spatial effects
                    this.coef %*% as.vector( df[[which.fixed]] ) %>% as.vector() ) 
    
    # breaks in continuous scale
    qs <- quantile( these.all, (0:8)/8, na.rm = TRUE )
    qs <- setNames( qs, 
                    rep( "", length(qs) ) )
    
    # find entry closest to 0 and assign "0.00" to the names
    names( qs )[ nearest( qs, 0 ) ] <- "0.00"
    
  }
  
  if( !share.legend ){
    # breaks in continuous scale
    qs <- quantile( df$rr, (0:8)/8, na.rm = TRUE )
    qs <- setNames( qs, c( sprintf( "%.2f", round( min(df$rr, na.rm = TRUE ), 2 ) ), rep( "", 8 ), sprintf( "%.2f", round( max(df$rr, na.rm = TRUE ), 2 ) ) ) )
  }
  
  # pull county-level shapefiles from `tigris`
  t <- counties( cb = TRUE ) # USA Census tract shapefiles
  
  t.df <- t %>%
    mutate( state.code = str_extract( GEOID, "^.{2}" ) ) %>%
    left_join( ., st_drop_geometry(df) %>% select( -state.code ), by = c( "GEOID" = "fips") ) %>%
    dplyr::filter( !state.code %in% state.codes.omit )
  
  # generate main plot (i.e., the fixed or random effect we are mapping)
  main.plot <- ggplot( data = t.df ) +
    geom_sf(           
      mapping = aes( fill = rr ), color = NA ) +
    coord_sf( crs = st_crs( 2163 ), expand = T ) +
    
    scale_fill_gradientn(colors = viridis(9), # we will use the viridis color pallete for the maps
                         space = "Lab",
                         na.value = "red",
                         breaks = c(0),
                         labels = c(paste0( "0.00" ) ), # show where zero is on color scale/bar
                         guide = "colourbar",
                         aesthetics = "fill") +
    theme_minimal( ) +
    theme( 
      legend.key.size = unit( legend.scale, 'cm' ),
      plot.background = element_rect( color = "transparent" ) ) +
    guides( fill = guide_colorbar( title = "Rel. Risk",
                                   nbin = 10,
                                   ticks.colour = NA,
                                   frame.colour =  "black",
                                   barwidth = 35,
                                   barheight = 0.7,
                                   direction = "horizontal",
                                   label.position = "bottom") )
  
  
  
  # plot of psi (the proportion of spatial variance compared to that of the unstructured random effect)
  
  df2 <- model$summary.random$state %>%
    rename( state = ID )  %>%
    select( state, u.sd = sd) %>%
    left_join( df, ., by = "state" ) %>%
    mutate( v.sd = model$summary.random$re.s$sd,
            e.sd = model$summary.random$fips$sd,
            psi = v.sd / (u.sd + v.sd + e.sd) ) %>%
    st_drop_geometry() %>% select( -state.code ) %>%
    left_join( t.df, ., by = c( "GEOID" = "fips") ) %>%
    dplyr::filter( !state.code %in% state.codes.omit ) # filter only mainland states
  
  
  # plot 
  qs.v <- quantile( df2$psi, (0:9)/9, na.rm = TRUE )
  qs.v <- setNames( qs.v, c( sprintf( "%.2f", round( min(df2$psi, na.rm = TRUE ), 2 ) ), rep( "", 8 ), sprintf( "%.2f",round( max(df2$psi, na.rm = TRUE ), 2  ) ) ) ) 
  
  psi.plot <- ggplot( data = df2 ) +
    geom_sf(           
      mapping = aes( fill = psi ),
      color = NA, size = 0.5 ) +
    coord_sf( crs = st_crs( 2163 ), expand = T ) +
    scale_fill_fermenter( palette = 6, direction = 1, guide = "colourbar", # `scale_fill_fermenter` for discretizing a continuous scale
                          breaks = qs.v, na.value = "red" ) + # if missing, the polygon will show up as red
    theme_minimal( ) +
    ggtitle( TeX( "$\\psi$" ) ) +
    theme( legend.title = element_blank(),
           legend.key.size = unit( legend.scale, 'cm' ),
           plot.background = element_rect( color = "transparent" ) ) +
    guides( fill = guide_colorbar( ticks.colour = NA,
                                   frame.colour =  "black",
                                   barwidth = 10,
                                   barheight = 0.6,
                                   direction = "horizontal",
                                   label.position = "bottom") )
  
  # generate NA grob
  leg.na <- generate_NA_legend( na.text = na.text, color = "red", text.size = 17 )
  
  # add legend layer 
  this.psi.leg <- ggdraw( psi.plot ) + 
    draw_plot( leg.na, 0.45, 0.03, 0.05, 0.05, scale = 0.05 )
  
  
  
  
  return( list( main.plot = main.plot, psi.plot = this.psi.leg, psi.plot.no.leg = psi.plot,
                nmiss = nrow( t.df ) - nrow( df ),
                range.vector = qs ) )
}


# Example: plot_rr( df = d.poly.so, model = res.jhu.so )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




####################################################################################################
############################## Forward Variable Selection w/ DIC or WAIC ###########################
####################################################################################################
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

forwards_z <- function( iter = 1, 
                        oth.cov = NULL, 
                        formula.full,
                        data,
                        base.value.criterion = NULL,
                        main.x = NULL, 
                        criterion.threshold, 
                        criterion = "dic",
                        out.model.list = list(),
                        log.out = vector(),
                        car.prior = NULL,
                        ur.prior = NULL,
                        model.next = NULL,
                        ... ){
  # iter = iteration counter. defaulted to `1`.
  # oth.cov = a character vector of other covariates to consider for addition in the stepwise algorithm. set to NULL and specified within the function using the argument to `formula.full`. more for internal purposes than anything else.
  # main.x = a character string of the main explanatory variable to keep in the model. NULL is defaulted meaning every variable has a chance of getting added
  # data = a dataframe
  # base.value.criterion = for internal use. establishes the criterion value to compare models to for the first iteration.
  # criterion = a character string of the criterion to use. can be either "dic" (deviance information criterion) or "waic" (Watanabe-Akaike information criterion)
  # criterion.threshold = a numeric entry with the value of the criterion to use as a cutoff for removing a variable
  # out.model.list = list() is the default. This is the list that holds the results. DO NOT change this defaulted option!
  # log.out = an empty vector to store the running log of additions to the model selection procedure. DO NOT change the defaulted option!
  # car.prior =  spatially structured random effect precision prior specification. should be in list format (e.g., list(prec = list(prior = "loggamma", param = c(1, 0.0005)),initial = 4, fixed = F)). default is NULL which specified the former prior for the log precision. for the vector in the `param` argument, the first number is the shape parameter and the second is the rate parameter
  # ur.prior =  unstructured random effect precision prior specification. should in be list format (e.g., list(prec = list(prior = "gaussian", param = c(0, 0.002))) where first number in the vector is the mean and second is the precision.
  # model.next = the model object carried over to the next iteration.for internal purposes only. would not recommend changing.
  
  #dependencies: INLA, spdep, dplyr
  ## fit base model
  print( paste0(  "predictor #", iter, " models" ) )
  
  # formula components
  
  # these needed for each iteration
  resp <- model_terms( formula.full )$response.var
  offs <- model_terms( formula.full )$offset.term
  reffs <- model_terms( formula.full )$random.effects
  feffs <- model_terms( formula.full )$fixed.effects
  
  if( iter == 1 ) oth.cov <- feffs
  
  if( !is.null( main.x ) ){
    
    if( sum( main.x %in% oth.cov ) >= 1 ) oth.cov <- oth.cov[ !oth.cov %in% main.x ]
    
  }   
  
  
  part.z.models <- list() # list to store model objects
  meta.table.z <- data.frame() # initialize data frame to hold meta data to later be returned
  for( i in seq_along( oth.cov ) ){
    
    start.time <- Sys.time()
    
    if( iter == 1 & i == 1){
      f.base <- as.formula( paste0( resp, "~",
                                    paste0( main.x, collapse = "+" ),"+",
                                    paste0( reffs, collapse = "+" ) %>%
                                      str_replace( ., '\\"', "'" ) ) )
      
      
      model.base <- NULL
      while( is.null( model.base ) ){ # condition for continuing to loop (iterations until model converges)
        
        model.base <- 
          tryCatch(
            inla( f.base, family = "nbinomial", data = data.frame(data) , # tryCatch to keep loop from breaking when error in INLA occurs
                  control.predictor = list( compute = TRUE ), scale = TRUE,
                  control.fixed = list( mean = 0.0, prec = 0.001 ), E = E,
                  control.compute = list( dic = TRUE, waic = TRUE, return.marginals.predictor = TRUE ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
            
            error = function( e )
              NULL ) # assign model1 as NULL if `inla` function breaks
      }
      
      
      # set base value from base model (dic or waic based on input to `criterion`)
      if( criterion == "dic" ){
        base.crit <- model.base$dic$dic
      }
      
      if( criterion == "waic" ){
        base.crit <- model.base$waic$waic
        
      }
    }
    
    # if null value for criterion is not NULL then assign the value based on supplied argument
    if( (!is.null( base.value.criterion )) & iter != 1 ){
      base.crit <- base.value.criterion
    }
    
    # fit thde model without the ith covariate and take note of criterion value
    f.z <- as.formula( paste0( resp, "~",
                               paste0( main.x,collapse = "+"), "+",
                               paste0( oth.cov[i], collapse = "+" ),"+",
                               paste0( reffs, collapse = "+" ) ) %>%
                         str_replace( ., '\\"', "'" ) )
    model.z <- NULL
    runs <- 0
    while( is.null( model.z ) ){ # condition for continuing to loop (iterations until model converges)
      
      model.z <- 
        tryCatch(
          inla( f.z, family = "nbinomial", data = data.frame(data) , # tryCatch to keep loop from breaking when error in INLA occurs
                control.predictor = list( compute = TRUE ), scale = TRUE,
                control.fixed = list( mean = 0.0, prec = 0.001 ), E = E,
                control.compute = list( dic = TRUE, waic = TRUE, return.marginals.predictor = TRUE ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
          
          error = function( e )
            NULL ) # assign model1 as NULL if `inla` function breaks
      
      # run tracker
      runs <- runs + 1
    }
    
    # criterion value from `model.z`
    if( criterion == "dic" ){
      model.crit <- model.z$dic$dic
    }
    
    if( criterion == "waic" ){
      model.crit <- model.z$waic$waic
    }
    
    end.time <- Sys.time()
    
    # meta data table
    meta.table.z <- rbind( meta.table.z,
                           data.frame( model = i,
                                       criterion.used = criterion,
                                       fixed.eff = oth.cov[i],
                                       model.criterion.value = model.crit,
                                       base.criterion.value = base.crit, 
                                       criterion.diff = base.crit - model.crit, 
                                       runs.required = runs ) )
    
    # store model in list for later use
    part.z.models[[oth.cov[i]]] <- model.z
    
    
    print( end.time - start.time )
    print( paste0( "STEP ", iter, ": iteration ", i, "/", length( oth.cov ), " complete. ", runs, " run(s) needed for convergence." ) ) 
    
  }
  
  # see which predictors was worst
  
  # BEST predictor that satisfies criterion threshold
  
  best.x.z <- NULL # set to null; if the predictor with the minimum criterion value satisfies the criterion threshold value, then `best.x.z` will not be NULL; otherwise it will
  
  if( max( meta.table.z$criterion.diff ) >= criterion.threshold ){
    
    best.x.z <- meta.table.z[ which( meta.table.z$criterion.diff == max( meta.table.z$criterion.diff ) & meta.table.z$criterion.diff >= criterion.threshold ), "fixed.eff" ]
    
    step.log.add <- print( paste0( "ADD PREDICTOR: ", best.x.z ))
    
    log.out[[paste0(" add step = ",iter )]] <- step.log.add
  }
  
  
  if( is.null( best.x.z ) & iter == 1 ){
    model.keep <- model.base
  }
  
  if( is.null( best.x.z ) & iter != 1 ){
    model.keep <- model.next
  }
  
  if( !is.null( best.x.z ) ){
    model.keep <- part.z.models[[best.x.z]]
  }
  
  # covariates we keep for the next iteration
  these.keep <- oth.cov[ !oth.cov %in% best.x.z ]
  
  # output list with all the information and deliverables from each iteration
  iter.out.list <- list( step = iter,
                         out.model = model.keep,
                         add.meta.table = meta.table.z,
                         new.worst.x = best.x.z,
                         remaining.x = these.keep,
                         log = log.out)
  
  # append the output lists from each iteration into one large list to return at the end
  out.model.list[[paste0( "step = ",iter )]] <- c( iter.out.list )
  
  # if is.null( best.x.z ) then it means no predictors satisfied the minimum threshold value and the function stops
  if( is.null( best.x.z ) ){
    
    print( "No other variables satisfied criterion minimums. Returning all objects corresponding to previous steps.")
    return( out.model.list )
    
  }
  
  # update iteration 
  iter <- iter + 1
  
  
  if( criterion == "dic" ) new.base.dic <- model.keep$dic$dic
  if( criterion == "waic" ) new.base.dic <- model.keep$dic$dic
  
  
  
  # recursively call the function until is.null( worst.x.z ) == TRUE
  forwards_z( iter = iter, 
              oth.cov = these.keep, 
              criterion.threshold = criterion.threshold,
              criterion = criterion,
              out.model.list = out.model.list,
              main.x = c( main.x, best.x.z ),
              data = data,
              base.value.criterion = new.base.dic,
              log.out = log.out,
              resp = resp,
              offs = offs,
              reffs = reffs,
              formula.full = formula.full,
              model.next = model.keep,
              car.prior = car.prior,
              ur.prior = ur.prior 
  )
}


# ---------------------------------------------------------------------------------------------------------------------------------------------------------




###########################################################################################
################################## Spatial Fractional Variance ############################
###########################################################################################
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# this function allows us to compute the % of geographic heterogeneity explained by the spatially structured and unstructured random effects in the INLA models

inla_var_explained <- function( model.fit ){
  # model.fit = an INLA model object
  
  # proportion of variance explained by spatial random effect see:  https://stats.stackexchange.com/questions/350235/how-to-convert-estimated-precision-to-variance
  
  # this code is reproduced from: http://inla.r-inla-download.org/r-inla.org/case-studies/Blangiardo-et-al-2012/Code/LondonSuicides.R
  
  n <- length( model.fit$summary.fitted.values$mean )
  
  mat.marg <- matrix( NA, nrow = n, ncol = 1000 )
  m <- model.fit$marginals.random$re.s
  
  # distributions of each marginal distribution of upsilon
  for ( i in 1:n ){
    u <- m[[i]]
    s <- inla.rmarginal( 1000, u )
    mat.marg[i,] <- s
  }
  # compute empirical variance
  var.space <- mean( apply( mat.marg, 2, sd ) )^2 # mean of the standard deviations (for each row) and then square
  
  # compute expectation of unstructured random effects' variance
  var.state <- inla.emarginal( function(x) 1/x,
                               model.fit$marginals.hyperpar$`Precision for state` )
  
  var.cty <- inla.emarginal( function(x) 1/x,
                             model.fit$marginals.hyperpar$`Precision for fips` )
  
  # estimate proportions of variance explained by each random component
  prop.var = var.space / ( var.space + var.state + var.cty )
  
  prop.var.state = var.state / ( var.space + var.state + var.cty )
  
  prop.var.cty = var.cty / ( var.space + var.state + var.cty )
  
  
  return( data.frame( prop.sp.explained = prop.var,
                      prop.state.explained = prop.var.state,
                      prop.cty.explained = prop.var.cty,
                      unstructured.random.var.state = var.state,
                      unstructured.random.var.cty = var.cty,
                      spatial.random.var = var.space ) )
  # returns the proportion of variance explained by the spatially structured random effect,
  # the unstructured random effect's variance posterior mean, and the spatially structured random effect's variance posterior mean
  # as well as the proportions of each the variance components on the residual variance
}

# example
# inla_var_explained( res.jhu )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



####################################################################################################
##################################### INLA CAR Model Fitting ######################################
####################################################################################################
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# d.jhu = the jhu dataset 
# d.cdc = the cdc dataset
# which.model = one of "jhu" or "cdc". it is the model ( cdc or jhu ) you choose to fit (they must be fit separately.)
# car.prior =  spatially structured random effect precision prior specification. should be in list format (e.g., list(prec = list(prior = "loggamma", param = c(1, 0.0005)),initial = 4, fixed = F)). default is NULL which specified the former prior for the log precision. for the vector in the `param` argument, the first number is the shape parameter and the second is the rate parameter
# ur.prior =  unstructured random effect precision prior specification. should in be list format (e.g., list(prec = list(prior = "gaussian", param = c(0, 0.002))) where first number in the vector is the mean and second is the precision.
# formula.jhu = default is NULL for jhu model ("final model"). if an alternative formula is desired, it can be put here. (see code below for default covariates). default model is the convolution model.
# formula.cdc = default is NULL for cdc model. if an alternative formula is desired, it can be put here. (see code below for default covariates). default model is the convolution model.
#   ***IMPORTANT***: DO NOT specify the `hyper` option in the formula for any random effects. instead, use the car.prior and ur.prior arguments to this function to specify priors for the random effects. otherwise, the function will break.

# E = a character string of the variable name to compute the relative risk for. this variable should be constant thorughout the dataset and should be specified as a a mean or other summary measure. see `help(inla)` for more information 
# term = a character string of the variable name to return the exponentiated fixed effect estimate and 95% credible interval for
# null.model.formula = an object of class `formula` that specifies the formula for the null model. necessary for creating the plot of explained geographic heterogeneity. default is NULL (i.e., no plot is returned but the percentage of explained variance by the structural and unstructured random effects is returned)
# model.2.formula = an object of class `formula` that specifies the formula for the "basic" model. necessary for creating the plot of explained geographic heterogeneity. default is NULL (i.e., no plot is returned but the percentage of explained variance by the structural and unstructured random effects is returned)
# model.3.formula = an object of class `formula` that specifies the formula for the "basic+state" model. necessary for creating the plot of explained geographic heterogeneity. default is NULL (i.e., no plot is returned but the percentage of explained variance by the structural and unstructured random effects is returned)
# model.4.formula = an object of class `formula` that specifies the formula for the "demographic" model. necessary for creating the plot of explained geographic heterogeneity. default is NULL (i.e., no plot is returned but the percentage of explained variance by the structural and unstructured random effects is returned)

# x.legend.pos = controls position of missing data legend along x axis. defaults to 0.73 
# na.text = a character string. text for denoting the missing or omitted data for the missing data legend. defaults to "Incomplete Data",

# dependencies: inla, ggplot2, ggpattern, ggpubr, dplyr, reshape2, viridis, cowplot, ggspatial

car_inla_analysis <- function( d.jhu = NULL, d.cdc = NULL , which.model, car.prior = NULL,
                               ur.prior = NULL,
                               formula.jhu = NULL, 
                               formula.cdc = NULL,
                               E = NULL,
                               null.model.formula = NULL,
                               model.2.formula = NULL,
                               model.3.formula = NULL,
                               model.4.formula = NULL,
                               term,
                               state.codes.omit,
                               na.text = "Incomplete Data",
                               x.legend.pos = 0.73 ){
  
  total.time.start <- Sys.time()
  
  ## hyperparameter priors ##
  
  # log-gamma distribution with specified hyperparameters for spatial variance
  if( is.null( car.prior ) ){
    
    car.prior <- list( prec = list( prior = "loggamma", param = c(1, 0.0005), initial = 4 ) )
    
  }
  
  
  # normal distribution with specified hyperparameters for unstructured random effect
  if( is.null( ur.prior ) ){
    
    ur.prior <- list( prec = list( prior = "loggamma", param = c(1, 0.0005), initial = 4 ) )
    
  }
  
  print( "Model-fitting...BEGIN")
  if( which.model == "jhu" ){
    
    # create neighborhood list
    nb <- poly2nb( d.jhu )
    
    td <- tempdir() # create temporary directory
    
    # INLA-formatted neighbors list (adjacency matrix)
    am <- nb2INLA( paste( td, "inla-mat.adj", sep="/" ), nb )
    
    # create graph object ( I don't like using super assignment operators within function environments but there is no other way to get the g vector to work)
    g <<- inla.read.graph( filename = paste( td, "inla-mat.adj", sep="/" ) )
    
    # create indices in data for the spatial random effects
    d.jhu$re.s <- 1:nrow( d.jhu )
  }
  
  
  if( which.model == "cdc" ){
    # create neighborhood list
    nb.cdc <- poly2nb( d.cdc )
    
    td <- tempdir() # create temporary directory
    
    # INLA-formatted neighbors list (adjacency matrix)
    am.cdc <- nb2INLA( paste( td, "inla-mat-cdc.adj", sep="/" ), nb.cdc )
    
    # create graph object ( I don't like using super assignment operators within function environments but there is no other way to get the g vector to work)
    g.cdc <<- inla.read.graph( filename = paste( td, "inla-mat-cdc.adj", sep="/" ) )
    
    d.cdc$re.s <- 1:nrow( d.cdc )
  }
  
  
  # create indices in data for the state random effects
  # d.jhu <- d.jhu %>% group_by( state ) %>% 
  #   dplyr::mutate( re.state = cur_group_id() ) %>%
  #   ungroup()
  # 
  # d.cdc <- d.cdc %>% group_by( state ) %>% 
  #   dplyr::mutate( re.state = cur_group_id() ) %>%
  #   ungroup()
  
  # specify formula
  if( is.null( formula.jhu ) ){
    
    f.jhu <- deaths.jhu ~ I(fi.perc.20/4) + log( offset(pop) ) + pct.emp.trans + 
      perc.black + perc.female + perc.hisp + perc.nh.white + perc.native + pop.density +
      ed.1less.than.hspct + ed.5college.plus.pct + pct.emp.trade + median.age + perc.asian +
      poverty.rate + perc.vaccinated + gini.index + avg.hhsize + elec.2020.margin +
      ratio.pop.edp + health.index + f( re.s, model = "besag", graph = g, scale.model = T, hyper = car.prior ) +
      f( state, model = "iid", hyper = ur.prior )
    
  }
  
  if( !is.null( formula.jhu ) ){
    
    f.parts <- model_terms( formula.jhu )
    
    f.char <- paste0( f.parts$response.var, "~", {if (length( f.parts$offset.term ) != 0) paste0(f.parts$offset.term , "+")}, 
                      paste0( f.parts$fixed.effects, collapse = "+" ),
                      "+", paste0( f.parts$random.effects, collapse = "+" ) ) %>%
      str_replace_all( ., "model = 'besag',",  "model = 'besag', hyper = car.prior,") %>%
      str_replace_all( ., "model = 'iid'",  "model = 'iid', hyper = ur.prior")
    
    f.jhu <- formula( f.char )
  }
  
  if( is.null( formula.cdc ) ){
    
    f.cdc <- deaths.cdc ~ I(fi.perc.20/4) + log( offset(pop) ) + pct.emp.trans + 
      perc.black + perc.female + perc.hisp + perc.nh.white + perc.native + pop.density +
      ed.1less.than.hspct + ed.5college.plus.pct + pct.emp.trade + perc.asian +
      poverty.rate + perc.vaccinated + gini.index + avg.hhsize + elec.2020.margin +
      ratio.pop.edp + health.index + f( re.s, model = "besag", graph = g.cdc, scale.model = T, hyper = car.prior ) +
      f( state, model = "iid", hyper = ur.prior )
  }
  
  if( !is.null( formula.cdc ) ){
    
    f.parts <- model_terms( formula.cdc )
    
    f.char <- paste0( f.parts$response.var, "~", {if (length( f.parts$offset.term ) != 0) paste0(f.parts$offset.term , "+")}, 
                      paste0( f.parts$fixed.effects, collapse = "+" ),
                      "+", paste0( f.parts$random.effects, collapse = "+" ) ) %>%
      str_replace_all( ., "model = 'besag',",  "model = 'besag', hyper = car.prior,") %>%
      str_replace_all( ., "model = 'iid'",  "model = 'iid', hyper = ur.prior")
    
    f.cdc <- formula( f.char )
    
  }
  
  ### call INLA ###
  
  print( "Fit Final model...BEGIN" )
  
  ## JHU data ##
  # no variable selection
  if( which.model == "jhu" ){
    
    if( !is.null( E ) ) d.jhu$E <- d.jhu[[E]] # set E value
    
    
    res.jhu <- NULL # initialize `res.jhu` so that while loop below runs
    runs <- 0 # initialize no. of runs
    while( is.null( res.jhu ) ){
      res.jhu <-  tryCatch( 
        
        inla( f.jhu, family = "nbinomial", data = data.frame( d.jhu ),
              control.predictor = list( compute = TRUE ), E = E, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
              control.compute = list( return.marginals.predictor = TRUE ), # compute marginals
              control.fixed = list( mean = 0.0, prec = 0.001 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
        
        error = function( e )
          NULL ) # assign model2 as NULL if `inla` function breaks
      
      # run tracker
      runs <- runs + 1
      
    }
    
    
    
    # save proportion of variance explained by random effects
    p.exp <- inla_var_explained( res.jhu )
    
    print( "Fit Final model...DONE" )
    
    
    ## model 4: demographic model ##
    
    # null model
    if( !is.null( model.4.formula ) ){
      
      print( "Fit Demographic model...BEGIN" )
      
      demo.model <- NULL # initialize `demo.model` so that while loop below runs
      runs <- 0 # initialize no. of runs
      while( is.null( demo.model ) ){
        demo.model <-  tryCatch( 
          
          inla( model.4.formula, family = "nbinomial", data = data.frame( d.jhu ),
                control.predictor = list( compute = TRUE ), E = E, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
                control.fixed = list( mean = 0.0, prec = 0.0001 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
          
          error = function( e )
            NULL ) # assign demo.model as NULL if `inla` function breaks
        
        # run tracker
        runs <- runs + 1
        
      }
      
      # save proportion of variance explained by random effects
      p.exp.demo.model <- inla_var_explained( demo.model )
      
      print( "Fit Demographic model...DONE" )
      
    }
    
    ## null model and model #2 if specified ##
    
    # null model
    if( !is.null( null.model.formula ) ){
      
      print( "Fit Null model...BEGIN" )
      
      null.model <- NULL # initialize `null.model` so that while loop below runs
      runs <- 0 # initialize no. of runs
      while( is.null( null.model ) ){
        null.model <-  tryCatch( 
          
          inla( null.model.formula, family = "nbinomial", data = data.frame( d.jhu ),
                control.predictor = list( compute = TRUE ), E = E, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
                control.fixed = list( mean = 0.0, prec = 0.0001 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
          
          error = function( e )
            NULL ) # assign null.model as NULL if `inla` function breaks
        
        # run tracker
        runs <- runs + 1
        
      }
      
      # save proportion of variance explained by random effects
      p.exp.null.model <- inla_var_explained( null.model )
      
      print( "Fit Null model...DONE" )
      
    }
    
    # model 2
    if( !is.null( model.2.formula ) ){
      
      print( "Fit Basic model...BEGIN" )
      
      model.2 <- NULL # initialize `model.2` so that while loop below runs
      runs <- 0 # initialize no. of runs
      while( is.null( model.2 ) ){
        model.2 <-  tryCatch( 
          
          inla( model.2.formula, family = "nbinomial", data = data.frame( d.jhu ),
                control.predictor = list( compute = TRUE ), E = E, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
                control.fixed = list( mean = 0.0, prec = 0.001 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
          
          error = function( e )
            NULL ) # assign model.2 as NULL if `inla` function breaks
        
        # run tracker
        runs <- runs + 1
        
      }
      
      # save proportion of variance explained by random effects
      p.exp.model.2 <- inla_var_explained( model.2 )
      
      print( "Fit Basic model...DONE" )
      
    }
    
    # model 3 (basic + state variables)
    if( !is.null( model.3.formula ) ){
      
      print( "Fit Basic + State model...BEGIN" )
      
      model.3 <- NULL # initialize `model.2` so that while loop below runs
      runs <- 0 # initialize no. of runs
      while( is.null( model.3 ) ){
        model.3 <-  tryCatch( 
          
          inla( model.3.formula, family = "nbinomial", data = data.frame( d.jhu ),
                control.predictor = list( compute = TRUE ), E = E, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
                control.fixed = list( mean = 0.0, prec = 0.1 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
          
          error = function( e )
            NULL ) # assign model.2 as NULL if `inla` function breaks
        
        # run tracker
        runs <- runs + 1
        
      }
      
      # save proportion of variance explained by random effects
      p.exp.model.3 <- inla_var_explained( model.3 )
      
      print( "Fit Basic + State model...DONE" )
      
    }
  }
  
  
  ## CDC data ##
  # no variable selection
  if( which.model == "cdc" ){
    
    if( !is.null( E ) ) d.cdc$E <- d.cdc[[E]] # set E value
      
      res.cdc <- NULL # initialize `res.cdc` so that while loop below runs
      runs <- 0 # initialize no. of runs
      while( is.null( res.cdc ) ){
        res.cdc <-  tryCatch( 
          
          inla( f.cdc, family = "nbinomial", data = data.frame( d.cdc ),
                control.predictor = list( compute = TRUE ), E = E,  scale = TRUE, # set scale = TRUE so that models with different priors can be compared
                control.compute = list( return.marginals.predictor = TRUE ), # compute marginals
                control.fixed = list( mean = 0.0, prec = 0.001 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
          
          error = function( e )
            NULL ) # assign model2 as NULL if `inla` function breaks
        
        # run tracker
        runs <- runs + 1
        
      }
    
    # save proportion of variance explained by random effects
    p.exp <- inla_var_explained( res.cdc )
    
    ## null model and model #2 if specified ##
    
    # null model
    if( !is.null( null.model.formula ) ){
      
      print( "Fit Null model...BEGIN" )
      null.model <- NULL # initialize `null.model` so that while loop below runs
      runs <- 0 # initialize no. of runs
      while( is.null( null.model ) ){
        null.model <-  tryCatch( 
          
          inla( null.model.formula, family = "nbinomial", data = data.frame( d.cdc ),
                control.predictor = list( compute = TRUE ), E = E, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
                control.compute = list( return.marginals.predictor = TRUE ), # compute marginals
                control.fixed = list( mean = 0.0, prec = 0.001 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
          
          error = function( e )
            NULL ) # assign null.model as NULL if `inla` function breaks
        
        # run tracker
        runs <- runs + 1
        
      }
      
      # save proportion of variance explained by random effects
      p.exp.null.model <- inla_var_explained( null.model )
      
      print( "Fit Null model...DONE" )
      
    }
    
    ## model 4: demographic model ##
    
    # null model
    if( !is.null( model.4.formula ) ){
      
      print( "Fit Demographic model...BEGIN" )
      
      demo.model <- NULL # initialize `demo.model` so that while loop below runs
      runs <- 0 # initialize no. of runs
      while( is.null( demo.model ) ){
        demo.model <-  tryCatch( 
          
          inla( model.4.formula, family = "nbinomial", data = data.frame( d.cdc ),
                control.predictor = list( compute = TRUE ), E = E, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
                control.fixed = list( mean = 0.0, prec = 0.0001 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
          
          error = function( e )
            NULL ) # assign demo.model as NULL if `inla` function breaks
        
        # run tracker
        runs <- runs + 1
        
      }
      
      # save proportion of variance explained by random effects
      p.exp.demo.model <- inla_var_explained( demo.model )
      
      print( "Fit Demographic model...DONE" )
      
    }
    
    
    # model 2
    if( !is.null( model.2.formula ) ){
      
      print( "Fit Basic model...BEGIN" )
      
      model.2 <- NULL # initialize `model.2` so that while loop below runs
      runs <- 0 # initialize no. of runs
      while( is.null( model.2 ) ){
        model.2 <-  tryCatch( 
          
          inla( model.2.formula, family = "nbinomial", data = data.frame( d.cdc ),
                control.predictor = list( compute = TRUE ), E = E, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
                control.fixed = list( mean = 0.0, prec = 0.1 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
          
          error = function( e )
            NULL ) # assign model.2 as NULL if `inla` function breaks
        
        # run tracker
        runs <- runs + 1
        
      }
      
      # save proportion of variance explained by random effects
      p.exp.model.2 <- inla_var_explained( model.2 )
      
      print( "Fit Basic model...DONE" )
      
    }
    
    # model 3 (basic + state variables)
    if( !is.null( model.3.formula ) ){
      
      print( "Fit Basic + State model...BEGIN" )
      
      model.3 <- NULL # initialize `model.2` so that while loop below runs
      runs <- 0 # initialize no. of runs
      while( is.null( model.3 ) ){
        model.3 <-  tryCatch( 
          
          inla( model.3.formula, family = "nbinomial", data = data.frame( d.cdc ),
                control.predictor = list( compute = TRUE ), E = E, scale = TRUE, # set scale = TRUE so that models with different priors can be compared
                control.fixed = list( mean = 0.0, prec = 0.1 ) ), # set mean (0) and precision (0.001 -- or SD = 1000) priors for fixed effects
          
          error = function( e )
            NULL ) # assign model.2 as NULL if `inla` function breaks
        
        # run tracker
        runs <- runs + 1
        
      }
      
      # save proportion of variance explained by random effects
      p.exp.model.3 <- inla_var_explained( model.3 )
      
      print( "Fit Basic + State model...DONE" )
      
    }
  }
  
  
  
  ### fixed effect information to return ###
  
  this.model <- if( which.model == "jhu" ) res.jhu else if( which.model == "cdc" ) res.cdc 
  
  fixed.return <- d.plot <- this.model$summary.fixed %>%
    mutate( row = rownames(.) ) %>%
    filter( str_detect( row, term ) ) %>%
    mutate( 
      y.short = round( exp( mean ), 2 ),
      y.min.short = round( exp( `0.025quant` ), 2 ),
      y.max.short = round( exp( `0.975quant` ), 2 ),
      wo = paste0( y.short, " (", y.min.short, " - ", y.max.short, ")" ) ) %>%
    select( variable = row, exp.estimate.95.ci = wo ) # returns variable and the exponentiated estimate of the posterior mean and the 95% credible interval
  
  rownames( fixed.return ) <- NULL # set rownames null to not be redundant
  
  print( "Model-fitting...DONE")
  print( "Map Decomposition...BEGIN" )
  
  ### map decomposition ###
  
  ## parameters to pass to `Map` which will make maps using `plot_rr` for each of the fixed, unstructured random, spatial effects, and psi plot
  these.plots <- c( "fixed.effect", "spatial.effect", "overall.risk", "unstructured.random.state", "unstructured.random.cty" )
  labss <- c( "$\\gamma_{FI}x_{ij}$", "$\\upsilon_{i}$", "$\\mu_{ij}/E_{i} \\ ^b$", "$u_{j}$", "$e_{i}$" ) # parameters we are plotting
  
  leg.pos <- c( rep( "bottom", 5 ) ) # legend position
  title.this <- c( rep( TRUE, 5 ) ) # title logical (for input into the `mapply` below)
  
  top.marg <- c( rep( 0.3, 5 ) )
  bottom.marg <- c( rep( 0.1, 3 ), rep( 0.5, 3 ))
  
  
  ## now apply with `Map` ##
  map.d <- Map( function( x, label, l, title.in, title.logic, top.margin, bottom.margin ) 
  { 
    
    # `plot_rr` function to generate each of the maps with some slight modifications
    plots <- plot_rr( df = st_as_sf( if ( which.model == "jhu" ) d.jhu else if (which.model == "cdc") d.cdc ), 
                      model = if( which.model == "jhu" ) res.jhu else if( which.model == "cdc" ) res.cdc, 
                      legend.scale = 0.4, 
                      share.legend = TRUE,
                      leg.position = l,
                      state.codes.omit = state.codes.omit,
                      plot.this = x, 
                      which.fixed = "fi.perc.20",
                      E = !is.null( E ) )
    
    # we modify the "main" plot and then return the psi plot along with the main plot at the end of the Map script
    out.rr.plot <- plots$main.plot + 
      labs( fill = unname( TeX( label ) ) ) +
      ggtitle( unname( TeX( title.in ) ) ) +
      ( this.theme <- theme( 
        title = element_text( size = 8.7,
                              color = "red"),
        legend.title = if( !title.logic ) element_text( color = "grey41", size = 18 ) else element_blank(),
        legend.text = element_text( family = "Avenir", size = 24 ),
        plot.title = element_text( color = "grey41", size = 34 ),
        axis.text.x = element_text( family = "Avenir", size = 14 ),
        axis.text.y = element_text( family = "Avenir", size = 14 ),
        plot.margin = unit( c(top.margin,0.1,bottom.margin,0.1), 'cm' ) ) )  #unit(c(top, right, bottom, left), units)
    
    return( list( main.plot = out.rr.plot, 
                  nmiss = plots$nmiss,
                  range.vector = plots$range.vector,
                  psi.plot = plots$psi.plot.no.leg +
                    ggtitle( unname( TeX( "$\\psi \\ ^a$" ) ) ) +
                    this.theme) )
    
  } , x = these.plots, label = labss, l = leg.pos, title.in = labss,
  title.logic = title.this, top.margin = top.marg, bottom.margin = bottom.marg
  )
  
  # arrange
  map.decomp <- ggarrange( map.d$fixed.effect$main.plot +
                             # add north-pointing compass rose
                             ggspatial::annotation_north_arrow( location = "tl", which_north = "true", 
                                                                pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                                                                style = north_arrow_fancy_orienteering, width = unit(1.3, "cm"), 
                                                                height = unit(1.7, "cm") ) +
                             annotation_scale( height = unit(0.14, "cm"), style = "bar" ), # adds map scale bar
                           map.d$spatial.effect$main.plot + 
                             # add north-pointing compass rose
                             ggspatial::annotation_north_arrow( location = "tl", which_north = "true", 
                                                                pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                                                                style = north_arrow_fancy_orienteering, width = unit(1.3, "cm"), 
                                                                height = unit(1.7, "cm") ) +
                             theme( axis.text.y = element_blank()), 
                           map.d$spatial.effect$psi.plot +
                             # add north-pointing compass rose
                             ggspatial::annotation_north_arrow( location = "tl", which_north = "true", 
                                                                pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                                                                style = north_arrow_fancy_orienteering, width = unit(1.3, "cm"), 
                                                                height = unit(1.7, "cm") ) +
                             theme( axis.text.y = element_blank()),
                           map.d$overall.risk$main.plot +
                             # add north-pointing compass rose
                             ggspatial::annotation_north_arrow( location = "tl", which_north = "true", 
                                                                pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                                                                style = north_arrow_fancy_orienteering, width = unit(1.3, "cm"), 
                                                                height = unit(1.3, "cm") ) +
                             annotation_scale( height = unit(0.14, "cm"), style = "bar" ),
                           map.d$unstructured.random.state$main.plot + 
                             # add north-pointing compass rose
                             ggspatial::annotation_north_arrow( location = "tl", which_north = "true", 
                                                                pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                                                                style = north_arrow_fancy_orienteering, width = unit(1.3, "cm"), 
                                                                height = unit(1.7, "cm") ) +
                             theme( axis.text.y = element_blank()),
                           map.d$unstructured.random.cty$main.plot + 
                             # add north-pointing compass rose
                             ggspatial::annotation_north_arrow( location = "tl", which_north = "true", 
                                                                pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                                                                style = north_arrow_fancy_orienteering, width = unit(1.3, "cm"), 
                                                                height = unit(1.7, "cm") ) + 
                             theme( axis.text.y = element_blank()), 
                           ncol = 3, nrow =2, common.legend = TRUE, legend = "bottom") +
    cowplot::draw_figure_label( label = if( which.model == "jhu" ) "Data Source: Johns Hopkins" else if ( which.model == "cdc" ) "Data Source: CDC", 
                                position = "bottom.left",
                                family = "Avenir", color = "grey44",
                                size = 17 )
  
  # generate missing value grob for legend
  leg.p <- generate_NA_legend( na.text = paste0( na.text, " (n = ", map.d$fixed.effect$nmiss, ")" ),
                               color = "red", text.size = 22 )
  
  leg.h <- get_legend(map.d$spatial.effect$psi.plot)
  
  # add missing value for main color bar
  map.decomp <- ggdraw( map.decomp ) +
    draw_plot( leg.p, x.legend.pos, 0.021, 0.05, 0.05, scale = 0.1 )
  
  # add legend for psi plot color bar
  map.decomp <- ggdraw( map.decomp ) +
    draw_plot( leg.h, 0.81, 0.935, 0.05, 0.05, scale = 0.1 )
  
  # add maximum and minimum to color bar
  map.decomp <- map.decomp + 
    annotate( "text",
              x = 0.354,
              y = 0.018,
              label = sprintf( "%.2f", min( map.d$spatial.effect$range.vector, na.rm = T ) ), 
              size = 9,family = "Avenir" ) +
    annotate( "text",
              x = 0.648,
              y = 0.018,
              label = sprintf( "%.2f", max( map.d$spatial.effect$range.vector, na.rm = T ) ), 
              size = 9,family = "Avenir" ) 
  
  map.decomp <- annotate_figure( map.decomp, bottom = "" ) # adds buffer space at bottom of plot
  
  
  print( "Map Decomposition...DONE" )
  
  
  print( "Other Plots...BEGIN" )
  
  ### no. of neighbors and spatial variance plot ###
  
  sp.var.df <- data.frame( sp.var = this.model$summary.random$re.s$sd, nnbs = if( which.model == "jhu" ) g$nnbs else if( which.model == "cdc" ) g.cdc$nnbs )
  
  sp.var.nnbs.plot <- ggplot( data = sp.var.df, aes( x = nnbs, y = sp.var ) ) +
    geom_point( shape = 6, size = 1, alpha = 0.5 ) +
    geom_smooth( method = loess, color = "red", fill = "khaki3" ) + # `color` controls line color and `fill` controls ribbon color 
    theme_minimal() +
    theme( text = element_text( family = "Avenir" ) ) +
    ylab( TeX( "$SD(\\upsilon_{i})$" ) ) +
    xlab( "No. of Neighbors" ) +
    scale_x_discrete( limits = c( 0,3,5,7,10, max( sp.var.df$nnbs ) ) ) # labels at the tick marks we want
  
  
  ## percent of geographic heterogeneity explained by random effects plot ##
  
  # bind results `inla_explained_var` fct
  df.exp.var <- bind_rows( p.exp.null.model, p.exp.model.2, p.exp.model.3,
                           p.exp.demo.model, p.exp ) %>%
    mutate( model = c( "Null Model",
                       "Basic Model",
                       "Basic + State Model",
                       "Demographic Model",
                       "Final Model" ) )
  
  df.resid <- cbind( melt( df.exp.var %>%
                             select( model, spatial.random.var, unstructured.random.var.state, unstructured.random.var.cty  ), 
                           id.var= "model", value.name = "variance" ) %>% select( model, variance ),
                     melt( df.exp.var %>%
                             select( model, prop.sp.explained, prop.state.explained, prop.cty.explained  ), 
                           id.var= "model", value.name = "perc.num" ) %>% select( -model ) ) %>%
    mutate( sp.var.perc = paste0( round( perc.num, 4 )*100, "%" ),
            perc.num = perc.num*100,
            sp.var.perc = ifelse( variable == "prop.sp.explained" & perc.num >= 10, sp.var.perc,
                                  ifelse( variable == "prop.sp.explained" & perc.num < 10, "",
                                          ifelse( variable == "prop.state.explained" & perc.num < 10, "",
                                                  ifelse(variable == "prop.state.explained" & perc.num >= 10, sp.var.perc,
                                                         ifelse( variable == "prop.cty.explained" & perc.num < 10, "",
                                                                 ifelse(variable == "prop.cty.explained" & perc.num >= 10, sp.var.perc, NA)))))),
            
            color = ifelse( variable == "prop.sp.explained" & perc.num >= 10, "black",
                            ifelse( variable == "prop.sp.explained" & perc.num < 10, "#FDE725FF",
                                    ifelse( variable == "prop.cty.explained" & perc.num < 90, "#440154FF",
                                            ifelse(variable == "prop.cty.explained" & perc.num >= 90, "white",
                                                   ifelse( variable == "prop.state.explained" & perc.num < 10, "#440154FF",
                                                           ifelse(variable == "prop.state.explained" & perc.num >= 10, "white", NA)))))))
  
  # generate plot #
  
  
  if( !is.null( null.model.formula ) & !is.null( model.2.formula ) & !is.null( model.3.formula ) &
      !is.null( model.4.formula ) ){
    
    summed.df <- df.resid[ 1:5, "variance" ] + df.resid[ 6:10, "variance" ] + df.resid[ 11:15, "variance" ]
    
    max.limit <- max( summed.df )
    
    
    # get limit of the y axis with following expression
    upper.x <- round( (max.limit*1.05) / 5, 2)*5
    
    # get limits vector and exclude zero
    these.limits <- signif( seq( 0, upper.x, length.out = 6 ), 2 ) %>% .[ !.== 0 ]
    
    # generate vector of y axis limits and the character version for the plot itself
    if ( upper.x < 3 ) {
      these.limits.character <- 
        paste0( sprintf( "%.2f", these.limits ) ) # use `sprintf` to ensure formatting of digits is correct
    } else if  ( upper.x >= 3 ) these.limits.character <- paste0( sprintf( "%.1f", these.limits ) )
    
    # now plot
    explained.var.plot <- ggplot( df.resid, aes (x = model, y = variance, fill = variable ) ) + 
      geom_bar( stat = "identity", color = "black" ) + 
      scale_fill_manual( breaks = levels(as.factor( df.resid$variable) ),
                         values = c("#FDE725FF", "#440154FF", "#1F968BFF" ),
                         labels = c( unname( TeX( "Spatially Structured\nRandom Effect\n(County-Level)")),
                                     "Unstructured Random\nEffect (State-Level)", 
                                     "Unstructured Random\nEffect (County-Level)")) +
      theme_classic() +
      theme( legend.title = element_blank(),
             text = element_text( family = "Avenir" ),
             axis.title.x = element_blank(),
             axis.title.y = element_text( size = 15),
             legend.text = element_text(margin = margin(t = 2.5), size = 14),
             axis.text.x = element_text( angle = -30, hjust = 0.1, vjust = 0.5,
                                         size = 14),
             axis.text.y = element_text( size = 13 ) ) +
      scale_x_discrete( limits = c("Null Model","Basic Model", "Basic + State Model", "Demographic Model", "Final Model" ) ) + # specify order
      scale_y_continuous( labels = these.limits.character, # ensure labels on y-axis go to 2 digits
                          breaks = these.limits )+
      ylab( "Residual Variance" ) +
      geom_text( aes( label = sp.var.perc ), color = df.resid$color,
                 position=position_stack( vjust = 0.5 ),
                 size = 4, family = "Avenir" )
  }
  
  print( "Other Plots...END" )
  
  
  total.time.end <- Sys.time()
  
  print( total.time.end - total.time.start )
  
  return( list( model.fit = this.model, 
                model.summary = summary( this.model ),
                runs.required = runs,
                map.decomp = map.decomp,
                psi.plot = map.d$fixed.effect$psi.plot +
                  guides( fill = guide_colorbar( ticks.colour = NA,
                                                 frame.colour =  "black",
                                                 barwidth = 0.6,
                                                 barheight = 9 ) ),
                fixed.term = fixed.return,
                nnbs.var.plot = sp.var.nnbs.plot,
                explained.geo.var = df.exp.var,
                rand.eff.var.plot = if( !is.null( null.model.formula ) & !is.null( model.2.formula ) ) explained.var.plot else NULL,
                total.time.elapsed = total.time.end - total.time.start ) )
  
}

# recommended sizing for map decomposition plot
# ggsave( "04-Tables-Figures/map-decomposition.png",
#         width = 12, height = 9 )

# Examples:
# car_inla_analysis( d.jhu = d.poly.so, 
#                    d.cdc = d.poly.so,
#                    which.model = "cdc" )

# car_inla_analysis( d.jhu = d.poly, 
#                    d.cdc = d.poly,
#                    which.model = "jhu", selection.criterion = "dic", criterion.threshold = 5)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




####################################################################################################
##################################### Results Saving Function ###########################################
####################################################################################################
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


inla_results_save <- function( car.inla.object, path, tag, map.decomp.width = 4124,
                               map.decomp.height = 2808, dpi = 600 ){
  # car.inla.object - the object that results when running car_inla_analysis
  # path = a string with the path to the destination folder
  # tag = a string that will be appended to the file names (e.g., "-main-model")
  # map.decomp.width = width (in pixels) for decomposition plots. default is 4124
  # map.decomp.height = height (in pixels) for decomposition plots. default is 2808
  # the upper two dimensions may need to be adjusted until display of the plot is suitable to liking.
  
  ## load map decomposition ggplot object ##
  print("BEGIN...save map decomposition" )
  
  # save it with dimensions fixed
  ggsave( paste0( path, "map-decomposition-", tag, ".png" ),
          width = 22, height = 12 , dpi = dpi,
          plot = car.inla.object$map.decomp )
  
  print("DONE...save map decomposition" )
  
  
  ## fixed effects table ##
  print("BEGIN...fixed effects table" )
  write.table( car.inla.object$fixed.term, 
               paste0( path, "fixed-effect-", tag, ".txt" ), 
               sep = "," )
  print("DONE...fixed effects table" )
  
  
  ## psi plot ##
  print("BEGIN...save psi plot" )
  
  ggsave( paste0( path, "psi-plot-", tag, ".png" ), dpi = dpi,
          units = "px", width = (4124), height = (2808) , plot = car.inla.object$psi.plot )
  
  print("DONE...save psi plot" )
  
  
  ## spatial variance/neighbors plot ##
  print("BEGIN...save neighbors variance plot" )
  
  suppressWarnings( {
    # save it with dimensions fixed ##
    ggsave( paste0( path, "nbs-var-plot-", tag, ".png" ), dpi = dpi,
            width = 7, height = 6, plot = car.inla.object$nnbs.var.plot )
  }
  
  )
  print("DONE...neighbors variance plot" )
  
  
  ## proportion of spatial variance plot and table ##
  print("BEGIN...save spatial variance proportions plot and table" )
  
  # save it with dimensions fixed ##
  if( !is.null( car.inla.object$rand.eff.var.plot ) ){
    
    ggsave( paste0( path, "prop-var-plot-", tag, ".png" ), dpi = dpi,
            width = 8*0.8, height = 12*0.8, plot = car.inla.object$rand.eff.var.plot )
    
  }
  
  write.table( car.inla.object$explained.geo.var, 
               paste0( path, "rand-eff-explained-var-", tag, ".txt" ), 
               sep = "," )
  
  print("DONE...save spatial variance proportions plot and table" )
  
  if( !is.null( car.inla.object$selection.meta.data ) ){
    
    print("BEGIN...save selection meta data" )
    write.table( car.inla.object$selection.meta.data$log, 
                 paste0( path, "model-selection-log.txt" ), 
                 sep = "," )
    
    write.table( car.inla.object$selection.meta.data$add.meta.table, 
                 paste0( path, "model-selection-add-meta.txt" ), 
                 sep = "," )
    
    print("DONE...save selection meta data" )
  }
  
}

# example
# inla_results_save(res.jhu.1, path = paste0( "/Users/Chris/Downloads/results-function-jhu-try/" ),
#                   tag = "jhu-main")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




####################################################################################################
##################################### Significant Digits ###########################################
####################################################################################################
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# not my function: 
sigfig <- function(vec, digits){
  return(gsub("\\.$", "", formatC(signif(vec,digits=digits), digits=digits, format="fg", flag="#")))
}

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




#####################################################################################################
############################# Extract Components of the Model Formula ###############################
#####################################################################################################
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

model_terms <- function( formula ){
  # formula = an object of class "formula"
  
  # response variable
  resp <- as.character( formula )[2]
  
  # character vector of terms in the specified formula
  char.f <- as.character( formula )[3] %>%
    str_replace_all( ., '\\"', "'") %>%
    str_split( ., "\\s\\+\\s" ) %>%
    unlist() %>%
    str_remove_all( ., "\\s" ) # removes all white space to keep discrepancies from being flagged that are not actual discrepancies
  
  # offset terms to always be included
  offs <- char.f[ which( str_detect( char.f, "offset" ) ) ]
  
  # random effects terms to always be included
  reffs <- char.f[ which( str_detect( char.f, "f\\(" ) ) ]
  
  # fixed effects terms (we will loop through these)
  feffs <- char.f[ which( char.f %notin% c( reffs, offs ) ) ]
  
  return( list( response.var = resp,
                fixed.effects = feffs,
                random.effects = reffs,
                offset.term = offs ) )
}
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


###############################################################################################
############################# Add Letter Labels to Base R Plots ###############################
###############################################################################################
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# **NOTE: I did not write these functions; they are taken from:
# https://stackoverflow.com/questions/38439211/figure-labels-add-text-on-graphs-in-the-same-location-despite-figure-size

# first function (internal helper for `addfiglab`)
line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
  y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
  switch(side,
         `1` = grconvertY(-line * y_off, 'npc', 'user'),
         `2` = grconvertX(-line * x_off, 'npc', 'user'),
         `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
         `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
         stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}

# second function
addfiglab <- function(lab, xl = par()$mar[2], yl = par()$mar[3]) {
  
  text(x = line2user(xl, 2), y = line2user(yl, 3), 
       lab, xpd = NA, font = 2, cex = 1.5, adj = c(0, 1))
  
}

# example
# addfiglab("A")
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
