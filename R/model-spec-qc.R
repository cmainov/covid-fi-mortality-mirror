library( dplyr )
library( stringr )
# model specs
source( "R/model-specs.R" )

# helper functions
source( "R/utils.R" )

## Model 5 (Final Model w/ Poverty) ##
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

jhu.5.pov <- model_terms( f.jhu.model.5.pov )
cdc.5.pov <- model_terms( f.cdc.model.5.pov )

# fixed effects
jhu.5.pov.f <- jhu.5.pov$fixed.effects
cdc.5.pov.f <- cdc.5.pov$fixed.effects

jhu.5.pov.f[ jhu.5.pov.f %notin% cdc.5.pov.f ]
cdc.5.pov.f[ cdc.5.pov.f %notin% jhu.5.pov.f ]

# random effects
jhu.5.pov.r <- jhu.5.pov$random.effects
cdc.5.pov.r <- cdc.5.pov$random.effects

jhu.5.pov.r[ jhu.5.pov.r %notin% cdc.5.pov.r ]
cdc.5.pov.r[ cdc.5.pov.r %notin% jhu.5.pov.r ]

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## Final Model w/o Poverty ##
## Model 5 (Final Model) ##
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

jhu.5 <- model_terms( f.jhu.model.5 )
cdc.5 <- model_terms( f.cdc.model.5 )

# fixed effects
jhu.5.f <- jhu.5$fixed.effects
cdc.5.f <- cdc.5$fixed.effects

jhu.5.f[ jhu.5.f %notin% cdc.5.f ]
cdc.5.f[ cdc.5.f %notin% jhu.5.f ]

# random effects
jhu.5.r <- jhu.5$random.effects
cdc.5.r <- cdc.5$random.effects

jhu.5.r[ jhu.5.r %notin% cdc.5.r ]
cdc.5.r[ cdc.5.r %notin% jhu.5.r ]

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## Final Model w/ Poverty Compared to w/o Poverty ##
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# JHU
# fixed effects
jhu.5.f[ jhu.5.f %notin% jhu.5.pov.f ]
jhu.5.pov.f[ jhu.5.pov.f %notin% jhu.5.f ]

# random effects
jhu.5.r[ jhu.5.r %notin% jhu.5.pov.r ]
jhu.5.pov.r[ jhu.5.pov.r %notin% jhu.5.r ]

# CDC
# fixed effects
cdc.5.f[ cdc.5.f %notin% cdc.5.pov.f ]
cdc.5.pov.f[ cdc.5.pov.f %notin% cdc.5.f ]

# random effects
cdc.5.r[ cdc.5.r %notin% cdc.5.pov.r ]
cdc.5.pov.r[ cdc.5.pov.r %notin% cdc.5.r ]

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## Model 4 (Demographic Model) ##
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

jhu.4 <- model_terms( f.jhu.model.4 )
cdc.4 <- model_terms( f.cdc.model.4 )

# fixed effects
jhu.4.f <- jhu.4$fixed.effects
cdc.4.f <- cdc.4$fixed.effects

jhu.4.f[ jhu.4.f %notin% cdc.4.f ]
cdc.4.f[ cdc.4.f %notin% jhu.4.f ]

# random effects
jhu.4.r <- jhu.4$random.effects
cdc.4.r <- cdc.4$random.effects

jhu.4.r[ jhu.4.r %notin% cdc.4.r ]
cdc.4.r[ cdc.4.r %notin% jhu.4.r ]

# ---------------------------------------------------------------------------------------------------------------------------------------------------------


## Model 3 (Basic + State Model) ##
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

jhu.3 <- model_terms( f.jhu.model.3 )
cdc.3 <- model_terms( f.cdc.model.3 )

# fixed effects
jhu.3.f <- jhu.3$fixed.effects
cdc.3.f <- cdc.3$fixed.effects

jhu.3.f[ jhu.3.f %notin% cdc.3.f ]
cdc.3.f[ cdc.3.f %notin% jhu.3.f ]

# random effects
jhu.3.r <- jhu.3$random.effects
cdc.3.r <- cdc.3$random.effects

jhu.3.r[ jhu.3.r %notin% cdc.3.r ]
cdc.3.r[ cdc.3.r %notin% jhu.3.r ]

# ---------------------------------------------------------------------------------------------------------------------------------------------------------


## Model 2 (Basic Model) ## #check this
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

jhu.2 <- model_terms( f.jhu.model.2 )
cdc.2 <- model_terms( f.cdc.model.2 )

# fixed effects
jhu.2.f <- jhu.2$fixed.effects
cdc.2.f <- cdc.2$fixed.effects

jhu.2.f[ jhu.2.f %notin% cdc.2.f ]
cdc.2.f[ cdc.2.f %notin% jhu.2.f ]

# random effects
jhu.2.r <- jhu.2$random.effects
cdc.2.r <- cdc.2$random.effects

jhu.2.r[ jhu.2.r %notin% cdc.2.r ]
cdc.2.r[ cdc.2.r %notin% jhu.2.r ]

# ---------------------------------------------------------------------------------------------------------------------------------------------------------


## Null Model ##
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

jhu.1 <- model_terms( f.jhu.model.null )
cdc.1 <- model_terms( f.cdc.model.null )

# fixed effects
jhu.1.f <- jhu.1$fixed.effects
cdc.1.f <- cdc.1$fixed.effects

jhu.1.f[ jhu.1.f %notin% cdc.1.f ]
cdc.1.f[ cdc.1.f %notin% jhu.1.f ]

# random effects
jhu.1.r <- jhu.1$random.effects
cdc.1.r <- cdc.1$random.effects

jhu.1.r[ jhu.1.r %notin% cdc.1.r ]
cdc.1.r[ cdc.1.r %notin% jhu.1.r ]

# ---------------------------------------------------------------------------------------------------------------------------------------------------------


## Compare Null to Basic Model ##
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# fixed effects
jhu.2.f[ jhu.2.f %notin% jhu.1.f ]
jhu.1.f[ jhu.1.f %notin% jhu.2.f ] # should be 0 

# random effects
jhu.2.r == jhu.1.r

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## Compare Basic to Basic + State Model ##
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# fixed effects
jhu.3.f[ jhu.3.f %notin% jhu.2.f ]
jhu.2.f[ jhu.2.f %notin% jhu.3.f ] # should be 0 

# random effects
jhu.3.r == jhu.2.r

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## Compare Demographic to Basic + State Model ##
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# fixed effects
jhu.4.f[ jhu.4.f %notin% jhu.3.f ]
jhu.3.f[ jhu.3.f %notin% jhu.4.f ] # should be 0 

# random effects
jhu.4.r == jhu.3.r

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## Compare Final to Demographic Model ##
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# fixed effects
jhu.5.f[ jhu.5.f %notin% jhu.4.f ]
jhu.4.f[ jhu.4.f %notin% jhu.5.f ] # should be 0 

# random effects
jhu.5.r == jhu.4.r

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
