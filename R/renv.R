## Set-up the reproducible environment with the `renv` package ##

# RUN ONCE ONLY
# renv::init()

# if needing to take a snapshot of current project environment
renv::snapshot()

# to mirror the project environment used for this analysis
renv::restore()

# NOTE: INLA was installed from the developer website 
# Version: INLA [23.04.24]