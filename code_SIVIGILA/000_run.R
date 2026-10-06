###-----------------------------------------------------------------------------
## Before running the entire pipeline with these two simple commands, verify that 
## you are working in the same directory (same folder); you can check with getwd().
## It is also important to have the raw data of the 3 databases and the 
## DANE file for the maps.

## Important: Before running the script, you must manually download the raw data 
## from the Dengue (210), Severe Dengue (220), and Dengue Mortality (580) databases 
## using the official SIVIGILA portal (https://portalsivigila.ins.gov.co/).

##------------------------------------------------------------------------------

source("run_all.R")
correr_pipeline()
