###############################################################
###############################################################
### CALCULATING ANNUAL EMISSION REDUCTIONS FOR NBS PATHWAYS ###
###  --- ANALYSIS FOR MANUSCRIPT: Robertson et al 2021 ---  ###
###############################################################
# SIMPLE MODEL-ONLY VERSION AS EXAMPLE #
########################################

rm(list=ls()) # Start with a blank slate if needed

# Libraries
library(tidyverse)

# General user-inputs
use_all_pathway_subtypes = F # Do you want to analyse all pathway subtypes separately? - otherwise group them? (T=Yes, F=No)
cost_threshold = "Max" # Which do you want to use a maximum biophysical potential, a $100/tCO2 threshold or $10/tCO2 threshold? - "Max", "100" or "10
use_average_ERintensities = T # Do you want to force the proportional averages to equal 1 overall (T=Yes, F=No)
starttime = 2000 # Default setting to start the scenario from 2000 but could choose current year if you wish
endtime = 2200 # By default the scenario runs from year 2000 but this is a user-defined input to assign the ending year (i.e. 2100 or 2200, or whatever!)
# Note - Of course end year must be after the start year for the scenario to work

#####
### INPUTS AND SET UP
#####
TopDir = getwd() # Assumes project is run from same location that this script is stored. This is denoted ShinyDir and contains two folders (data and figures) as well as this script
setwd(TopDir) # Set the working directory to be explicit
DataDir = file.path(TopDir, "data") # Essentially the input folder (where input data is stored)

# Load functions
source(file.path(TopDir, "NBS_ET_functions.R")) # Load functions used to calculate proportional adoption over time and derive ERs over time

# Load necessary datafiles - Note: newer files have uncertainty columns as well
adopt_in = read.csv(file.path(DataDir, "adoption_uncertALL.csv"), header=T) # Read in the inputs required to calculate annual adoption rates of each NBS type
extent_in = read.csv(file.path(DataDir, "extent_uncertALL.csv"), header=T) # Read in the inputs required to convert adoption rates to absolute values of hectares and emissions reductions

# Choose whether to do the full breakdown of all activities or the grouped 23 pathways
if(!use_all_pathway_subtypes){
  adopt_in = adopt_in %>%
    group_by(NBS_subgroup) %>%
    summarise(NBS_desc = head(NBS_subgroup,1),
              NBS_short_name = head(NBS_short_name,1),
              time_at_threat = max(time_at_threat),
              function_type = head(function_type,1),
              years_to_max = max(years_to_max),
              year_start = max(year_start),
              year_50 = max(year_50),
              saturation = max(saturation))
  extent_in = extent_in %>%
    group_by(NBS_subgroup) %>%
    summarise(NBS_desc = head(NBS_subgroup,1),
              NBS_short_name = head(NBS_short_name,1),
              extent_rate = max(extent_rate),
              extent = sum(extent),
              global_extent = sum(global_extent),
              global_ter = sum(global_ter),
              ters_100 = sum(ters_100),
              ters_10 = sum(ters_10),
              restore_type = head(restore_type,1),
              active_remove = max(active_remove),
              prop_co2 = max(prop_co2),
              prop_ch4 = max(prop_ch4),
              prop_n2o = max(prop_n2o))
  extent_in$ter_intensity = extent_in$global_ter/extent_in$global_extent
}

# Choose which cost threshold to apply to the simulation/analysis
if(cost_threshold == "100") {
  extent_in$extent = extent_in$ters_100/extent_in$global_ter*extent_in$global_extent
} else {
  if(cost_threshold == "10") {
    extent_in$extent = extent_in$ters_10/extent_in$global_ter*extent_in$global_extent
  }
}

#####
### RUN BASIC MODEL AND PLOT
#####

# Combine the inputs specific to that scenario
NBS_in = merge(adopt_in, extent_in) # Create one dataframe that can be used as inputs to the function

# Run the model and write the output to the global environment
results = Annual_NBS_potentials(run_name = "Example simulation",
                                inputs = NBS_in,
                                timeframe_start = starttime,
                                timeframe_end = endtime,
                                output_area_data = F,
                                force_average_ERprofile = use_average_ERintensities)

# Quick plotcheck
ggplot(merge(results, extent_in, by="NBS_short_name"), aes(x=year, y=TERs_MtCO2/1000, fill=NBS_desc)) +
  geom_area(colour='black') +
  ylab(expression(paste("Global annual emission reductions (Gt ", CO[2], "e ", yr^-1, ")"))) +
  theme_bw() + 
  theme(legend.title = element_blank(),
        legend.position = "bottom")

#####
### CLEAN UP AND SAVE AS NECESSARY
#####

# Remove all unnecessary objects
rm(list=ls(pattern = "^tmp"))
