##########################################################
##########################################################
### FUNCTIONS USED BY THE NBS POTENTIALS SCENARIO WORK ###
##########################################################

### NOTES: How to read the comment syntax of this script --
# A Five-Three-Five # pattern refers to the sections within each function (where relevant)
# A One-Three-One (# === #) pattern describe the R-environment input/output required/created by that function
# Two #s (i.e. ##) below the One-Three-One refers to the detail (typically column names) within the One-Three-One objects
# A single # just provides additional info about that data or object (typically valid options for that parameter)

#####################################################################
### Function to calculate proportional adoption of an NBS pathway ###
#####################################################################

### INPUTS:

# === # timeframe_df: dataframe format - requires the column 'year':
## 'year' - An integer vector that states the years that proportional adoptiong should be calculated for (e.g. 2000 - 2100)
# === # year_start: An integer stating in which year the NBS pathway reaches 0.1% gross enrollment of maximal extent (i.e. if the maximum extent is 1,000,000 ha, then the year in which 1,000 ha has been enrolled) (e.g. 2025)
# === # year_50: An integer stating in which year the NBS pathway reaches 50% gross enrollment of maximal extent (i.e. if the maximum extent is 1,000,000 ha, then the year in which 500,000 ha has been enrolled) (e.g. 2060)
# === # saturation: The number of years that the NBS pathway is deemed to continually generate emission reductions (e.g. once 'enrolled' a reforestation hectare is expected to generate additional ERs for 30 years before that hectare is unenrolled)

### OUTPUTS:

# === # Dataframe format with multiple columns that correspond to the years timeframe column that was input:
## 'enrol' - A numeric vector that states the absolute cumulative area that has been enrolled, expressed as a proportion of given extent (i.e. 0-1)
## 'unenrol' - A numeric vector that states the absolute cumulative area that has been unenrolled, expressed as a proportion of given extent (i.e. 0-1)
## 'net_enrol' - A numeric vector that states the difference between the area enrolled in year X minus the area unenrolled in year X, expressed as a proportion of given extent (i.e. 0-1 though it will never be close to 0 or 1)
## 'diff_enrol' - A numeric vector that states the annual change in gross enrollment from year to year, expressed as a proportion of given extent (i.e. 0-1) (e.g. if this is 0.1 in a year, then 10% of the available extent was enrolled that year)
## 'diff_unenrol' - A numeric vector that states the annual change in gross unenrollment from year to year, expressed as a proportion of given extent (i.e. 0-1) (e.g. if this is 0.1 in a year, then 10% of the available extent was unenrolled that year, regardless of how much was enrolled)
## 'diff_net' - A numeric vector that states the annual change in net enrollment from year to year, expressed as a proportion of given extent (i.e. 0-1) (e.g. if this is 0.1 in a year, then the difference between enrollment and unenrollment in that year was 10% of the available extent - i.e. it could be that 15% was enrolled but 5% was unenrolled)

proportional_adoption = function(timeframe_df,
                                 year_start,
                                 year_50,
                                 saturation) {
  enrol = 1/(1+exp(-((log(1/0.001-1)-log(1/0.5-1))/(year_50-year_start))*(timeframe_df$year-(log(1/0.001-1)/((log(1/0.001-1)-log(1/0.5-1))/(year_50-year_start))+year_start)))) # Calculate enrollment of each NBS
  unenrol = 1/(1+exp(-((log(1/0.001-1)-log(1/0.5-1))/(year_50-year_start))*((timeframe_df$year-saturation)-(log(1/0.001-1)/((log(1/0.001-1)-log(1/0.5-1))/(year_50-year_start))+year_start)))) # Calculate unenrollment of each NBS
  net_enrol = enrol - unenrol # Calculate the net change
  diff_enrol = c(0,diff(enrol)) # Calculate annual enrollment difference over time
  diff_unenrol = c(0,diff(unenrol)) # Calculate annual unenrollment difference over time
  diff_net = c(0,diff(net_enrol)) # Calculate annual net enrollment difference when accounting for unenrollment
  return(cbind(timeframe_df, enrol, unenrol, net_enrol, diff_enrol, diff_unenrol, diff_net)) # Output the resulting dataframe
} # End of function


##################################################################################################################
### Function to calculate 'staggered' Emission Reductions based on temporal profiles of the different pathways ###
##################################################################################################################

### INPUTS:

# === # area_timeframe_data: dataframe format - requires the columns 'year' and 'diff_net':
## 'year' - An integer vector that states the years that emission reductions should be calculated for (e.g. 2000 - 2100)
## 'area_Mha' - A numeric vector that states the annual net area enrolled (i.e. diff_net output from proportional_adoption function) to generate emission reductions for the aforementioned years (in hectares or equivalent)
# === # saturation: integer vector that states the number of years that a hectare is enrolled before being unenrolled - equates to 'number of years that emission reductions are generated by that hectare'
# === # function_type: categorical character vector type with 4 factor levels (linear, lognormal, constant, pulse) refers to temporal profile in annual emission reductions for that NBS type. Details on options below:
# linear = linear ramp up from 0 tCO2/ha/yr in year 0 to maximum tCO2/ha/yr in 'years_to_max' (see additional input parameter below) and then sustained for the remaining portion of the 'saturation timeline' (as defined by 'saturation' parameter above)
# lognormal = lognormal ramp up from 0 tCO2/ha/yr in year 0 to maximum tCO2/ha/yr in 'years_to_max' (see additional input parameter below) and then reduced down in accordance with lognormal distribution so that the area under the curve is consistent with the original data from maximum potentials (uses the 'saturation timeline' (as defined by 'saturation' parameter above) for the final year)
# constant = given annual ER intensity (tCO2/ha/yr; ter_intensity parameter below) is generated every year during entire 'saturation timeline'
# pulse = all emission reductions are generated immediately in year 1 but land must be protected for saturation timeline with 0 emissions generated each year after year 1 (note - should align with a complimentary ter_intensity designation)
# === # years_to_max: integer vector that states the number of years before the maximum emission reduction intensity is reached (used when function_type is linear, as described above)
# === # ter_intensity: numeric vector of annual emission reduction intensity (i.e. tCO2/ha/yr) of an NBS type that is used to convert the area (Ha) into total emission reductions (TERs)
# === # force_average: boolean vector to determine whether the total proportion should average to 1 or not (i.e. if proportional ER intensity is ramped up over 10 years then the final average will be < 1, which translates to fewer emission reductions over time.
# === # truncate_ERs: boolean vector to determine whether the total emission reductions (ERs) generated should be limited to the timeframe of the 'saturation horizon' or not. Only applies when a 'function_type' is lognormal
# === # prob_lognormal: numeric vector that must be > 0 and < 1 that defines the 'spread' of the lognormal distribution when a pathway is using that 'function_type'. Default is set to 0.75 but higher values will narrow the ERs generated around the year which reaches maximum sequestration and a lower value will broaden that distribution (longer tail)
# Setting this to TRUE (the default) forces the remaining years (after the ramp up) to compensate for this and increase remaining ERs/year)

### OUTPUTS:

# === # Dataframe format equivalent to 'area_timeframe_data' with one additional column called 'TERs_MtCO2'
## 'TERs_MtCO2' - A numeric vector that states the absolute emission reductions (tCO2e) generated by all active areas (i.e. excluding unenrolled area) for each year specified (e.g. 2000 - 2100)

annual_TER_calculator = function(area_timeframe_data,
                                 saturation,
                                 function_type,
                                 years_to_max,
                                 ter_intensity, # annual_ER_intensity in tCO2/ha/yr
                                 force_average = T,
                                 truncate_ERs = F,
                                 prob_lognormal = 0.75) {
  tmp = data.frame(time = 1:saturation) # Create simple dataframe of timeline from 1 to end of 'saturation horizon'
  if(function_type == "lognormal") { # Lognormal refers to the common distribution CO2 sequestration rates over time for trees/vegetation in the pathways where this applies - see the adoption input file
    E = qnorm(((2*prob_lognormal-1) + 1)/2)/sqrt(2) # Inverse error function using the probability distribution set
    Sigma = (-sqrt(2)*E+(2*E^2+4*log(saturation/years_to_max))^0.5)/2 # First parameter used to define the lognormal distribution (uses both saturation and years to max specific to this pathway)
    Mu=log(years_to_max)+Sigma^2 # Second parameter used to define the lognormal distribution (uses the Sigma and years to max specific to this pathway)
    total_ERs = ter_intensity*saturation # The area under the curve that the lognormal distribution needs to match (same as original 'constant' equivalent). Note the truncation comes from the below difference in calculation...
    if(truncate_ERs) {
      tmp$annualERs_tCO2_ha = total_ERs/prob_lognormal*dlnorm(tmp$time, Mu, Sigma) # When the calculation is truncated the last year (saturation) is still part of the tail but the total area under this truncated curve should match the constant equivalent
    } else {
      tmp = data.frame(time = 1:200) # When not truncating the area under the lognormal curve you include the entire area under the curve but with some of the ERs distributed well into the tail of the distribution. Therefore the time constraints is no longer informed by the saturation - note that the issue then becomes that the total simulation length may not be long enough to include all the ERs expected (e.g. if you just run 2020 to 2100)
      tmp$annualERs_tCO2_ha = total_ERs*dlnorm(tmp$time, Mu, Sigma) # When the calculation is not truncated the total_ERs are adjusted but otherwise the calculation is the same 
    }
  } else {
    if(function_type == "linear") { # Linear refers to a phased increase in emission reductions per hectare over time after enrollment
      if(force_average) {
        tmp$prop_annERs = ifelse((tmp$time/years_to_max)>1, 1+(((years_to_max+1)/2)-1)/(saturation-years_to_max), tmp$time/years_to_max)
      } else {
        tmp$prop_annERs = ifelse((tmp$time/years_to_max)>1, 1, tmp$time/years_to_max) # Linear regression to increase ERs up to maximum after X years, after max is reached it is sustained until unenrollment. NOTE - Could be a logistic function instead (e.g. (1/(1+exp(-0.4*annualERs$time))^5)) but need to link coefficients to 'years_to_max'
      }
    } else {
      if(function_type == "constant") { # Constant refers to a static number of emission reductions per hectare over time after enrollment
        tmp$prop_annERs = 1 # Fixed at 1 as a proportion of maximum annual ERs
      } else {
        if(function_type == "pulse") { # Pulse refers to those pathways where emission reductions are generated in the year they are avoided but not again afterwards
          tmp$prop_annERs = c(1, rep(0, length(tmp$time)-1)) # Proportion starts at 1 but is set to 0 for every year after that. Assumption is that area needs to remain protected.
        } else {
          cat("The 'function_type' that was input is not recognized. Current support for 'linear, constant or pulse only.") # In case there is an error thrown, print a message to console
        } # End of error statement
      } # End of pulse calculation 
    } # End of constant calculation 
    tmp$annualERs_tCO2_ha = tmp$prop_annERs*ter_intensity # Convert proportion to total using intensity
  } # End of all if statements
  annual_ERs = area_timeframe_data # Create replicated dataframe to put calculated data into
  annual_ERs$TERs_MtCO2 = NA # Empty column to fill
  for(yr in area_timeframe_data$year) { # Loop over all years in time range
    TERs = 0 # Reset TERs object to 0
    for(i in tmp$time) { # Loop over years of annual ERs temporal profile
      if(yr+1-i>=min(area_timeframe_data$year)) { # Needed to stop trying to calculate a year below the minimum of the time range (e.g. 2000)
        TERs = TERs + area_timeframe_data$area_Mha[area_timeframe_data$year==(yr+1-i)]*tmp$annualERs_tCO2_ha[tmp$time==i] # Cumulate emission reductions for that year (new area plus previously enrolled areas)
      } # End of if statement
    } # End of annual ER calculation loop
    annual_ERs$TERs_MtCO2[annual_ERs$year==yr] = TERs # Allocate the cumulated TERs to the created dataframe
  } # End of full time-frame loop
  return(annual_ERs) # Output the resulting calculation
} # End of function


############################################################################################################################################
### Functionalized version of the whole S-curve approach (uses functions above) resulting in annual emission reductions for all pathways ###
############################################################################################################################################

### INPUTS:

# === # run_name: Character vector used to give this run a unique identifier for further analysis
# Recommended naming convention: "Geographic scale - NBS pathways included - Adoption parameter scenario - Maximal extent scenario - DD/MM/YYYY of run"
# === # inputs: Dataframe format - requires the columns 'NBS_short_name', 'year_start', 'year_50', 'saturation', 'function_type', 'years_to_max''NBS_short_name', 'extent_rate', 'extent' and 'ter_intensity'
# - Most of the columns are used by the two embedded functions (proportional_adoption and annual_TER_calculator) but in summary:
## 'NBS_short_name' - A character vector with no spaces that gives short names of all the NBS pathways in this simulation
## year_start - Integer stating in which year the NBS pathway reaches 0.1% gross enrollment of maximal extent (i.e. if the maximum extent is 1,000,000 ha, then the year in which 1,000 ha has been enrolled) (e.g. 2025)
## year_50 - Integer stating in which year the NBS pathway reaches 50% gross enrollment of maximal extent (i.e. if the maximum extent is 1,000,000 ha, then the year in which 500,000 ha has been enrolled) (e.g. 2060)
## saturation - Integer vector that states the number of years that a hectare is enrolled before being unenrolled - equates to 'number of years that emission reductions are generated by that hectare'
## function_type - Categorical character vector type with 3 factor levels (linear, constant, pulse) refers to temporal profile in annual emission reductions for that NBS type. Details on options below:
# linear = linear ramp up from 0 tCO2/ha/yr in year 0 to maximum tCO2/ha/yr in 'years_to_max' (see additional input parameter below) and then sustained for the remaining portion of the 'saturation timeline' (as defined by 'saturation' parameter above)
# lognormal = lognormal ramp up from 0 tCO2/ha/yr in year 0 to maximum tCO2/ha/yr in 'years_to_max' (see additional input parameter below) and then reduced down in accordance with lognormal distribution so that the area under the curve is consistent with the original data from maximum potentials (uses the 'saturation timeline' (as defined by 'saturation' parameter above) for the final year)
# constant = given annual ER intensity (tCO2/ha/yr; ter_intensity parameter below) is generated every year during entire 'saturation timeline'
# pulse = all emission reductions are generated immediately in year 1 but land must be protected for saturation timeline with 0 emissions generated each year after year 1 (note - should align with a complimentary ter_intensity designation)
## years_to_max - Integer vector that states the number of years before the maximum emission reduction intensity is reached (used when function_type is linear, as described above)
## extent_rate - Binary integer vector (0 or 1) to indicate whether extent data (i.e. the extent column values) for that NBS pathway are provided in Mha or Mha/yr
## extent - Numeric vector reporting the maximal extent for that pathway reported in total Millions of hectares (Mha) or Mha per year (i.e. 678 Mha maximum available for reforestation but 8.97 Mha/yr maximum for avoiding deforestation)
## ter_intensity - Numeric vector of average annual emission reduction intensity (i.e. tCO2/ha/yr) of an NBS type that is used to convert the area (Ha) into total emission reductions (TERs)
# === # timeframe_start: Integer value to represent the first year that you want to run the simulation over (defaults to year 2000)
# === # timeframe_end: Integer value to represent the last year that you want to run the simulation over (defaults to year 2200)
# === # force_average_ERprofile: Boolean vector to determine whether the total proportion should average to 1 or not (i.e. if proportional ER intensity is ramped up over 10 years then the final average will be < 1, which translates to fewer emission reductions over time.
# Setting this to TRUE (the default) forces the remaining years (after the ramp up) to compensate for this and increase remaining ERs/year)
# === # output_area_data: Boolean vector to state whether the function should also output the entire adoption data to the R environment (in addition to the standard annual TER data returned)
# Setting this to TRUE (the default) means that all forms (i.e. area enrolled each year, area unenrolled, cumulative, etc. - see 'proportional_adoption function for all) of the area/adoption/enrolment data is also written to the R environment
# Note that the adoption/area data used to calculate the annual TERs are included in the standard function output, this extra option just provides more info

### OUTPUTS:

# === # Dataframe format that includes the annual Total Emission Reductions (TERs) and areas generating them, for each NBS pathway simulated over the whole timeframe (i.e. between 2000 and 2200). Includes the columns:
## 'year' - Integer vector showing the years over which the simulation was run
## 'NBS_short_name' - Character vector with a short name for each NBS pathway (no spaces) -- tied to the inputs provided by the inputs dataframe 
## 'adoption_datatype' - A character string vector that states which 'adoption data' was used to calculate the annual emission reductions. Relates to the proportional_adoption function output columns
# Note that 'extent_rate' NBS pathways (i.e. those with maximal extent expressed as Mha/yr) use 'net change in cumulative area to date' to calculate emission reductions, but others just use absolute annual gross enrolment.
## 'area_proportion' - A numeric vector showing change in enrollment from year to year, expressed as a proportion of given extent (i.e. 0-1) (e.g. if this is 0.1 in a year, then 10% of the available extent was enrolled that year)
## 'area_Mha' - A numeric vector that converts the 'area_proportion' to an absolute ADDITIONAL area (in Mha) that is enrolled to generate emission reductions that year  
## 'TERs_MtCO2' - A numeric vector that states the amount of emission reductions generated in that each year in million tonnes CO2e per year
## 'area_Mha_todate' - A numeric vector that cumulates the area_Mha column to show the total area of each NBS pathway that has been enrolled to date (doesn't account for any unenrolment -- tied to maximal extent, not emissions being generated)
## 'TERs_MtCO2_todate' - A numeric vector that cumulates the TERs_MtCO2 column to show the total emission reductions of each NBS pathway that has been generated to date
## 'run_ID' - Simple character vector with the provided "run_name" used by the function
# === # (OPTIONAL) Dataframe format that shows the full annual and cumulative enrolment and unenrolment (and net change) areas for each NBS pathway for every year over the simulated timeframe. Includes teh columns:
# Same columns as above dataframe output, see descriptions. Only difference is the 'adoption_datatype' factor levels. These are:
# Area enrolled - Cumulative 'area' enrolled to date
# Area unenrolled" - Cumulative 'area' unenrolled to date
# Net change in area - Difference between cumulative 'area' enrolled to date and cumulative 'area' unenrolled to date
# Annual gross enrollment - Annual change in 'area' enrolled
# Annual gross unenrollment - Annual change in 'area' unenrolled
# Annual net enrollment - Net annual change in 'area' between enrolled and unenrolled

Annual_NBS_potentials = function(run_name = NULL,
                                 inputs,
                                 timeframe_start = 2000,
                                 timeframe_end = 2200,
                                 force_average_ERprofile = T,
                                 output_area_data = F,
                                 print_messages = T) {
  #####
  ### SETUP AND INITIAL PARAMETER CHECKS
  #####
  
  # If the run hasn't been given a unique name, throw an error and report as much
  if(!is.character(run_name)) { 
    cat("RUN FAILED! \n")
    stop("Specify a valid character name for this run - use the run_name input") # If NULL, stop the function
  } 
  
  # Requires packages - could potentially be stripped back to base but not worth it at this point
  if(!require(tidyverse)) install.packages("tidyverse")
  library(tidyverse) # Needed for simple data-wrangling so that outputs are simple and easy to understand (mainly column headers and factor levels)
  
  # Get the full list of NBS pathways that you want to analyse (should be complete in input file)
  TNCfull_NBSlist = as.character(inputs[["NBS_short_name"]]) 
  
  #####
  ### CALCULATING ADOPTION RATES
  #####
  
  # Loop over pathways to calculate adoption for each one using the provided input data for S-curves
  adoption_perc = data.frame() # Create empty dataframe to put the adoption data into
  
  for(i in TNCfull_NBSlist) { # Start of adoption calculation loop
    timeframe_df = data.frame(year=timeframe_start:timeframe_end) # Create simple dataframe with column of continuous years running from start to finish
    tmp_df = inputs[inputs$NBS_short_name==i,] # Limit adoption input data to the pathway of interest (in the loop)
    
    # Use proportional_adoption function (user-defined) to calculate adoption rates for the provided years
    tmp = proportional_adoption(timeframe_df = timeframe_df, 
                                year_start = tmp_df$year_start,
                                year_50 = tmp_df$year_50,
                                saturation = tmp_df$saturation)
    
    tmp$NBS_short_name = i # Assign the pathway name to a column
    adoption_perc = rbind(adoption_perc, tmp) # Bind all pathway adoption data together
  } # End of adoption calculation loop
  
  # To allow for simple plotting, convert dataframe from wide to long format and rename for plot ease/aesthetics
  adoption_perc_alldata = adoption_perc %>%
    rename("Area enrolled" = "enrol", # Cumulative 'area' enrolled to date
           "Area unenrolled" = "unenrol", # Cumulative 'area' unenrolled to date
           "Net change in area" = "net_enrol", # Difference between cumulative 'area' enrolled to date and cumulative 'area' unenrolled to date
           "Annual gross enrollment" = "diff_enrol", # Annual change in 'area' enrolled
           "Annual gross unenrollment" = "diff_unenrol", # Annual change in 'area' unenrolled
           "Annual net enrollment" = "diff_net") %>% # Net annual change in 'area' between enrolled and unenrolled)
    gather(variable, value, -year, -NBS_short_name) # Create long-form, plot-friendly dataframe
  
  # Rename factor levels
  colnames(adoption_perc_alldata)[colnames(adoption_perc_alldata) %in% c("variable")] = "adoption_datatype" # Rename the variable column
  colnames(adoption_perc_alldata)[colnames(adoption_perc_alldata) %in% c("value")] = "area_proportion" # Rename the value column
  
  #####
  ### CONVERT ADOPTION PERCENT TO HECTARES FOR MAXIMUM EXTENT
  #####
  
  ### Convert the percentage adoption into hectares
  adoption_alldata = data.frame() # Create empty dataframe to put data into
  for(i in TNCfull_NBSlist) { # Loop over each NBS and allocate based on maximum hectare extent
    tmp = adoption_perc_alldata[adoption_perc_alldata$NBS_short_name == i,] # Area proportion timeline - Subset based on pathway
    tmp$area_Mha = tmp$area_proportion * inputs$extent[inputs$NBS_short_name == i] # Simply multiply proportions by maximum extent
    adoption_alldata = rbind(adoption_alldata, tmp) # Combine new dataframes
  }
  
  #####
  ### CONVERT HECTARES TO EMISSIONS REDUCTIONS
  #####
  
  ### As above, but convert the hectares into maximal TERs
  TERs_alldata = data.frame() # Create empty dataframe to fill with Total Emission Reduction data
  for(i in TNCfull_NBSlist) { # Loop over each NBS pathway
    tmp = adoption_alldata[adoption_alldata$NBS_short_name == i,] # Area absolute timeline - Subset based on pathway
    tmp2 = inputs[inputs$NBS_short_name == i,] # Inputs needed for TER calculator function - Subset based on pathway
    ### Importantly, the global extent for certain avoidance NBS pathways is expressed as a rate (i.e. Mha/yr) instead of absolute total (i.e. Mha) so annual ER calculations must use the 'cumulative' net change in area instead of 'annual' net change
    if(inputs$extent_rate[inputs$NBS_short_name==i] == 1) { # If the pathway expresses extent as a rate (i.e. Mha/yr) then use 'cumulative' area data
      tmp = tmp[tmp$adoption_datatype=="Net change in area",] # Subset to only include datatype that is relevant - NOTE: 'Net change in area' for those pathways with extents expressed as a rate actually means the gross new enrollment that year
      
      # Use annual_TER_calculator function (user-defined) to calculate lag of emission reductions specific to that pathway over time - uses inputs from inputs
      tmp_annual_ERs = annual_TER_calculator(area_timeframe_data = tmp, # Proportion enrollment rate - note, uses the 'diff_net' column as this is the annual proportional change (i.e. enrolled that year - unenrolled that year)
                                             saturation = tmp2$saturation, # 'Saturation' is actually the length of time that credits are generated or hectare must be maintained to generate ERs
                                             function_type = tmp2$function_type, # Function type refers to the temporal profile of how emission reductions are generated for that NBS pathway
                                             years_to_max = tmp2$years_to_max, # Years to max is used by linear function types to calculate the slope of 'ramp up' in ERs (i.e. how long it takes to reach maximum ERs generated per year)
                                             ter_intensity = inputs$ter_intensity[inputs$NBS_short_name==i], # Maximum emission reduction intensity used to scale annual ERs
                                             force_average = force_average_ERprofile) # Use previously set input as to whether you want to force the calculated overall average ER intensity to equal the value in inputs (default to TRUE)
      
      # Because the TER calculations are based off the number of new hectares that year, pathways are expressed as a rate. As a result it's also useful to have a cumulative column showing total number of hectares enrolled to date (i.e. have Mha as well as Mha/yr)
      tmp_annual_ERs$area_Mha_todate = cumsum(tmp_annual_ERs$area_Mha) # Cumulate the areas to that date - note that when area is unenrolled the actually means less and less is protected (assumes after 'saturation' - i.e. 100 years then that area is no longer at threat)
      tmp_annual_ERs$TERs_MtCO2_todate = cumsum(tmp_annual_ERs$TERs_MtCO2) # Cumulate the emission reductions generated to that date
    } else { # The pathway expresses extent as a global absolute value (i.e. Mha) then use 'annual' area data
      tmp = tmp[tmp$adoption_datatype=="Annual gross enrollment",] # Subset to only include datatype that is relevant - NOTE: For these we use 'gross' enrollment instead of 'net' because the calculator uses the same 'saturation' value to limit the years which generate ERs (and therefore we don't need to account for unenrollment as well)
      
      # Use annual_TER_calculator function (user-defined) to calculate lag of emission reductions specific to that pathway over time - uses inputs from inputs
      tmp_annual_ERs = annual_TER_calculator(area_timeframe_data = tmp, # Proportion enrollment rate - note, uses the 'diff_net' column as this is the annual proportional change (i.e. enrolled that year - unenrolled that year)
                                             saturation = tmp2$saturation, # 'Saturation' is actually the length of time that credits are generated or hectare must be maintained to generate ERs
                                             function_type = tmp2$function_type, # Function type refers to the temporal profile of how emission reductions are generated for that NBS pathway
                                             years_to_max = tmp2$years_to_max, # Years to max is used by linear function types to calculate the slope of 'ramp up' in ERs (i.e. how long it takes to reach maximum ERs generated per year)
                                             ter_intensity = inputs$ter_intensity[inputs$NBS_short_name==i], # Maximum emission reduction intensity used to scale annual ERs
                                             force_average = force_average_ERprofile) # Use previously set input as to whether you want to force the calculated overall average ER intensity to equal the value in inputs (default to TRUE)
      
      tmp_annual_ERs$area_Mha_todate = cumsum(tmp_annual_ERs$area_Mha) # Cumulate the areas to that date so this 'area' column is consistent with 'avoidance' pathways. Assumption is that once a hectare is enrolled it must be protected indefinitely
      tmp_annual_ERs$TERs_MtCO2_todate = cumsum(tmp_annual_ERs$TERs_MtCO2) # Cumulate the emission reductions generated to that date
    } # End of if statement
    TERs_alldata = rbind(TERs_alldata, tmp_annual_ERs) # Combine all pathways
  } # End of loop
  
  # Add the given run name to the output 
  TERs_alldata$run_ID = run_name
  
  ##### 
  # IMPORTANT DATA CHANGE - ADOPTION AREA DATA FOR 'EXTENT RATE' PATHWAYS NEED CHANGING TO ABSOLUTE AREA (from Mha/yr)
  #####
  
  tmp_output = data.frame() # Empty dataframe to put data from loop into
  for(NBS in TNCfull_NBSlist) { # Loop over each NBS pathway
    if(inputs$extent_rate[inputs$NBS_short_name==NBS]==0) { # Use input parameter of 'extent_rate' to determine whether anything needs changing
      tmpDF = adoption_alldata # Work with duplicate dataset to avoid overwriting issues
      tmpDF = tmpDF[tmpDF$NBS_short_name == NBS,] # If the NBS pathway does not have its potential expressed as a rate (i.e. absolute Mha used instead of Mha/yr) then just use the existing data
    } else { # If NBS pathway IS an extent_rate pathway (i.e. potential is expressed as Mha/yr) then:
      tmpDFx = adoption_alldata # Work with duplicate dataset to avoid overwriting issues
      tmpDFx = tmpDFx[tmpDFx$NBS_short_name == NBS,] # Subset based on the pathway
      tmpDF = data.frame() # Create empty dataframe to put data from loop into
      for(i in levels(as.factor(tmpDFx$adoption_datatype))) { # Loop over each adoption datatype to cumulate for each 
        tmpDFx2 = tmpDFx[tmpDFx$adoption_datatype == i,] # Subset based on adoption_datatype
        tmpDFx2$area_Mha = cumsum(tmpDFx2$area_Mha) # Cumulate the Mha/yr values so they now actually refer to the Mha absolute enrolled to-date
        tmpDF = rbind(tmpDF, tmpDFx2) # Bind with the 'non' extent_rate pathways
      } # End of loop
    } # End of else statement
    tmp_output = rbind(tmp_output, tmpDF) # Bind each NBS pathway to the larger dataframe
  } # End of loop
  
  # If properly calculated in the above loop the new dataframes should replace the existing adoption dataframes
  # NOTE - because of this, TERs cannot be calculated in the same way as before using these new dataframes
  adoption_alldata = tmp_output
  
  # Add the given run name to the output 
  adoption_alldata$run_ID = run_name
  
  ##### 
  # FINAL OUTPUTS AND MESSAGES
  #####
  
  if(print_messages) cat(paste0("- RUN COMPLETED SUCCESSFULLY \n    Run Name: '", run_name, "'\n\n"))
  if(output_area_data){ # If that section was run, save to environment
    if(print_messages) cat("- ADOPTION DATASET SAVED TO THE GLOBAL ENVIRONMENT IN R \n  Important note: \n    All values in this adoption dataset are in Millions of hectares (Mha). But, some have been converted from Mha/yr. \n    As a result, this 'converted' dataset cannot be used to back-calculate the main TER output dataset.\n\n")
    assign(paste0("adoption_alldata_", run_name), adoption_alldata, envir = .GlobalEnv) # Output the resulting calculation
  } # End of if statement
  return(TERs_alldata) # Output the dataset from this analysis
  
} # End of function
