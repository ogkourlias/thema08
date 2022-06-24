########## BASE PARAMETERS ##########

# Defining parameters that stay the same
base_parameters <- c(
  # Parameters governing cell growth
  rH = 0.00275, #The per-capita growth rate of normal cells. Unit: hr^-1
  rC = 0.003, #	The per-capita growth rate of tumor cells. Unit: hr^−1
  KH = 10^11, #		The carrying capacity of normal cells. Unit: cells
  KC = 1.47*10^12 #		The carrying capacity of tumor cells. Unit: cells
)

########## ADENO VIRUS ##########

# Defining adeno virus specific parameters for tumour cells
adeno_parms <- append(base_parameters, c(
  lambdaC = 1/48, # Lysing rate of tumor cells. Unit: cell^-1hr^-1
  bC = 1000, # Burst size from lysing infected tumor cells. Unit: -
  omega = 1, # Viral clearance rate. Unit: virus^-1hr^-1
  betaC = 5*10^-12 # Uptake/encounter/infection rate of tumor cells. Unit: viruses cell^−1hr^−1
))

# Appending adeno virus parameters with values for normal cells
adeno_parms <- append(adeno_parms, c(
  bH = (adeno_parms["bC"] * 0.1), # Burst size of normal (Healthy) cells. Unit: -
  lambdaH = (adeno_parms["lambdaC"] * 0.1), # Lysing rate of normal (Healthy) cells
  betaH = (adeno_parms["betaC"] * 0.0025) #(4*10^4) # Infection rate of normal cells. Normal cell to tumor cell infection rate ratio is 0.0025
))

# Replace weird "bH.bC" type names with just "bH"
names(adeno_parms) <- gsub("\\..*", "", names(adeno_parms))

########## HSV ##########

# Defining hsv specific parameters for tumour cells
hsv_parms <- append(base_parameters, c(
  lambdaC = 1/18, # Lysing rate of tumor cells. Unit: cell^-1hr^-1
  bC = 50, # Burst size from lysing infected tumor cells. Unit: -
  omega = 0.025, # Viral clearance rate. Unit: virus^-1hr^-1
  betaC = 2.5*10^-12 # Uptake/encounter/infection rate of tumor cells. Unit: viruses cell^−1hr^−1
))

# Appending hsv parameters with values for normal cells
hsv_parms <- append(hsv_parms, c(
  bH = (hsv_parms["bC"] * 0.1), # Burst size of normal (Healthy) cells. Unit : -
  lambdaH = (hsv_parms["lambdaC"] * 0.1), # Lysing rate of normal (Healthy) cells
  betaH = (hsv_parms["betaC"] * 0.0025) # Infection rate of normal cells. Normal cell to tumor cell infection rate ratio is 0.0025
))

# Replace weird "bH.bC" type names with just "bH"
names(hsv_parms) <- gsub("\\..*", "", names(hsv_parms))

########## VSV ##########

# Defining vsv specific parameters for tumour cells
vsv_parms <- append(base_parameters, c(
  lambdaC = 1/24, # Lysing rate of tumor cells. Unit: cell^-1hr^-1
  bC = 1350, # Burst size from lysing infected tumor cells. Unit: -
  omega = 0.244, # Viral clearance rate. Unit: virus^-1hr^-1
  betaC = 5*10^-13 # Uptake/encounter/infection rate of tumor cells. Unit: viruses cell^−1hr^−1
))

# Appending vsv parameters with values for normal cells
vsv_parms <- append(vsv_parms, c(
  bH = (vsv_parms["bC"] * 0.1), # Burst size of normal (Healthy) cells. Unit : -
  lambdaH = (vsv_parms["lambdaC"] * 0.1), # Lysing rate of normal (Healthy) cells
  betaH = (vsv_parms["betaC"] * 0.0025) # Infection rate of normal cells. Normal cell to tumor cell infection rate ratio is 0.0025
))

# Replace weird "bH.bC" type names with just "bH"
names(vsv_parms) <- gsub("\\..*", "", names(vsv_parms))

########## STATE VARIABLES ##########

# Defining the states
state <- c(
  Ci = 0, # Amount of tumor cells infected
  Hi = 0, # Amount of normal cells infected
  v0 = 10^9, # Initial free virions.
  Cs = 4.2*10^9, # Amount of susceptible tumor cells,
  Hs = 10^11 # Amount of susceptible normal cells
  
)

########## TIME PARAMETERS ##########

# Defining the 7 days in amount of hours.
times <- seq(0, 168)
