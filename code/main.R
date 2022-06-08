# Remodeling oncolytic virotherapy, R code.
# Loading prerequisites
library(deSolve)

# Parameterss & state will need to be defined first.
parms_adeno <- c(
  # Parameters governing cell growth
  rH = 0.00275, #The per-capita growth rate of normal cells. Unit: hr^-1
  rC = 0.003, #	The per-capita growth rate of tumor cells. Unit: hr^−1
  KH = 10^11, #		The carrying capacity of normal cells. Unit: cells
  KC = 1.47*10^12, #		The carrying capacity of tumor cells. Unit: cells
  
  # Parameters for the Adenovirus.
  lambdaC = 1/48, # Lysing rate of tumor cells. Unit: cell^-1hr^-1
  bC = 1000, # Burst size from lysing infected tumor cells. Unit: -
  omega = 1, # Viral clearance rate. Unit: virus^-1hr^-1
  betaC = range(5.25*10^-12, 5.25*10^-13.5) # Uptake/encounter/infection rate of tumor cells. Unit: viruses cell^−1hr^−1
)

# Parameterss & state will need to be defined first.
parms_hsv <- c(
  # Parameters governing cell growth
  rH = 0.00275, #The per-capita growth rate of normal cells. Unit: hr^-1
  rC = 0.003, #	The per-capita growth rate of tumor cells. Unit: hr^−1
  KH = 10^11, #		The carrying capacity of normal cells. Unit: cells
  KC = 1.47*10^12, #		The carrying capacity of tumor cells. Unit: cells
  
  # Parameters for the HSV.
  lambdaC = 1/18, # Lysing rate of tumor cells. Unit: cell^-1hr^-1
  bC = 50, # Burst size from lysing infected tumor cells. Unit: -
  omega = 0.025, # Viral clearance rate. Unit: virus^-1hr^-1
  betaC = range(2.5*10^-12, 2.5*10^-13.5) # Uptake/encounter/infection rate of tumor cells. Unit: viruses cell^−1hr^−1
)

# Parameterss & state will need to be defined first.
parms_vsv <- c(
  # Parameters governing cell growth
  rH = 0.00275, #The per-capita growth rate of normal cells. Unit: hr^-1
  rC = 0.003, #	The per-capita growth rate of tumor cells. Unit: hr^−1
  KH = 10^11, #		The carrying capacity of normal cells. Unit: cells
  KC = 1.47*10^12, #		The carrying capacity of tumor cells. Unit: cells
  
  # Parameters for VSV.
  lambdaC = 1/24, # Lysing rate of tumor cells. Unit: cell^-1hr^-1
  bC = 1350, # Burst size from lysing infected tumor cells. Unit: -
  omega = 0.244, # Viral clearance rate. Unit: virus^-1hr^-1
  betaC = range(5*10^-12.5, 5*10^-14) # Uptake/encounter/infection rate of tumor cells. Unit: viruses cell^−1hr^−1
)


state <- c(
  test1 = 4.74,
  test2 = 267,
)
