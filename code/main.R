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
  betaC = 5.25*10^-12 # Uptake/encounter/infection rate of tumor cells. Unit: viruses cell^−1hr^−1
)

parms_adeno["bH"] = (parms_adeno["bC"] * 0.1) # Burst size of normal (Healthy) cells. Unit : -
parms_adeno["lambdaH"] = (parms_adeno["lambdaC"] * 0.1) # Lysing rate of normal (Healthy) cells
parms_adeno["betaH"] = (parms_adeno["betaC"] * 0.0025) #(4*10^4) # Infection rate of normal cells. Normal cell to tumor cell infection rate ratio is 0.0025

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
  betaC = 2.5*10^-12 # Uptake/encounter/infection rate of tumor cells. Unit: viruses cell^−1hr^−1
)

parms_hsv["bH"] = parms_hsv["bC"] * 0.1 # Burst size of normal (Healthy) cells. Unit : -
parms_hsv["lambdaH"] = parms_hsv["lambdaC"] * 0.1 # Lysing rate of normal (Healthy) cells
parms_hsv["betaH"] = parms_hsv["betaC"] * 0.0025 # Infection rate of normal cells. Normal cell to tumor cell infection rate ratio is 0.0025


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
  betaC = 5*10^-12.5 # Uptake/encounter/infection rate of tumor cells. Unit: viruses cell^−1hr^−1
)

parms_vsv["bH"] = parms_vsv["bC"] * 0.1 # Burst size of normal (Healthy) cells. Unit : -
parms_vsv["lambdaH"] = parms_vsv["lambdaC"] * 0.1 # Lysing rate of normal (Healthy) cells
parms_vsv["betaH"] = parms_vsv["betaC"] * 0.0025 # Infection rate of normal cells. Normal cell to tumor cell infection rate ratio is 0.0025

# Defining the states
state <- c(
  Ci = 0, # Amount of tumor cells infected
  Hi = 0, # Amount of normal cells infected
  v0 = 10^9, # Initial free virions.
  Cs = 4.2*10^9, # Amount of susceptible tumor cells,
  Hs = 10^11 # Amount of susceptible normal cells
  
)

# Defining the 7 days in amount of hours.
times <- seq(0, 168)

# Defining the functions
vir <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    # Function governing healthy susceptible cells.
    dHs <- rH * Hs * (1 - ( (Hs+Hi)/KH ) ) - Hs * betaH * v0
    
    # Function governing cancerous susceptible cells.
    dCs <- rC * Cs * (1 - ( (Cs + Ci) / KC) ) - Cs * betaC * v0
    
    # Function governing cancerous infected cells.
    dCi <- betaC * Cs * v0 - lambdaC * Ci
    
    # Function governing healthy infected cells.
    dHi <- betaH * Hs * v0 - lambdaH * Hi
    
    # Function governing amount of virus particles.
    dv <- betaC * lambdaC * Ci + bH * bC *lambdaH * Hi - betaH * Hs * v0 - betaC * Cs * v0 - omega * v0
    
    return(list(c(dCi, dHi, dv, dCs, dHs)))
  }
  )
}

# Defining the desolve outputs to their respective variables.
adeno_res<- ode(y = state, times = times, func = vir, parms = parms_adeno)
hsv_res <- ode(y = state, times = times, func = vir, parms = parms_hsv)
vsv_res <- ode(y = state, times = times, func = vir, parms = parms_vsv)

# Applying post processing Log10 transformations for all values but time.
adeno_res[,2:6] <- log10(adeno_res[,2:6] + 1)
hsv_res[,2:6] <- log10(hsv_res[,2:6] + 1)
vsv_res[,2:6] <- log10(vsv_res[,2:6] + 1)