
# Remodeling oncolytic virotherapy, R code.
# Loading prerequisites
library(deSolve)
library(ggplot2)
library(gridExtra)

# Parameters & state will need to be defined first.
# Defining parameters that stay the same
base_parameters <- c(
  # Parameters governing cell growth
  rH = 0.00275, #The per-capita growth rate of normal cells. Unit: hr^-1
  rC = 0.003, #	The per-capita growth rate of tumor cells. Unit: hr^−1
  KH = 10^11, #		The carrying capacity of normal cells. Unit: cells
  KC = 1.47*10^12 #		The carrying capacity of tumor cells. Unit: cells
)

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
    # Function governing cancerous infected cells.
    dCi <- betaC * Cs * v0 - lambdaC * Ci
    
    # Function governing healthy infected cells.
    dHi <- betaH * Hs * v0 - lambdaH * Hi
    
    # Function governing amount of virus particles.
    dv <- betaC * lambdaC * Ci + bH * bC *lambdaH * Hi - betaH * Hs * v0 - betaC * Cs * v0 - omega * v0
    
    # Function governing cancerous susceptible cells.
    dCs <- rC * Cs * (1 - ( (Cs + Ci) / KC) ) - Cs * betaC * v0
    
    # Function governing healthy susceptible cells.
    dHs <- rH * Hs * (1 - ( (Hs+Hi)/KH ) ) - Hs * betaH * v0
    
    return(list(c(dCi, dHi, dv, dCs, dHs)))
  }
  )
}

# Defining the desolve outputs to their respective variables.
adeno_res<- ode(y = state, times = times, func = vir, parms = adeno_parms)
hsv_res <- ode(y = state, times = times, func = vir, parms = hsv_parms)
vsv_res <- ode(y = state, times = times, func = vir, parms = vsv_parms)

# Applying post processing Log10 transformations for all values but time.
adeno_res[,2:6] <- log10(adeno_res[,2:6] + 1)
hsv_res[,2:6] <- log10(hsv_res[,2:6] + 1)
vsv_res[,2:6] <- log10(vsv_res[,2:6] + 1)

adeno_res_df <- data.frame(adeno_res)
hsv_res_df <- data.frame(hsv_res)
vsv_res_df <- data.frame(vsv_res)

ylim.prim <- c(0, 10)
ylim.sec <- c(8, 16)

# Adeno virus plotting

V0.temp <- adeno_res_df$v0

fit = lm(b ~ . + 0,
         tibble::tribble(
           ~a, ~s, ~b,
           1, (ylim.sec[1] - mean(V0.temp)) / sd(V0.temp), ylim.prim[1],
           1, (ylim.sec[2] - mean(V0.temp)) / sd(V0.temp), ylim.prim[2]))

a <- fit$coefficients["a"]
s <- fit$coefficients["s"]

reg_adeno <- ggplot(adeno_res_df, aes(x = time, y = Cs)) +
  geom_line(aes(colour = "CS")) +
  geom_line(aes(y = (a + ((v0 - mean(V0.temp)) / sd(V0.temp)) * s ), colour = "V0")) +
  scale_y_continuous(
    name = "Cs", 
    limits = ylim.prim,
    #breaks = c(0, 2, 4, 6, 8, 10),
    sec.axis = sec_axis(name = "v0", 
                        trans = ~ (. - a) / s * sd(V0.temp) + mean(V0.temp))
    #breaks = c(9.0, 9.5, 10.0, 10.5, 11.0))
  ) +
  labs(x = "Time")

# HSV Virus Plotting
V0.temp <- hsv_res_df$v0

fit = lm(b ~ . + 0,
         tibble::tribble(
           ~a, ~s, ~b,
           1, (ylim.sec[1] - mean(V0.temp)) / sd(V0.temp), ylim.prim[1],
           1, (ylim.sec[2] - mean(V0.temp)) / sd(V0.temp), ylim.prim[2]))

a <- fit$coefficients["a"]
s <- fit$coefficients["s"]

reg_hsv <- ggplot(hsv_res_df, aes(x = time, y = Cs)) +
  geom_line(aes(colour = "CS")) +
  geom_line(aes(y = (a + ((v0 - mean(V0.temp)) / sd(V0.temp)) * s ), colour = "V0")) +
  scale_y_continuous(
    name = "Cs", 
    limits = ylim.prim,
    #breaks = c(0, 2, 4, 6, 8, 10),
    sec.axis = sec_axis(name = "v0", 
                        trans = ~ (. - a) / s * sd(V0.temp) + mean(V0.temp))
    #breaks = c(9.0, 9.5, 10.0, 10.5, 11.0))
  ) +
  labs(x = "Time")

# VSV Virus Plotting
V0.temp <- vsv_res_df$v0

fit = lm(b ~ . + 0,
         tibble::tribble(
           ~a, ~s, ~b,
           1, (ylim.sec[1] - mean(V0.temp)) / sd(V0.temp), ylim.prim[1],
           1, (ylim.sec[2] - mean(V0.temp)) / sd(V0.temp), ylim.prim[2]))

a <- fit$coefficients["a"]
s <- fit$coefficients["s"]

reg_vsv <- ggplot(vsv_res_df, aes(x = time, y = Cs)) +
  geom_line(aes(colour = "CS")) +
  geom_line(aes(y = (a + ((v0 - mean(V0.temp)) / sd(V0.temp)) * s ), colour = "V0")) +
  scale_y_continuous(
    name = "Cs", 
    limits = ylim.prim,
    #breaks = c(0, 2, 4, 6, 8, 10),
    sec.axis = sec_axis(name = "v0", 
                        trans = ~ (. - a) / s * sd(V0.temp) + mean(V0.temp))
    #breaks = c(9.0, 9.5, 10.0, 10.5, 11.0))
  ) +
  labs(x = "Time")

reg_full_adeno <- ggplot(data = adeno_res_df, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = Hi, color = "Hi")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  geom_line(aes(y = Hs, color = "Hs")) +
  scale_color_manual(values =  c("red", "blue", "brown", "green", "purple"))

reg_full_hsv <- ggplot(data = hsv_res_df, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = Hi, color = "Hi")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  geom_line(aes(y = Hs, color = "Hs")) +
  scale_color_manual(values =  c("red", "blue", "brown", "green", "purple"))

reg_full_vsv <- ggplot(data = vsv_res_df, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = Hi, color = "Hi")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  geom_line(aes(y = Hs, color = "Hs")) +
  scale_color_manual(values =  c("red", "blue", "brown", "green", "purple"))

# Increasing the burst and lysing rates of normal cells
adeno_parms["bH"] = (adeno_parms["bC"] * 0.3) # Burst size of normal (Healthy) cells. Unit : -
adeno_parms["lambdaH"] = (adeno_parms["lambdaC"] * 0.3) # Lysing rate of normal (Healthy) cells
hsv_parms["bH"] = (hsv_parms["bC"] * 0.3) # Burst size of normal (Healthy) cells. Unit : -
hsv_parms["lambdaH"] = (hsv_parms["lambdaC"] * 0.3) # Lysing rate of normal (Healthy) cells
vsv_parms["bH"] = (vsv_parms["bC"] * 0.3) # Burst size of normal (Healthy) cells. Unit : -
vsv_parms["lambdaH"] = (vsv_parms["lambdaC"] * 0.3) # Lysing rate of normal (Healthy) cells

# Defining the desolve outputs to their respective variables.
adeno_res_2<- ode(y = state, times = times, func = vir, parms = adeno_parms)
hsv_res_2 <- ode(y = state, times = times, func = vir, parms = hsv_parms)
vsv_res_2 <- ode(y = state, times = times, func = vir, parms = vsv_parms)

# Applying post processing Log10 transformations for all values but time.
adeno_res_2[,2:6] <- log10(adeno_res_2[,2:6] + 1)
hsv_res_2[,2:6] <- log10(hsv_res_2[,2:6] + 1)
vsv_res_2[,2:6] <- log10(vsv_res_2[,2:6] + 1)

# Making dataframes for plotting
adeno_res_df_2 <- data.frame(adeno_res_2)
hsv_res_df_2 <- data.frame(hsv_res_2)
vsv_res_df_2 <- data.frame(vsv_res_2)

# Plotting
des_adeno <- ggplot(data = adeno_res_df_2, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = Hi, color = "Hi")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  geom_line(aes(y = Hs, color = "Hs")) +
  scale_color_manual(values =  c("red", "blue", "brown", "green", "purple")) +
  guides(color = FALSE, size = FALSE)

des_hsv <- ggplot(data = vsv_res_df_2, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = Hi, color = "Hi")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  geom_line(aes(y = Hs, color = "Hs")) +
  scale_color_manual(values =  c("red", "blue", "brown", "green", "purple")) +
  guides(color = FALSE, size = FALSE)

des_vsv <- ggplot(data = hsv_res_df_2, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = Hi, color = "Hi")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  geom_line(aes(y = Hs, color = "Hs")) +
  scale_color_manual(values =  c("red", "blue", "brown", "green", "purple")) +
  guides(color = FALSE, size = FALSE)


# Increasing the burst and lysing rates of normal cells
adeno_parms["bH"] = (adeno_parms["bC"] * 0.1) # Burst size of normal (Healthy) cells. Unit : -
adeno_parms["lambdaH"] = (adeno_parms["lambdaC"] * 0.1) # Lysing rate of normal (Healthy) cells
hsv_parms["bH"] = (hsv_parms["bC"] * 0.1) # Burst size of normal (Healthy) cells. Unit : -
hsv_parms["lambdaH"] = (hsv_parms["lambdaC"] * 0.1) # Lysing rate of normal (Healthy) cells
vsv_parms["bH"] = (vsv_parms["bC"] * 0.1) # Burst size of normal (Healthy) cells. Unit : -
vsv_parms["lambdaH"] = (vsv_parms["lambdaC"] * 0.1) # Lysing rate of normal (Healthy) cells
adeno_parms["betaH"] = 0 # Infection rate of normal cells. Normal cell to tumor cell infection rate ratio is 0.0025
hsv_parms["betaH"] = 0 # Infection rate of normal cells. Normal cell to tumor cell infection rate ratio is 0.0025
vsv_parms["betaH"] = 0 # Infection rate of normal cells. Normal cell to tumor cell infection rate ratio is 0.0025

# Defining the desolve outputs to their respective variables.
adeno_res_3<- ode(y = state, times = times, func = vir, parms = adeno_parms)
hsv_res_3 <- ode(y = state, times = times, func = vir, parms = hsv_parms)
vsv_res_3 <- ode(y = state, times = times, func = vir, parms = vsv_parms)

# Applying post processing Log10 transformations for all values but time.
adeno_res_3[,2:6] <- log10(adeno_res_3[,2:6] + 1)
hsv_res_3[,2:6] <- log10(hsv_res_3[,2:6] + 1)
vsv_res_3[,2:6] <- log10(vsv_res_3[,2:6] + 1)

# Making dataframes for plotting
adeno_res_df_3 <- data.frame(adeno_res_3)
hsv_res_df_3 <- data.frame(hsv_res_3)
vsv_res_df_3 <- data.frame(vsv_res_3)

# Plotting
spec_adeno <- ggplot(data = adeno_res_df_3, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  scale_color_manual(values =  c("red", "brown", "purple")) +
  guides(color = FALSE, size = FALSE)

spec_hsv <- ggplot(data = hsv_res_df_3, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  scale_color_manual(values =  c("red", "brown", "purple")) +
  guides(color = FALSE, size = FALSE)

spec_vsv <- ggplot(data = vsv_res_df_3, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  scale_color_manual(values =  c("red", "brown", "purple")) +
  guides(color = FALSE, size = FALSE)

