
# Remodeling oncolytic virotherapy, R code.
# Loading prerequisites
library(deSolve)
library(ggplot2)
library(gridExtra)

#setwd("./thema08/code")
source("~/thema08/code/functions.R")

########## DESOLVE W/ DATA PROCESSING ##########

# Defining the desolve outputs to their respective variables.
adeno_res<- ode(y = state, times = times, func = vir, parms = adeno_parms)
hsv_res <- ode(y = state, times = times, func = vir, parms = hsv_parms)
vsv_res <- ode(y = state, times = times, func = vir, parms = vsv_parms)

# Applying post processing Log10 transformations for all values but time.
adeno_res[,2:6] <- log10(adeno_res[,2:6] + 1)
hsv_res[,2:6] <- log10(hsv_res[,2:6] + 1)
vsv_res[,2:6] <- log10(vsv_res[,2:6] + 1)

# Converting results to dataframes
adeno_res_df <- data.frame(adeno_res)
hsv_res_df <- data.frame(hsv_res)
vsv_res_df <- data.frame(vsv_res)

# reg_hsv

adeno_dyp <- double_y_plot(res_df = adeno_res_df, x = adeno_res_df$time, y.prim = adeno_res_df$Cs, 
                     y.sec = adeno_res_df$v0, ylim.prim = c(0, 10), ylim.sec = c(8, 16))
adeno_dyp

hsv_dyp <- double_y_plot(res_df = hsv_res_df, x = hsv_res_df$time, y.prim = hsv_res_df$Cs, 
                           y.sec = hsv_res_df$v0, ylim.prim = c(0, 10), ylim.sec = c(8, 16))
hsv_dyp

vsv_dyp <- double_y_plot(res_df = vsv_res_df, x = vsv_res_df$time, y.prim = vsv_res_df$Cs, 
                         y.sec = vsv_res_df$v0, ylim.prim = c(0, 10), ylim.sec = c(8, 16))
vsv_dyp

reg_full_adeno <- ggplot(data = adeno_res_df, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = Hi, color = "Hi")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  geom_line(aes(y = Hs, color = "Hs")) +
  scale_color_manual(values =  c("red", "blue", "brown", "green", "purple")) +
  ylab("")

reg_full_hsv <- ggplot(data = hsv_res_df, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = Hi, color = "Hi")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  geom_line(aes(y = Hs, color = "Hs")) +
  scale_color_manual(values =  c("red", "blue", "brown", "green", "purple")) +
  ylab("")

reg_full_vsv <- ggplot(data = vsv_res_df, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = Hi, color = "Hi")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  geom_line(aes(y = Hs, color = "Hs")) +
  scale_color_manual(values =  c("red", "blue", "brown", "green", "purple")) +
  ylab("")

# Increasing the burst and lysing rates of normal cells
adeno_parms["bH"] = (adeno_parms["bC"] * 0.5) # Burst size of normal (Healthy) cells. Unit : -
adeno_parms["lambdaH"] = (adeno_parms["lambdaC"] * 0.5) # Lysing rate of normal (Healthy) cells
hsv_parms["bH"] = (hsv_parms["bC"] * 0.5) # Burst size of normal (Healthy) cells. Unit : -
hsv_parms["lambdaH"] = (hsv_parms["lambdaC"] * 0.5) # Lysing rate of normal (Healthy) cells
vsv_parms["bH"] = (vsv_parms["bC"] * 0.5) # Burst size of normal (Healthy) cells. Unit : -
vsv_parms["lambdaH"] = (vsv_parms["lambdaC"] * 0.5) # Lysing rate of normal (Healthy) cells

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
  guides(color = "none", size = "none") +
  ylab("10Log + 1")

des_hsv <- ggplot(data = vsv_res_df_2, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = Hi, color = "Hi")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  geom_line(aes(y = Hs, color = "Hs")) +
  scale_color_manual(values =  c("red", "blue", "brown", "green", "purple")) +
  guides(color = "none", size = "none") +
  ylab("10Log + 1")

des_vsv <- ggplot(data = hsv_res_df_2, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = Hi, color = "Hi")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  geom_line(aes(y = Hs, color = "Hs")) +
  scale_color_manual(values =  c("red", "blue", "brown", "green", "purple")) +
  guides(color = "none", size = "none") +
  ylab("10Log + 1")


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
  scale_color_manual(values =  c("red", "blue", "purple")) +
  guides(color = "none", size = "none") +
  ylab("10Log + 1")

spec_hsv <- ggplot(data = hsv_res_df_3, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  scale_color_manual(values =  c("red", "blue", "purple")) +
  guides(color = "none", size = "none") +
  ylab("10Log + 1")

spec_vsv <- ggplot(data = vsv_res_df_3, aes(x = time)) +
  geom_line(aes(y = Ci, color = "Ci")) +
  geom_line(aes(y = v0, color = "v0")) +
  geom_line(aes(y = Cs, color = "Cs")) +
  scale_color_manual(values =  c("red", "blue", "purple")) +
  guides(color = "none", size = "none") +
  ylab("10Log + 1")

