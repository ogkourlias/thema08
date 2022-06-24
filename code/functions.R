source("parameters.R")

########## MODEL FUNCTION ##########

vir <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    # Function governing cancerous infected cells.
    dCi <- betaC * Cs * v0 - lambdaC * Ci
    
    # Function governing healthy infected cells.
    dHi <- betaH * Hs * v0 - lambdaH * Hi
    
    # Function governing amount of virus particles.
    dv <- betaC * lambdaC * Ci + bH * bC * lambdaH * Hi - betaH * Hs * v0 - betaC * Cs * v0 - omega * v0
    
    # Function governing cancerous susceptible cells.
    dCs <- rC * Cs * (1 - ( (Cs + Ci) / KC) ) - Cs * betaC * v0
    
    # Function governing healthy susceptible cells.
    dHs <- rH * Hs * (1 - ( (Hs+Hi)/KH ) ) - Hs * betaH * v0
    
    return(list(c(dCi, dHi, dv, dCs, dHs)))
  }
  )
}

########## PLOTTING FUNCTIONS ##########

double_y_plot <- function(res_df, x, y.prim, y.sec, ylim.prim, ylim.sec) {
  #' @description Function which plots two series using ggplot2, where the second series has it's own scaled y-axis.
  #' @param res_df A data frame containing data from a DeSolve ode function call.
  #' @param x Data to be displayed on the x-axis.
  #' @param y.prim Data to be displayed on the primary y-axis.
  #' @param y.sec Data to be displayed on the secondary y-axis.
  #' @param ylim.prim ylim values for the primary y-axis.
  #' @param ylim.sec ylim values for the secondary y-axis.

  y.sec.data <- y.sec
  
  fit = lm(b ~ . + 0,
           tibble::tribble(
             ~a, ~s, ~b,
             1, (ylim.sec[1] - mean(y.sec.data)) / sd(y.sec.data), ylim.prim[1],
             1, (ylim.sec[2] - mean(y.sec.data)) / sd(y.sec.data), ylim.prim[2]))
  
  a <- fit$coefficients["a"]
  s <- fit$coefficients["s"]
  
  dyp <- ggplot(res_df, aes(x = x, y = y.prim)) +
    geom_line(aes(colour = "CS")) +
    geom_line(aes(y = (a + ((y.sec - mean(y.sec.data)) / sd(y.sec.data)) * s ), colour = "V0")) +
    scale_y_continuous(
      name = "Cs", 
      limits = ylim.prim,
      sec.axis = sec_axis(name = "v0", 
                          trans = ~ (. - a) / s * sd(y.sec.data) + mean(y.sec.data))
    ) +
    labs(x = "Time")
  
  return(dyp)
  
}
