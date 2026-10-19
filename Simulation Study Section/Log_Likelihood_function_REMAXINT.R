Log_Likelihood_function_REMAXINT <- function(DC,I,J,initM) {
  E <- DC-initM   ##### The is the difference between DC and reconstructed matrix M
  LL <- -1*sum(E^(2)) 
  return(LL) 
  }