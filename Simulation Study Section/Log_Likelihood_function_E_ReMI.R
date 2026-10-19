Log_Likelihood_function_E_ReMI <- function(DC,I,J,initM,initOmega_hat) {
  E <- DC-initM   ##### The is the difference between DC and reconstructed matrix M
  #LL <- J*sum(initOmega_hat*I*log(initOmega_hat))-(I*J)/2*log(2*pi*sum(E^(2)))+(I*J)/2*(log(I*J)-1)
  # dropped the constant
  #LL <- J*sum(initOmega_hat*I*log(initOmega_hat))-(I*J)/2*log(2*pi*sum(E^(2))) 
  LL <- 1*sum(initOmega_hat*I*log(initOmega_hat))-(I*J)/2*log(2*pi*sum(E^(2))) 
  return(LL)
}

