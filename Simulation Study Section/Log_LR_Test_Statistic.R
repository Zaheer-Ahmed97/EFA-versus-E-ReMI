Log_LR_Test_Statistic <- function(DC,I,J,initM,initOmega_hat) {
  E <- DC-initM   ##### The is the difference between DC and reconstructed matrix M
  Log_LR <- 1*sum(initOmega_hat*I*log(initOmega_hat))-(I*J)/2*log(2*pi*sum(E^(2)))+(I*J)/2*log(2*pi*sum(DC^(2))) 
  return(Log_LR)
}