##### Function for intial Unequal random start for rows 
Unequal_Randompartition_Function <- function(i,p,firstp) {
  #firstp is size of the first cluster
  A <- zeros (i,p)
  I <- diag(p)
  Prob_vector <- zeros (1,p)
  Prob_vector[1] <- firstp
  for(k in 2:p) {
    Prob_vector[k] <- (1 - firstp)/(p-1)
  }
  A[1:p, ] <- I                                       ##### for first p rows we generate identity matrix of order p
  A[c(p+1):i, ] <- t(rmultinom(i-p,1,Prob_vector))    ##### for the remaining i-p rows we generate matrix of cluster membership
  A <- A[sample(i),]                           
  #while(sum(colSums(A)==0)>0){
  #  A <- rmultinom(i,1,Prob_vector)
  #}
  return(A)
}

