##### Function to perfrom M step 
Update_G_Omega <- function(DC, I,J, Updated_C, Updated_R){
  
  DC_solved_initC<- DC%*%(Updated_C)%*%pinv(t(Updated_C)%*%(Updated_C))
  Updated_omegahat <- colsums(Updated_R)/I    ### corresponding cluster sizes
  Updated_G <- pinv(t(Updated_R)%*%(Updated_R))%*%t(Updated_R)%*%DC_solved_initC ### equivalent to line 45 but faster
  Updated_M <- Updated_R%*%Updated_G%*%t(Updated_C)
  Updated_sigma <- sum((DC-Updated_M)^2)/(I*J)
  my_list <- list("Omega" = Updated_omegahat, "G" = Updated_G, "M"=Updated_M, "Sigma"=Updated_sigma) 
  return(my_list)
}









