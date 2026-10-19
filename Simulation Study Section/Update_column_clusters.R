Update_column_clusters <- function(DC, I, J, initR, initC, initG, Q){
  possible_Q <- diag(Q) 
  updated_C <- matrix(0, J, Q)
  new_LL <- matrix(-1*Inf, J, Q)
  max_LL <- matrix(-1*Inf, J, 1)
  #####The following part of the computation of initG only has to be computed once, since it is constant 
  #####so we do that to save time (since it invloves taking an inverse)
  DC_solved_initR <- pinv(t(initR)%*%(initR))%*%t(initR)%*%DC
  initM <- initR%*%initG%*%t(initC) 
  ##### Compute Initial Loglikelihood
  initLL <- Log_Likelihood_function_REMAXINT(DC,I,J,initM)
  for (i in 1:J) { 
    for (j in 1:Q){
      initC[i,] <- possible_Q[j,]             
      new_C <- as.matrix(initC[i,])
      #G <- DC_solved_initR%*%(new_C)%*%pinv(t(new_C)%*%(new_C)) ##### equivalent to line 45 but faster
      #cat(dim(new_C)) ## new_C is a column vector
      #cat("\n")
      M <- initR%*%initG%*%(new_C)
      new_LL[i,j]<- Log_Likelihood_function_REMAXINT(DC[,i],I,J,M)
    } ##### end of for loop across column clusters (NOTICE THAT I AM ENDING IT HERE)
    
    ##### Now we check which assignment was the best and choose that assignment 
    max_LL[i] <- which.max(new_LL[i,])  ##### find which assignment was best
    initC[i,] <- possible_Q[max_LL[i],] ##### change row i in initC to this best assignment
    ##### print best LL on screen (useful for debugging)
    #####cat(new_LL[i,max_LL])  print best LL on screen (useful for debugging)
    #####cat("\n") print any next value on a new line
    updated_C <- initC
  } ##### end of for loop across row elements (NOTICE THIS WAS WHERE YOU ENDED THE FIRSt FOR LOOP )     
  
  ### NOW CHECK WHETHER THERE ARE EMPTY CLUSTERS
  
  empties <- which(colSums(updated_C) == 0)
  
  
  sortfits <- sort.int(max_LL,index.return= TRUE) # sort elements from worst to best in terms of post prob
  worst <- as.matrix(sortfits$ix) # get index for sorted elements
  #worst <- 1:J
  
  ####now put in each of the empty clusters the first elements, as much as there are needed
  if (length(empties)>0){
    for (e in 1:length(empties)){
      
      # below has to be adjusted to deal with 2 or more empty clusters
      #Updated_R[t(worst[1:(1+ceiling(I*threshold_omega))]),] <- ones((1+ceiling(I*threshold_omega)),1)%*%possible_P[empties[e],]
      updated_C[worst[1],] <- possible_Q[empties[e],]
      
      
      # remove the first set of elements from worst before going to another step in for loop
      #remove <- 1:(1+ceiling(J*threshold_omega))
      #cat(length(remove))
      #cat("\n")
      worst <- worst[-1]
      
    }# end empty cluster loop
    
  }# end empty cluster fix
  
  updated_G <-DC_solved_initR%*%(updated_C)%*%pinv(t(updated_C)%*%(updated_C)) # equivalent to line 45 but faster
  my_list <- list("C" = updated_C, "G" = updated_G)
  return(my_list)
} ##### end of fun
