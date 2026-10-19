Update_row_clusters_E_ReMI <- function(DC, I,J, P, Updated_R, Updated_C, Updated_omegahat, Updated_G, Updated_sigma){
  ## function for updating row cluster matrix not assuming equal cluster sizes
  
  # define a matrix containg all possible cluster assignments for a single row
  possible_P <- diag(P)
  
  # compute LL for input solution
  M <- Updated_R%*%Updated_G%*%t(Updated_C)
  maxLL <- Log_Likelihood_function_E_ReMI(DC,I,J,M,Updated_omegahat)
  
  Temp_G <- Updated_G # temporary G matrix, will be changed inside the for loops
  
  # initialize 1 x P vector of LL values
  # that will be used to store LL values of all possible reassignment for a single row element
  TempLL <- matrix(1, 1, P)*(-Inf) # initial value -Inf beacuse LL's can be negative
  
  # Do the necessary computations on the column side of the data once, no need to do it inside the loop
  DC_solved_initC <- DC%*%Updated_C%*%inv(t(Updated_C)%*%(Updated_C))
  # DC_solved_initC is an I x Q matrix where each row contains the row mean of DC per column cluster
  
  indexclus <- 1:P # this vector will be used to define the inner for loop (across the clusters)
  
  #### now start evaluating reassignments for each row i = 1:I
  for (i in 1:I) {
    rowR = Updated_R[i,] # obtain current assignment of row i
    
    # compute current cluster cardinalities
    # (only for checking whether there is need to prevent empty clusters)
    cluster_card <- colSums(Updated_R) 
    
    singletons <- vector()# define singeltons as an empty vector
    ## because an object 'singletons' is needed in one of the if statements below
    
    # if there are empty clusters, then check which clusters these are
    if (any(cluster_card == 1)){
      # find out which clusters are singletons
      singleton_clus <- which(cluster_card == 1)
      # and identify which row elements currently form a singleton cluster
      singletons <- which(Updated_R[,singleton_clus]==1) 
    }
    
    # now only evaluate cluster reassignments if row i is NOT in a singleton cluster
    if (!i %in% singletons){ 
      
      # place current maxLL (i.e., from step i-1) in the position for current cluster assignment
      TempLL[rowR == 1] <- maxLL
      
      ## Next, only evaluate reassignments to other clusters, NOT the current cluster assignment
      ## (this is merely to save computation time)
      for (p in indexclus[rowR==0]){  ## rowR == 0 is TRUE for all clusters except the current one
        
        Temp_R <- Updated_R # define temporary partition matrix R
        Temp_R[i,] <- possible_P[p,] # reassign row i to cluster p
        
        ##option with updating G (computer intensive because of the inverse)
        #Temp_G <- inv(t(Temp_R)%*%(Temp_R))%*%t(Temp_R)%*%DC_solved_initC ### equivalent to line 45 but faster
        ## The following computation of G is faster, avoiding the inverse
        Temp_G <- diag(1/colSums(Temp_R))%*%t(Temp_R)%*%DC_solved_initC
        TempM <- Temp_R%*%Temp_G%*%t(Updated_C)
        
        ## temporarily store LL value of this reassignment in 1 x P vector TempLL
        ## the position in TempLL for the current cluster is already filled in,
        ## with the maxLL value of the previous evaluation (i.e., step i-1)
        TempLL[p] <- Log_Likelihood_function_E_ReMI(DC,I,J,TempM,colSums(Temp_R)/I)
        
      } ##### end for loop across p (the possible cluster assignments)
      
      ##### find index corresponding to best possible assignment (i.e. highest LL)
      max_LL_index  <- which.max(TempLL) 
      maxLL <- TempLL[max_LL_index] # update current maxLL for use in evaluations of cluster assignment of the next row (i.e., row i+1)
      
      ##### change row i in Updated_R to the best assignment
      Updated_R[i,] <- possible_P[max_LL_index,]
      
    }# end if(!i %in% singletons) and go to next row (i.e., step i+1)
    
  } ##### end for loop across i (the row elements)
  
  #  ### NOW CHECK WHETHER THERE ARE EMPTY CLUSTERS
  #empties <- which(colSums(Updated_R) <= 0)
  #
  #
  ####now put in each of the empty clusters the worst fitting elements, as much as there are needed
  #if (length(empties)>0){
  #  for (e in 1:length(empties)){
  #    
  #    # below has to be adjusted to deal with 2 or more empty clusters
  #    #Updated_R[t(worst[1:(1+ceiling(I*threshold_omega))]),] <- ones((1+ceiling(I*threshold_omega)),1)%*%possible_P[empties[e],]
  #    Updated_R[e,] <- possible_P[empties[e],]
  #  
  #  }# end empty cluster loop
  #}# end empty cluster fix
  
  #cat(colSums(Updated_R))
  #cat("\n")
  return(Updated_R)
  
} ##### function end
