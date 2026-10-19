E_ReMI <- function(DC, P, Q, Nruns, initR, initC, initG) { 
  I <- dim(DC)[1]
  J <- dim(DC)[2]
  
  # make sure data is doubly centered
  Xi <- as.matrix(rowMeans(DC))%*% ones(1,J)   ### matrix of row means 
  Xj <- ones(I,1)%*%t(as.matrix(colMeans(DC))) ### matrix of column means
  DC <- DC - Xi - Xj + ones(I,J)*mean(DC)
  
  M <- initR%*%initG%*%t(initC)
  ##### Compute omega_hat, given initial R
  initOmega_hat <- colSums(initR)/I
  init_LL <- Log_Likelihood_function_E_ReMI(DC,I,J,M,initOmega_hat)
  Updated_LL <- init_LL
  ##cat(init_LL) We used this for first simualtion study
  ##cat("\n")
  ##cat("\n")
  
  
  ##### Initialize final estimates, these will be replaced during the algorithm
  BestR <- as.matrix(zeros(I,P))
  BestC <- as.matrix(zeros(J,Q))
  BestG <- as.matrix(zeros(P,Q))
  Best_omega_hat <- as.matrix(zeros(1,P))
  Best_sigma_square <- 0   ##### we can update scalar in this way (choose any value because it will replace during iteration)
  BestLL <- -1*Inf         ##### Maximum badness of fit 
  
  ##### Algorithm Start
  for (i in 1:Nruns){
    ##### Now we generate initial random row partition matrix R
    if (i > 1){ # for i = 1 we user the user input initR, initC and initG, otherwise random
      
      if (P < I) {
        initR <- as.matrix(Unequal_Randompartition_Function (I, P, 0.5))}
      else { initR <- diag(I) }       ##### Otherwise take identity matrix for R
      
      ##### Now we generate initial random column partition matrix C  
      if (Q < J){
        initC <- as.matrix(Randompartition_function (J, Q))}
      else { initC <- diag(J) }      ##### Otherwise take identity matrix for C
      
      ##### Now we sample initial Gamma matrix G
      initG <- as.matrix(solve(t(initR)%*%(initR))%*%t(initR)%*%DC%*%initC%*%solve(t(initC)%*%(initC)))
      #sampleI <- sample(I)
      #sampleJ <- sample(J)
      #initG <- DC[sampleI[1:P],sampleJ[1:Q]]
      
      ##### generate initial omega_hat based on uniform distribution
      #u <- runif(P,0,1)
      #initOmega_hat <- u/sum(u)
      initOmega_hat <- colSums(initR)/I
      
    } # end if i > 1
    
    ##### Now we compute initial reconstructed matrix M, based on initial R, C and G
    initM <- initR%*%initG%*%t(initC)
    
    ###
    #initOmega_hat <- as.matrix(zeros(1,P))+1/P # this option for equal sized clusters
    ##### Now we compute loglikelihood (Below we define a function for that)
    initLL <- Log_Likelihood_function_E_ReMI(DC,I,J,initM,initOmega_hat)
    #cat(initLL)
    #cat("\n")
    #cat("\n")
    ##### Updating Phase
    Updated_R <- initR
    Updated_C <- initC
    Updated_G <- initG
    ################################
    ####CHECK THESE ADJUSTMENTS ####
    ################################
    #Updated_omegahat <- 0 # THIS IS A WRONG CHOICE, AS IT WILL BE USED IN THE FIRST ITERATION OF WHILE LOOP
    Updated_omegahat <- initOmega_hat # THIS IS A CORRECT CHOICE
    #Updated_sigma <- # NO VALUE SPECIFIED HERE??
    Updated_sigma <- sum(sum((DC - initM)^2))/(I*J)
    ################################
    ################################
    ################################
    
    prevLL <- initLL
    iterationLL <- vector() ##### Empty vector (for now, will eventually include all LL iterations)
    diffLL <- 1 
    ##### Now we will go while loop
    while (diffLL > 0.0001){
      ##### First step is to update R and omega_hat keeping C, G fixed (Below we define E_C and M functions for that)
      if( (P > 1) & (P != I) ){
        ##### Function to perfrom step E_C
        
        
        #cat(Updated_G)
        #cat("\n")
        Temp_R <- Update_row_clusters_E_ReMI(DC, I,J, P, Updated_R, Updated_C, Updated_omegahat, Updated_G, Updated_sigma)
        #cat(colSums(Temp_R))
        #cat("\n")
        # check impact on LL
        Temp_M <- Temp_R%*%Updated_G%*%t(Updated_C)
        Temp_omegahat <- colSums(Temp_R)/I
        #Temp_omegahat[Temp_omegahat < 0.15] <- 0.15
        #Temp_omegahat[Temp_omegahat > 0.85] <- 0.85
        
        Temp_LL <- Log_Likelihood_function_E_ReMI(DC,I,J,Temp_M,Temp_omegahat) ###Discuss this step 
        #cat(Temp_LL)
        #cat("\n")
        ## WHAT YOU DID BELOW WAS WRONG BECAUSE Output_E_C_step IS NOT A LIST
        #Updated_R       <- Output_E_C_step$Updated_R
        
        ################################
        ################################
        ################################
        
        ##### Function to perfrom M step
        
        ### YOU ARE NOT USING Updated_C BELOW
        #Output_M_step <- Update_row_cluster_assignments_M_step(DC, I, initC, Updated_R)
        Output_M_step <- Update_G_Omega(DC, I,J, Updated_C, Temp_R)
        
        Temp_G     <- Output_M_step$G
        Temp_M     <- Output_M_step$M
        Temp_sigma <- Output_M_step$Sigma
        Temp_omegahat <- Output_M_step$Omega
        
        #check impact on LL
        Temp_M <- Temp_R%*%Temp_G%*%t(Updated_C)
        Temp_LL <- Log_Likelihood_function_E_ReMI(DC,I,J,Temp_M,Temp_omegahat) ###Discuss this step 
        #cat(Temp_LL)
        #cat("\n")
        
      }  # end if
      else if (P == 1){
        Updated_R <- ones(I, 1)
      }  else    {(P == I)
        Updated_R <- diag(I)}
      
      #####################################################
      ##### Second step is to update C using Updated_R, Updated_G
      if( (Q > 1) & (Q != J) ){
        Updates_col_clus <- Update_column_clusters(DC, I, J, Temp_R, Updated_C, Temp_G, Q)
        ##Updates_col_clus <- Update_column_clusters_adjusted (DC, I, J, Q, Temp_R, Updated_C, Temp_G)
        Temp_C <- Updates_col_clus$C
        Temp_G <- Updates_col_clus$G}
      else if (Q == 1){
        Temp_C <- ones(J, 1)
      }        else {  (Q == J)
        Temp_C <- diag(J)
      }
      
      #cat(colSums(Temp_C))
      #cat("\n")
      #check impact on LL
      Temp_M <- Temp_R%*%Temp_G%*%t(Temp_C)
      Temp_LL <- Log_Likelihood_function_E_ReMI(DC,I,J,Temp_M,Temp_omegahat) ###Discuss this step 
      #cat(Temp_LL)
      #cat("\n")
      ##########################################################
      Output_M_step <- Update_G_Omega(DC, I,J, Temp_C, Temp_R)
      
      Temp_G     <- Output_M_step$G
      Temp_M     <- Output_M_step$M
      Temp_sigma <- Output_M_step$Sigma
      Temp_omegahat <- Output_M_step$Omega
      
      #######################################
      
      
      ##### Third step is to compute current loglikelihood (using that function)
      Temp_M <- Temp_R%*%Temp_G%*%t(Temp_C)
      Temp_LL <- Log_Likelihood_function_E_ReMI(DC,I,J,Temp_M,Temp_omegahat) ###Discuss this step 
      #cat(Temp_LL)
      #cat("\n")
      #cat("\n")
      iterationLL <- c(iterationLL,Temp_LL) ### Save likelihood value at each iteration
      ##### Compute difference between current and previous LL
      diffLL <- Temp_LL - prevLL
      if (diffLL >= 0){
        Updated_R <- Temp_R
        Updated_C <- Temp_C
        Updated_G <- Temp_G
        Updated_omegahat <- Temp_omegahat
        #Best_PP_mat <-Updated_Post_prob # THIS IS THE FIRST OCCURRENCE OF Updated_Post_prob in this function 
        Updated_sigma <- Temp_sigma
        Updated_LL <- Temp_LL
      }
      prevLL <- Temp_LL # set prevLL to current (for use in the next iteration)
    }  ##### ending the while loop
    #cat("\n")
    
    ##### The following is useful for doing multiple runs (e.g., 20).
    if (Temp_LL >= BestLL){     
      BestR <- Updated_R
      BestC <- Updated_C
      BestG <- Updated_G
      Best_omega_hat <- Updated_omegahat
      Best_sigma <- Updated_sigma
      #Best_PP_mat <-Updated_Post_prob # THIS IS THE FIRST OCCURRENCE OF Updated_Post_prob in this function 
      BestLL <- Updated_LL
      ##cat(BestLL) We used this for first simualtion study 
      ##cat("\n")
    }  ##### ending of if statement
    #cat(iterationLL) 
    
  } ##### ending of for loop for Nruns
  ##cat("\n")
  ##cat(BestLL) We used this for first simualtion study
  ##cat("\n")
  my_list <- list("BestR" = BestR, "BestC" = BestC, "BestG"=BestG, "LL" = BestLL, "Best_O"=Best_omega_hat)
  return(my_list)
} ##### ending of function



