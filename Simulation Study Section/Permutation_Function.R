Permutation_Function <- function(D, P, Q, Nruns, permutations, alpha_level){
  I<- nrow(D)
  J<- ncol(D)
  Obs_Log_LR_REMAXINT      <- vector("numeric", 1L)           ### To store results for observed sample based on REMAXINT.
  Obs_Log_LR_EReMI         <- vector("numeric", 1L)           ### To store results for observed sample based on EReMI.
  Perm_Log_LR_REMAXINT     <-vector("numeric", permutations)  ### To store results for permuted distribution for REMAXINT.
  Perm_Log_LR_EReMI        <-vector("numeric", permutations)  ### To store results for permuted distribution for EReMI.
  Crit_Value_Perm_REMAXINT <-vector("numeric", 1L)            ### REMAXINT
  Crit_Value_Perm_EReMI    <-vector("numeric", 1L)            ### EReMI

  ### Fitting observed data with the REMAXINT.
  Output_REMAXINT     <- REMAXINT (D, P, Q, Nruns)
  Z_REMAXINT          <- Output_REMAXINT$BestR
  K_REMAXINT          <- Output_REMAXINT$BestC
  Gamma_REMAXINT      <- Output_REMAXINT$BestG
  Constant_omegahat   <- ones(1,P)*(1/P)
  T_REMAXINT          <- Z_REMAXINT%*%Gamma_REMAXINT%*%t(K_REMAXINT)
  ### Doubly center is required to compute Obs_Log_LR_REMAXINT
  Xi <- as.matrix(rowMeans(D))%*% ones(1,J)   ### matrix of row means 
  Xj <- ones(I,1)%*%t(as.matrix(colMeans(D))) ### matrix of column means
  DC <- D - Xi - Xj + ones(I,J)*mean(D)
  Obs_Log_LR_REMAXINT  <- Log_LR_Test_Statistic (DC,I,J,T_REMAXINT,Constant_omegahat)
  
  ### Fitting observed data with the EReMI
  Output_EReMI        <- E_ReMI (D, P, Q, Nruns, Z_REMAXINT, K_REMAXINT, Gamma_REMAXINT)
  Z_EReMI             <- Output_EReMI$BestR
  K_EReMI             <- Output_EReMI$BestC
  Gamma_EReMI         <- Output_EReMI$BestG
  Omegahat_EReMI      <- Output_EReMI$Best_O
  T_EReMI             <- Z_EReMI%*%Gamma_EReMI%*%t(K_EReMI)
  Obs_Log_LR_EReMI    <- Log_LR_Test_Statistic (DC,I,J,T_EReMI,Omegahat_EReMI)
  
  for (perm in 1:permutations){
    
    ###perm1 <- DC[sample(nrow(DC),replace=F),] ### Changing positions of rows, since rows are unordered catagorical variables. Therefore, the position of a row is arbitrary and so this permutation does not change the data set. 
    ###perm2 <- DC[,sample(ncol(DC),replace=F)] ### Changing positions of columns, since columns are unordered catagorical variables. Therefore, the position of a column is arbitrary and so this permutation does not change the data set.
    ###perm3 <- DC[sample(nrow(DC),replace=F),sample(ncol(DC),replace=F)] ### Changing positions of rows and columns, since rows/columns are unordered catagorical variables. Therefore, the position of a row/col is arbitrary and so this permutation does not change the data set.  
    ###perm4 <- randomizeMatrix(DC,null.model = "frequency",iterations = 1000) ### The option null.model="frequency" maintains column sums and null.model="richness" maintains row sums (Note: Given matrix may be square or non-squared). 
    ###perm5 <- rmperm(DC) ### This command is used when we do not maintain rows and columns sums but matrix should be square which is not always the case.
    
    Perm_D  <- matrix(sample(D, replace = FALSE), nrow=nrow(D), ncol=ncol(D))
   
     ### Fitting each permuted data with the REMAXINT.
    Output_Perm_REMAXINT   <- REMAXINT (Perm_D, P, Q, Nruns)
    Z_Perm_REMAXINT        <- Output_Perm_REMAXINT$BestR
    K_Perm_REMAXINT        <- Output_Perm_REMAXINT$BestC
    Gamma_Perm_REMAXINT    <- Output_Perm_REMAXINT$BestG
    Constant_omegahat      <- ones(1,P)*(1/P)
    T_Perm_REMAXINT        <- Z_Perm_REMAXINT%*%Gamma_Perm_REMAXINT %*%t(K_Perm_REMAXINT)
    ### Doubly center Perm_D is required to compute Perm_Log_LR_REMAXINT
    Xi <- as.matrix(rowMeans(Perm_D))%*% ones(1,J)   ### matrix of row means 
    Xj <- ones(I,1)%*%t(as.matrix(colMeans(Perm_D))) ### matrix of column means
    Perm_DC <- Perm_D - Xi - Xj + ones(I,J)*mean(Perm_D)
    Perm_Log_LR_REMAXINT[perm] <- Log_LR_Test_Statistic(Perm_DC,I,J,T_Perm_REMAXINT,Constant_omegahat)
    
    ### Fitting each permuted sample with the EReMI.
    Output_Perm_EReMI   <- E_ReMI (Perm_D, P, Q, Nruns, Z_Perm_REMAXINT, K_Perm_REMAXINT, Gamma_Perm_REMAXINT)
    Z_Perm_EReMI        <- Output_Perm_EReMI$BestR
    K_Perm_EReMI        <- Output_Perm_EReMI$BestC
    Gamma_Perm_EReMI    <- Output_Perm_EReMI$BestG
    Omegahat_Perm_EReMI <- Output_Perm_EReMI$Best_O
    T_Perm_EReMI        <- Z_Perm_EReMI%*%Gamma_Perm_EReMI%*%t(K_Perm_EReMI)
    Perm_Log_LR_EReMI[perm] <- Log_LR_Test_Statistic(Perm_DC,I,J,T_Perm_EReMI,Omegahat_Perm_EReMI)

  } ### end for loop for rep
    P_value_Perm_REMAXINT     <- sum(Perm_Log_LR_REMAXINT > Obs_Log_LR_REMAXINT)/permutations
    P_value_Perm_EReMI        <- sum(Perm_Log_LR_EReMI    > Obs_Log_LR_EReMI)/permutations
    Perm_distn_REMAXINT       <- sort(Perm_Log_LR_REMAXINT)
    Crit_Value_Perm_REMAXINT  <- Perm_distn_REMAXINT[ceiling((1-alpha_level)*permutations)]
    Perm_distn_EReMI          <- sort(Perm_Log_LR_EReMI)
    Crit_Value_Perm_EReMI     <- Perm_distn_EReMI[ceiling((1-alpha_level)*permutations)]

    
    my_list <- list("Obs_Log_LR_REMAXINT" = Obs_Log_LR_REMAXINT, "Obs_Log_LR_EReMI" = Obs_Log_LR_EReMI, 
                    "Perm_Log_LR_REMAXINT"=Perm_Log_LR_REMAXINT, "Perm_Log_LR_EReMI "=Perm_Log_LR_EReMI , 
                    "Crit_Value_Perm_REMAXINT" = Crit_Value_Perm_REMAXINT, "Crit_Value_Perm_EReMI" = Crit_Value_Perm_EReMI, 
                    "P_value_Perm_REMAXINT"= P_value_Perm_REMAXINT, "P_value_Perm_EReMI"= P_value_Perm_EReMI,
                    "Z_EReMI"=Z_EReMI, "K_EReMI"=K_EReMI, "Z_REMAXINT"=Z_REMAXINT, "K_REMAXINT"=K_REMAXINT,
                    "Gamma_REMAXINT"=Gamma_REMAXINT, "Gamma_EReMI"=Gamma_EReMI)
    return(my_list)
}   ### end function












