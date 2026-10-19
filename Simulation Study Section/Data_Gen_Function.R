##### In this function, we generate data from true structure.First of all, we generate True_R of size I times P 
##### and True_C of size J times Q and then define True_G. Using True_R, True_C 
##### and True_G we generate True_DC of size I times J and add some Gaussian noise to True_DC in order to generate DC.

Data_Gen_Function <- function(I, J, P, Q){     
  ##### SNR stands for how many times noise variance is larger than true variance
if (P < I) {
  True_R <- as.matrix(Unequal_Randompartition_Function (I, P, 0.5))
} else {
  True_R <- diag(I)        ##### Otherwise take identity matrix for R
}
##### Now we generate initial random column partition matrix C  
if (Q < J){
  True_C <- as.matrix(Randompartition_function (J, Q))
} else {
  True_C <- diag(J)      ##### Otherwise take identity matrix for C
}
  CardC <- colSums(True_C) # column cluster cardinalities
  int_con <- I*J # a constant simply to avoid small gamma values
  
  # define gamma's for (P,Q) = (2,2), (2,3) and (3,2) (3,3), respectively
  if (P==2 && Q==2){
  ### We want to restrict entries of gamma matrix bewteen -0.5 to +0.5, so we multiplied a constant factor to this matrix.
    True_G =  matrix(c(+1*int_con/(0.5*I*CardC[1]),-1*int_con/(0.5*I*CardC[2]),
                                      -1*int_con/(0.5*I*CardC[1]),+1*int_con/(0.5*I*CardC[2])),
                                                              nrow = 2, ncol = 2,byrow = TRUE)}
  if (P==2 && Q==3){
    ### We want to restrict entries of gamma matrix bewteen -0.5 to +0.5, so we multiplied a constant factor to this matrix.
    True_G = matrix(c(+1*int_con/(0.5*I*CardC[1]), 0, -1*int_con/(0.5*I*CardC[3]),
                                      -1*int_con/(0.5*I*CardC[1]), 0, +1*int_con/(0.5*I*CardC[3])),
                                    nrow = 2, ncol = 3,byrow = TRUE)}
  if (P==3 && Q==2){
    True_G = matrix(c(-1*int_con/(0.5*I*CardC[1]),+1*int_con/(0.5*I*CardC[3]),
                      0,0,
                      +1*int_con/(0.25*I*CardC[1]),-1*int_con/(0.25*I*CardC[3])),
                    nrow = 3, ncol = 2,byrow = TRUE)}
  
  if (P==3 && Q==3){
    True_G = matrix(c(-1*int_con/(0.5*I*CardC[1]),0,+1*int_con/(0.5*I*CardC[3]),
                                       0,0,0,
                                      +1*int_con/(0.25*I*CardC[1]),0,-1*int_con/(0.25*I*CardC[3])),
                                                                 nrow = 3, ncol = 3,byrow = TRUE)}

  True_T <- True_R%*%True_G%*%t(True_C)  ### This is Gamma in model equation (i.e. equation (1))
  mu         <- 0 ### We could change this value and its affects proportion of zeros and ones in the binary data 
  alpha_i    <- matrix(rnorm(I*J, mean=0, sd=0.5), nrow = I, ncol =J)
  beta_j     <- matrix(seq(from = -0.5, to = 0.5, length.out = J), I, J, byrow = T)    
  epsilon_ij <- matrix(rnorm(I*J, mean=0, sd=1), nrow = I, ncol =J)
  data       <- mu + alpha_i + beta_j +  True_T + epsilon_ij
  ### Note that if we make data doubly centered then row, column and overall mean are zero so we don't use doubly centered data.
  ### Now we make dichotomous data using the procedure described by Dr Jan
  median_data <- median(data)
  binary_data <- 1*(data > median_data) ### data > 0 gives a matrix of logical and multiplying with one gives a binary data matrix
  ### For making a vector of group label's using the binary data
  overall_mean  <- mean(binary_data)
  ###For IRT based tests, we need binary grouping of rows.
  row_means     <- rowMeans(binary_data)
  group_labels  <- 1*(row_means <= overall_mean)
  ### Now we count number of rows or columns have constants 1's or 0's.
  num_of_rows_constant_1 <- sum(apply(binary_data == "1", 1, all))  ### Here 1 is for row dimension of the matrix
  num_of_rows_constant_0 <- sum(apply(binary_data == "0", 1, all))
  num_of_cols_constant_1 <- sum(apply(binary_data == "1", 2, all))  ### Here 2 is for col dimension of the matrix
  num_of_cols_constant_0 <- sum(apply(binary_data == "0", 2, all))
  ### Now we ensure to remove rows or columns have constants 1's or 0's.
  # binary_data <- binary_data[!apply(binary_data == "1", 1, all), ] ### remove rows with constant 1's 
  # binary_data <- binary_data[!apply(binary_data == "0", 1, all), ] ### remove rows with constant 0's 
  # binary_data <- binary_data[, !apply(binary_data == "1", 2, all)] ### remove cols with constant 1's 
  # binary_data <- binary_data[, !apply(binary_data == "0", 2, all)] ### remove cols with constant 0's 
  my_list <- list("binary_data"=binary_data, "group_labels"=group_labels,
                  "num_of_rows_constant_1"=num_of_rows_constant_1, "num_of_rows_constant_0"=
                    num_of_rows_constant_0, "num_of_cols_constant_1"=num_of_cols_constant_1,
                  "num_of_cols_constant_0"=num_of_cols_constant_0)
  return (my_list)
} ##### ending of function
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  










