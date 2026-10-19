Randompartition_function <- function(i,p) {
  A <- zeros (i,p)
  I <- diag(p)
  A[1:p, ] <- I                                         ##### for first p rows we generate identity matrix of order p
  A[c(p+1):i, ] <- I[sample(p,i-p,replace=TRUE),]       ##### for the remaining i-p rows we generate discrete uniform numbers 
  A <- A[sample(i),]                                    ##### ranging between p and i-p using sample command and assign rows of
  return(A)}                                            ##### I according to these randomly generated uniform numbers 
                                                        ##### Random permutation of rows for A on line# 06



