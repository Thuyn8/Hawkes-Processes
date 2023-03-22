
## 2. permutation
  permutation_test <- function(dat1, dat2, B, f,obs_stat){
  # B: number of permutation
  #combine two lists, called data
  # f : function to obtain the test statistic
  n1 <- nrow(dat1)
  n2 <- nrow(dat2)
  index <- 1: (n1 + n2) 
  data <- c(dat1, dat2)

  # create permuted data
  folds <- lapply(1: B, function(i){
    sample(index)
  })
  
  permuted_statistics <- sapply(folds, function(i){
    dat1 <- data[,i[1:n]]
    dat2 <- data[,i[(n+1): (n+m)]]
    f(avg_CIF(dat1), avg_CIF(dat2))
  })
  
  # p.value 
  mean(permuted_statistics >= obs_stat)
  }
