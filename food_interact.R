library(stelfi)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(VIM)
library(hawkesbow)

# 1. Loading  the two feeding datasets
S2 <- read.csv("~/Documents/Hawkes Processes/report/feeding/S2Data.csv") %>% 
  select(-c("Date", "Time", "MSec", "Sample")) 


S20 <- read.csv("~/Documents/Hawkes Processes/report/feeding/S20Data.csv") %>% 
  select(-c("Date", "Time", "MSec", "Sample")) 

## a.  Two flies didn't interact with S2 food : D6W5, D19W10. Delete them from S2
sapply(S2, max)
S2 <- S2[, !(colnames(S2) %in% c("D6W5", "D19W10"))]

## b. All flies interacted at least once with S20
sapply(S20, max)

## c. Creating data with threshold >=30
thres <- 30
S2_thres <- sapply(S2, function(i){as.integer(i>=thres)})%>%
            as.data.frame() %>% 
            mutate(Time = 1: nrow(S2))
S2_thres$row_sum <- rowSums(S2_thres[,1:17])

S20_thres <- sapply(S20, function(i){as.integer(i>=thres)}) %>%
            as.data.frame() %>%
            mutate(Time = 1: nrow(S20)) # Create Time column
S20_thres$row_sum <- rowSums(S20_thres[,1:17])

## d. creating L2 and L20 data sets including spike times for each fly
### sugar 2% (s2)
L2 <- lapply(S2_thres[,-ncol(S2_thres)], function(i){S2_thres$Time[i == 1]})
### sugar 20% (s20)
L20 <-lapply(S20_thres[,-ncol(S20_thres)], function(i){S20_thres$Time[i ==1]})

# 2. Independent Hawkes processes
# a. Initial parameters
params <- c(mu = 0.05, alpha = 0.05, beta = 1)

# b. estimate parameters
## 1. baseline intensity: mu
## 2. repr:  reproduction mean
  ## 3. rate: rate of the exponential fertility function
  
    est_S2 <-t(do.call(cbind, 
                       lapply(L2, function(i){
                         mle(i, end = nrow(S2_thres), init = params, 'Exponential')$par})))
    
    
    colnames(est_S2) <-c('mu','repr','rate') 
  
    
    est_S20 <-t(do.call(cbind, 
                        lapply(L20, function(i){
                          mle(i, end = nrow(S20_thres), init = params, 'Exponential')$par})))
    
    
    colnames(est_S20) <-c('mu','repr','rate') 



# c. calculate the Conditional intensity of the two datasets 
# using intensity function from the Hawkesbow package 
# inputs are estimated parameters, and events vectors


    CIF_2 <- vector("list", length(names(S2)))
    names(CIF_2) <- names(S2)
    
    
    CIF_2 <-lapply(1:17, function(i){
      intensity(L2[[i]], 1: max(L2[[i]]),
                fun = est_S2[i,1], 
                repr = est_S2[i,2] , 
                family = 'exp',
                rate = est_S2[i,3])
    })
    
    names(CIF_2) <- names(S2)
    
    CIF_20 <- vector("list", length(names(S20)))
    CIF_20 <- lapply(1:17, function(i){
      intensity(L20[[i]], 1: max(L20[[i]]),
                fun = est_S20[i,1], 
                repr = est_S20[i,2] , 
                family = 'exp',
                rate = est_S20[i,3])
    })
    
    names(CIF_20) <- names(S20)



  ## Aggregate Analysis
    ### estimate the aggregate intensities and plots
    Agg_S2_thres <- S2_thres$Time[S2_thres$row_sum >0]
    Agg_S20_thres <- S2_thres$Time[S20_thres$row_sum >0]



#  Read the Social interaction data sets
  ## Control male data
    control <- read.csv("~/Documents/Hawkes Processes/report/social_interact/ControlMale.csv") 
    
  ## NPF chrimson Male data
    npf <- read.csv("~/Documents/Hawkes Processes/report/social_interact/NPFChrimsonMale.csv") 
  

# Exploring data: checking missing data, non-intereacting flies
  ## check if any flies have no social interaction
    sapply(control, max, na.rm = T)  >=8
    sapply(npf, max, na.rm = T)      >=8    

  ## Check missing data and plot the Missing pattern 
    vis_miss(npf)
    vis_miss(control)
    
  ## convert data to 0 and 1 using  given threshold 8mm
    thres_dist <- 8
    
    control_h <- sapply(control, function(i){ifelse(is.na(i),0, as.integer(i>=thres_dist))}) %>% 
      as.data.frame() %>% 
      mutate(Time = 1:nrow(control))

    npf_h <- sapply(npf, function(i){ifelse(is.na(i),0,as.integer(i>=thres_dist))}) %>% 
      as.data.frame() %>% 
      mutate(Time = 1:nrow(npf))
    
  
  ## create hawkes processes's objects
  control_hwk <- lapply(control_h[,!(colnames(control_h) %in% 'Time')], function(i){control_h$Time[i ==1]})
  
  npf_hwk <- lapply(npf_h[, !(colnames(npf_h) %in% 'Time')], function(i){npf_h$Time[i == 1]})


## # Permutation test

## 1. Average function: finds the avg of CIFs 

  avg_CIF <- function(CIF_list, weights) {
    # Find the maximum time points for all CIFs
      max_time <- max(sapply(CIF_list, function(CIF){length(CIF)}))
      
      CIF_max <- which.max(sapply(CIF_list, function(CIF){length(CIF)})) 
      
      # Interpolate each CIF to the common time points
      for (i in names(CIF_list)[-CIF_max]){
        CIF_list[[i]] <- c(CIF_list[[i]], rep(max(CIF_list[[i]]), max_time- length(CIF_list[[i]])))
      }
      
    # Calculate the weighted average of the interpolated CIFs
    CIF_avg <- rowMeans(do.call(cbind, CIF_list) * weights) 
    return(CIF_avg)
    
  }

#test

  avg_CIF(CIF_2, 1) 
  avg_CIF(CIF_20,1) 

## observed_statistics: MSE of two curves

# observed statisitics:
  obs_statistics <- function(x,y){
    n1 = length(x)
    n2 = length(y)
    if (n1> n2) {y <- c(y, rep(0, n1-n2))}
    else {x <- c(x, rep(0, n2-n1))}
    mse = (x-y)^2
    return(mse)
}

## 2. permutation

  obs_statistics <- obs_statistics(dat1, dat2)
  #combine two lists, called data
  index <- 1: (ncol(dat1) + ncol(dat2)) 
  data <- c(dat1, dat2)
  num_permutations <- 1000
  
  # create permuted data
  permuted_dat <- lapply(1: num_permutations, function(i){
    sample(index)
  })
  
  permuted_stat <- sapply(permuted_dat, function(i){
    dat1 <- data[,i[1:length(i)/2]]
    dat2 <- data[,-i[1:length(i)/2]]
    avg_CIF(dat1) -avg_CIF(dat2)
  })
  
  p-value <- mean(permuted_stat >= obs_statistics)


