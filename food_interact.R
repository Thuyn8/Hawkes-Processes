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

## a.  Two flies didn't interact with S2 food : D6W5, D19W10. 
## b. All flies interact at least once with S20
sapply(S2, max)

## c. Creating data with threshold >=30
thres <- 30
S2_thres <- sapply(S2[, !(colnames(S2) %in% c("D6W5", "D19W10"))],  # rm two flies D6W5, D19W10
                                 function(i){as.integer(i>=thres)})%>%
            as.data.frame() %>% 
            mutate(Time = 1: nrow(S2)) # create Time column
S20_thres <- sapply(S20, function(i){as.integer(i>=thres)}) %>%
            as.data.frame() %>%
            mutate(Time = 1: nrow(S20)) # Create Time column

## d. creating L2 and L20 data sets including spike times for each fly
### sugar 2% (s2)
L2 <- lapply(S2_thres[,-ncol(S2_thres)], function(i){S2_thres$Time[i == 1]})
### sugar 20% (s20)
L20 <-lapply(S20_thres[,-ncol(S20_thres)], function(i){S20_thres$Time[i ==1]})

# 2. Independent Hawkes processes
## a. Initial parameters
params <- c(mu = 0.05, alpha = 0.05, beta = 1)


## b. estimate mu_0, alpha, beta using LH 

est_S2 <-t(do.call(cbind, 
                   lapply(L2, function(i){
                     get_coefs(fit_hawkes(times = i, parameters = params))[1:3,1]}))) %>% 
          as.data.frame()

est_S20 <-t(do.call(cbind,
                    lapply(L20, function(i){
                      get_coefs(fit_hawkes(times = i, parameters = params))[1:3,1]}))) %>%
          as.data.frame()

## c. obtain intensity rates for each fly in each data set

rate2 <- lapply(L2, function(i){show_hawkes_GOF(fit_hawkes(times = i,parameters = params), 
                                                plot = F,
                                                return_values = T)$interarrivals})
rate2 <-lapply(L20, function(i){show_hawkes_GOF(fit_hawkes(times = i, parameters = params),
                                                plot = F,
                                                return_values = T)$interarrivals})

## d. plots 
show_hawkes_GOF(fit_hawkes(times = L2[[1]], parameters = params))

# 3. Mutually exciting Hawkes process
## a. 



#  Read the Social interaction data sets
## Control male data
control <- read.csv("~/Documents/Hawkes Processes/report/social_interact/ControlMale.csv") 
dim(control)
## NPF chrimson Male data
npf <- read.csv("~/Documents/Hawkes Processes/report/social_interact/NPFChrimsonMale.csv") 
dim(npf)

## convert data to 0 and 1 using  given threshold 8mm
thres_dist <- 8

control_h <- sapply(control, function(i){as.integer(i>=thres_dist)}) %>% 
            as.data.frame() %>% 
            mutate(Time = 1:nrow(control))
head(control_h)

npf_h <- sapply(npf, function(i){as.integer(i>=thres_dist)}) %>% 
        as.data.frame() %>% 
        mutate(Time = 1:nrow(npf))
head(npf_h)

## create hawkes data
control_hwk <- lapply(control_h[,-ncol(control_h)], function(i){control_h$Time[i ==1]})
## Check missing data
sapply(control_hwk, function(i){100*mean(is.na(i))}) %>%signif(3)  # T2L,T4L have NAs
sapply(npf_h, function(i){100*mean(is.na(i))})%>% signif(3)  # T3R has missing
## missing pattern plots --> monotone missing pattern

aggr(control, numbers = T, gap = 2,
     cex.axis = .5,
     labels = colnames(control),
     sortVars = T,
     ylab = c("Histogram of NAs-Control data", "Pattern of missing data"), only.miss = T)
mtext("Figure 4: Pattern of missing data of Control Male", side =1, line = 7, outer= T )

aggr(npf, numbers = T, gap = 2,
     cex.axis = .5,
     labels = colnames(npf),
     sortVars = T,
     ylab = c("Histogram of NAs-npf data", "Pattern of missing data"), only.miss = T)
mtext("Figure 4: Pattern of missing data of Control Male", side =1, line = 7, outer= T )

## check if any flies have no social interaction
sapply(control, max, na.rm = T) # All flies interacted at least  once 
sapply(npf, max, na.rm = T)     # All flies interacted at least  once 

## create 