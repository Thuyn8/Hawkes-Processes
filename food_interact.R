library(stelfi)
library(dplyr)
library(gridExtra)
library(ggplot2)
# Loading  the two feeding datasets
S2 <- read.csv("~/Documents/Hawkes Processes/report/feeding/S2Data.csv") %>% 
  select(-c("Date", "Time", "MSec", "Sample")) 
 

S20 <- read.csv("~/Documents/Hawkes Processes/report/feeding/S20Data.csv") %>% 
  select(-c("Date", "Time", "MSec", "Sample")) 

# Two flies didn't interact with S2 food : D6W5, D19W10. 
## All flies interact at least once with S20
sapply(S2, max)

# Creating data with threshold >=30
thres <- 30
S2_thres <- as.data.frame(sapply(S2[, !(colnames(S2) %in% c("D6W5", "D19W10"))], 
                                 function(i){as.integer(i>=thres)}))%>%
  mutate(Time = 1: nrow(S2))
S20_thres <- as.data.frame(sapply(S20, function(i){as.integer(i>=thres)})) %>%
  mutate(Time = 1: nrow(S20))

# creating L2, L20 (time points for each fly)
L2 <- lapply(S2_thres[,-ncol(S2_thres)], function(i){S2_thres$Time[i == 1]})
L20 <-lapply(S20_thres[,-ncol(S20_thres)], function(i){S20_thres$Time[i ==1]})

params <- c(mu = 0.05, alpha = 0.05, beta = 1)


## estimate mu_0, alpha, beta
est_S2 <-as.data.frame(t(do.call(cbind, lapply(L2, function(i){get_coefs(fit_hawkes(times = i, parameters = params))[1:3,1]}))))

est_S20 <-as.data.frame(t(do.call(cbind, lapply(L20, function(i){get_coefs(fit_hawkes(times = i, parameters = params))[1:3,1]}))))

## obtain intensity rates

rate2 <- lapply(L2, function(i){show_hawkes_GOF(fit_hawkes(times = i,parameters = params), 
                                                plot = F,
                                                return_values = T)$interarrivals})
rate2 <-lapply(L20, function(i){show_hawkes_GOF(fit_hawkes(times = i, parameters = params),
                                                plot = F,
                                                return_values = T)$interarrivals})

## plots

#show_hawkes_GOF(fit_hawkes(times = L2[[i]], parameters = params))


      
