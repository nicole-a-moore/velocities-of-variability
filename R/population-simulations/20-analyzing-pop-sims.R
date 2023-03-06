## analyzing pop simulations 
## data from script 17
library(tidyverse)
library(ggplot2)

########################################################
###            calculate extinction risk              ## 
########################################################
## read in 100 simulations for each noise colour and calculate time to extinction 
## analyze change in noise colour over time windows 
base_col = seq(0, 3, by = 0.1)
lambda = c(0.1, 0.5, 0.9)

ER_all <- data.frame()

## loop through colours
col = 1
while (col <= length(base_col)) {
  start_colour = base_col[col]
  
  ## loop through icps
  icp = 1
  while (icp <= length(lambda)) {
    
    filename =  paste("data-processed/pop-sims_semi-stable/tvfdpopdynam_icp-", lambda[icp], "_base-col-", 
                      start_colour, "_20-steps.csv", sep = "")
    
    if(file.exists(filename)) {
      
      ## read in data:
      sims <- read.csv(filename)
      
      ## calculate time to extinction for each time series 
      ER <- sims %>%
        group_by(sim) %>%
        summarise(N_stable = first(which(N_stable == 0)),
                  N_inc = first(which(N_inc == 0)),
                  N_dec = first(which(N_dec == 0))) %>%
        mutate(lambda = sims$lambda[1],
               start_colour = start_colour)
      
      ## bind to rest of the data:
      ER_all <- rbind(ER_all, ER)
      
      print(paste("On icp: ", lambda[icp], " and colour: ", base_col[col],
                  sep = ""))
    }
    
    icp = icp + 1
  }
  col = col + 1
}

#write.csv(ER_all, "data-processed/tvfd-pop-sims/17_ER-metrics_20-steps_semi-stable.csv", row.names = F)
ER <- read.csv("data-processed/tvfd-pop-sims/17_ER-metrics_20-steps_semi-stable.csv")

length(which(!is.na(ER$N_stable)))
length(which(!is.na(ER$N_inc)))
length(which(!is.na(ER$N_dec)))

## reformat data for plotting 
ER <- ER %>%
  gather(key = "noise_type", value = "time_to_extinction", c(N_stable, N_inc, N_dec)) 

## replace NAs (cases where population doesn't go extinct) with 200000
ER <- ER %>%
  mutate(time_to_extinction = ifelse(is.na(time_to_extinction), 
                                     200000,
                                     time_to_extinction)) 

## make start noise colour column accurate for decreasing noise
ER <- ER %>%
  mutate(start_colour = ifelse(noise_type == "N_dec",
                               start_colour + 0.1,
                               start_colour))

## calculate mean time to extinction, coefficient of variation in time to extinction 
ER %>%
  group_by(start_colour, noise_type, lambda) %>%
  mutate(mean_pers_time = mean(time_to_extinction)) %>%
  unique(.) %>%
  ggplot(., aes(y = mean_pers_time, x = start_colour, colour = noise_type)) + geom_point() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Starting noise colour", y = "Mean persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) +
  geom_line()+
  scale_color_manual(values = c("red", "green", "black"))

ggsave("figures/all_types.png", width = 10, height = 3)


ER %>%
  group_by(start_colour, noise_type, lambda) %>%
  mutate(mean_pers_time = mean(time_to_extinction)) %>%
  unique(.) %>%
  filter(noise_type == "N_stable") %>%
  ggplot(., aes(y = mean_pers_time, x = start_colour, colour = noise_type)) + geom_point() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Starting noise colour", y = "Mean persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) +
  geom_line() +
  scale_color_manual(values = c("black"))

ggsave("figures/stable.png", width = 10, height = 3)

ER %>%
  group_by(start_colour, noise_type, lambda) %>%
  mutate(cv = sd(time_to_extinction)/mean(time_to_extinction)) %>%
  unique(.) %>%
  ggplot(., aes(y = cv, x = start_colour, colour = noise_type)) + geom_point() +
  theme_light() +
  labs(x = "Starting noise colour", y = "CV persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) +
  geom_line() +
  scale_color_manual(values = c("red", "green", "black"))


########################################################
###            calculate extinction risk              ## 
########################################################
## read in 100 simulations for each noise colour and calculate time to extinction 
## analyze change in noise colour over time windows 
base_col = seq(0, 3, by = 0.1)
lambda = c(0.1, 0.5, 0.9)

ER_all <- data.frame()

## loop through colours
col = 1
while (col <= length(base_col)) {
  start_colour = base_col[col]
  
  ## loop through icps
  icp = 1
  while (icp <= length(lambda)) {
    
    filename =  paste("data-processed/pop-sims_semi-stable/tvfdpopdynam_icp-", lambda[icp], "_base-col-", 
                      start_colour, "_20-steps_1000.csv", sep = "")
    
    if(file.exists(filename)) {
      
      ## read in data:
      ER <- read.csv(filename) %>%
        rename("start_colour" = base_col)
      
      ## bind to rest of the data:
      ER_all <- rbind(ER_all, ER)
      
      print(paste("On icp: ", lambda[icp], " and colour: ", base_col[col],
                  sep = ""))
    }
    
    icp = icp + 1
  }
  col = col + 1
}

ER = ER_all

length(which(!is.na(ER$N_stable)))
length(which(!is.na(ER$N_inc)))
length(which(!is.na(ER$N_dec)))

## reformat data for plotting 
ER <- ER %>%
  gather(key = "noise_type", value = "time_to_extinction", c(N_stable, N_inc, N_dec)) 

## replace NAs (cases where population doesn't go extinct) with 200000
ER <- ER %>%
  mutate(time_to_extinction = ifelse(is.na(time_to_extinction), 
                                     200000,
                                     time_to_extinction)) 

## make start noise colour column accurate for decreasing noise
ER <- ER %>%
  mutate(start_colour = ifelse(noise_type == "N_dec",
                               start_colour + 0.1,
                               start_colour))

## calculate mean time to extinction, coefficient of variation in time to extinction 
ER %>%
  group_by(start_colour, noise_type, lambda) %>%
  mutate(mean_pers_time = mean(time_to_extinction)) %>%
  unique(.) %>%
  ggplot(., aes(y = mean_pers_time, x = start_colour, colour = noise_type)) + geom_point() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Starting noise colour", y = "Mean persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) +
  geom_line()+
  scale_color_manual(values = c("red", "green", "black"))

ggsave("figures/all_types.png", width = 10, height = 3)


ER %>%
  group_by(start_colour, noise_type, lambda) %>%
  mutate(mean_pers_time = mean(time_to_extinction)) %>%
  unique(.) %>%
  filter(noise_type == "N_stable") %>%
  ggplot(., aes(y = mean_pers_time, x = start_colour, colour = noise_type)) + geom_point() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Starting noise colour", y = "Mean persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) +
  geom_line() +
  scale_color_manual(values = c("black"))

ggsave("figures/stable.png", width = 10, height = 3)

ER %>%
  group_by(start_colour, noise_type, lambda) %>%
  mutate(cv = sd(time_to_extinction)/mean(time_to_extinction)) %>%
  unique(.) %>%
  ggplot(., aes(y = cv, x = start_colour, colour = noise_type)) + geom_point() +
  theme_light() +
  labs(x = "Starting noise colour", y = "CV persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) +
  geom_line() +
  scale_color_manual(values = c("red", "green", "black"))

ggsave("figures/all_types_cv.png", width = 10, height = 3)
