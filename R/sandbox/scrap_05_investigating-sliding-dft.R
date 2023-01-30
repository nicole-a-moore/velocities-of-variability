## trying sliding window dft and seeing effects of window width and time step on spectral exponent 
library(tidyverse)
library(ncdf4)
library(abind)
library(sp)
library(sf)
library(raster)
library(broom)
library(evobiR)
select <- dplyr::select

## read in the data:
slopes_spec_exp <- read.csv("./data-processed/cmcc-cesm_spectral-exponent-across-windows_year-time-step.csv")

## create subset that can be easily looked at:
slopes_sub <- slopes_spec_exp[1:10000,]

## calcualte global average slope and intercept of change in spectral exponent using sliding window across window widths:
slidingwindow_estimates <- slopes_spec_exp %>%
  group_by(lat, long, time_window_width) %>%
  do(tidy(lm(., formula = spec_exp ~ time_window_start))) %>%
  ungroup() %>%
  select(time_window_width, term, estimate, p.value, std.error) %>%
  group_by(time_window_width, term) %>%
  do(mutate(., estimate = mean(.$estimate), p.value = mean(.$p.value), 
            std.error = mean(.$std.error))) %>%
  unique() 
  
## compare with time hop estimates:
time_hop_estimates <- slopes_spec_exp %>%
  mutate(width = as.numeric(str_split_fixed(as.character(slopes_spec_exp$time_window_width), pattern = " ", n = 2)[,1])) %>%
  filter((time_window_start-1)%%(365*width) == 0) %>% ## make data as if it used a time hop beginning at 1850
  group_by(lat, long, time_window_width) %>%
  do(tidy(lm(., formula = spec_exp ~ time_window_start))) %>%
  ungroup() %>%
  select(time_window_width, term, estimate, p.value, std.error) %>%
  group_by(time_window_width, term) %>%
  do(mutate(., estimate = mean(.$estimate), p.value = mean(.$p.value), 
            std.error = mean(.$std.error))) %>%
  unique() 


## split data to show what would happen if you started the time hop at 1850, 1851, 1852, 1853, and 1854
## 1850:
time_hop <- slopes_spec_exp %>%
  mutate(width = as.numeric(str_split_fixed(as.character(slopes_spec_exp$time_window_width), pattern = " ", n = 2)[,1])) %>%
  filter((time_window_start-1)%%(365*width) == 0) %>%
  mutate(start_year = "1850") 

## 1851:
one <- slopes_spec_exp %>%
  mutate(width = as.numeric(str_split_fixed(as.character(slopes_spec_exp$time_window_width), pattern = " ", n = 2)[,1])) %>%
  filter((time_window_start-1)%%(365*width) == 365) %>%
  mutate(start_year = "1851") 

## 1852:
two <- slopes_spec_exp %>%
  mutate(width = as.numeric(str_split_fixed(as.character(slopes_spec_exp$time_window_width), pattern = " ", n = 2)[,1])) %>%
  filter((time_window_start-1)%%(365*width) == 365*2) %>%
  mutate(start_year = "1852") 

## 1853:
three <- slopes_spec_exp %>%
  mutate(width = as.numeric(str_split_fixed(as.character(slopes_spec_exp$time_window_width), pattern = " ", n = 2)[,1])) %>%
  filter((time_window_start-1)%%(365*width) == 365*3) %>%
  mutate(start_year = "1853")

## 1854:
four <- slopes_spec_exp %>%
  mutate(width = as.numeric(str_split_fixed(as.character(slopes_spec_exp$time_window_width), pattern = " ", n = 2)[,1])) %>%
  filter((time_window_start-1)%%(365*width) == 365*4) %>%
  mutate(start_year = "1854") 


## plot start year against average spectral exponent to see how time hop start time affects results 
## expect: if no effect on results, should be a straight horizontal line
time_hops <- rbind(time_hop, one, two, three, four)

time_hops$time_window_width <- factor(time_hops$time_window_width, 
                                      levels = c("5 years", "6 years",
                                                 "7 years", "8 years",
                                                 "9 years", "10 years"))

## plot global average change in spectral exponent as function of start year
global_avg <- time_hops %>%
  group_by(lat, long, time_window_width, start_year) %>%
  do(tidy(lm(., formula = spec_exp ~ time_window_start))) %>%
  filter(term == "time_window_width")
  ungroup() %>%
  select(time_window_width, start_year, estimate) %>%
  group_by(time_window_width, start_year) %>%
  do(mutate(., estimate = mean(.$estimate))) %>%
  unique() 


global_avg %>%
  ggplot(., aes(x = start_year, y = estimate, colour = time_window_width, group = time_window_width)) + geom_point() + geom_path() +
  labs(x = "Time hop start year", colour = "Window width", y = "Global average change in spectral exponent")



sensitivity <- time_hops %>% 
  group_by(start_year, time_window_start, time_window_width) %>%
  do(mutate(., spec_exp = mean(.$spec_exp))) %>% ## calculate global average spectral exponent for each time hop across start years
  ungroup() %>%
  select(start_year, time_window_start, time_window_width, spec_exp) %>%
  unique() %>%
  mutate(ylab = (time_window_start-1)/365+1850)

windowsensitivity <- sensitivity %>%
  filter(start_year == "1850") %>%
  ggplot(., aes(x = ylab, y = spec_exp, colour = time_window_width)) +
  geom_point() + geom_path() + stat_smooth(method = "lm", se = FALSE) +
  labs(x = "Year", y = "Spectral Exponent", colour = "Window width")


sensitivity_5 <- sensitivity %>%
  filter(time_window_width == "5 years") %>%
  ggplot(., aes(x = ylab, y = spec_exp, colour = start_year, group = start_year))+
  geom_point() + geom_path() + stat_smooth(method = "lm", se = FALSE) +
  labs(x = "Year", y = "Spectral Exponent", colour = "Time hop start year", title = "5 year window width")

sensitivity_6 <- sensitivity %>%
  filter(time_window_width == "6 years") %>%
  ggplot(., aes(x = ylab, y = spec_exp, colour = start_year, group = start_year))+
  geom_point() + geom_path() + stat_smooth(method = "lm", se = FALSE) +
  labs(x = "Year", y = "Spectral Exponent", colour = "Time hop start year", title = "6 year window width")

sensitivity_7 <- sensitivity %>%
  filter(time_window_width == "7 years") %>%
  ggplot(., aes(x = ylab, y = spec_exp, colour = start_year, group = start_year))+
  geom_point() + geom_path() + stat_smooth(method = "lm", se = FALSE) +
  labs(x = "Year", y = "Spectral Exponent", colour = "Time hop start year", title = "7 year window width")

sensitivity_8 <- sensitivity %>%
  filter(time_window_width == "8 years") %>%
  ggplot(., aes(x = ylab, y = spec_exp, colour = start_year, group = start_year))+
  geom_point() + geom_path() + stat_smooth(method = "lm", se = FALSE) +
  labs(x = "Year", y = "Spectral Exponent", colour = "Time hop start year", title = "8 year window width")

sensitivity_9 <- sensitivity %>%
  filter(time_window_width == "9 years") %>%
  ggplot(., aes(x = ylab, y = spec_exp, colour = start_year, group = start_year))+
  geom_point() + geom_path() + stat_smooth(method = "lm", se = FALSE)+
  labs(x = "Year", y = "Spectral Exponent", colour = "Time hop start year", title = "9 year window width")

sensitivity_10 <- sensitivity %>%
  filter(time_window_width == "10 years") %>%
  ggplot(., aes(x = ylab, y = spec_exp, colour = start_year, group = start_year))+
  geom_point() + geom_path() + stat_smooth(method = "lm", se = FALSE)+
  labs(x = "Year", y = "Spectral Exponent", colour = "Time hop start year", title = "10 year window width")


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(sensitivity_5)

sensitivity_plot <- grid.arrange(arrangeGrob(sensitivity_5 + theme(legend.position="none"),
                                             sensitivity_6 + theme(legend.position="none"),
                                             sensitivity_7 + theme(legend.position="none"),
                                             sensitivity_8 + theme(legend.position="none"),
                                             sensitivity_9 + theme(legend.position="none"),
                                             sensitivity_10 + theme(legend.position="none"),
                                             nrow=2), mylegend, ncol=2,widths=c(100, 20))


## see how much variation is caused by changing window width vs changing start year
sens <- sensitivity %>%
  group_by(time_window_width, start_year) %>%
  do(tidy(lm(., formula = spec_exp ~ time_window_start))) %>%
  filter(term == "time_window_start")  

sens_ww <- sens %>%
  group_by(start_year) %>%
  do(mutate(., sd = sd(.$estimate))) %>%
  mutate(parameter_varied = "window_width") %>%
  rename(held_constant = "start_year") %>%
  select(parameter_varied, held_constant, sd) %>%
  unique()

sens_sy <- sens %>%
  group_by(time_window_width) %>%
  do(mutate(., sd = sd(.$estimate)))  %>%
  mutate(parameter_varied = "start_year")  %>%
  rename(held_constant = "time_window_width") %>%
  select(parameter_varied, held_constant, sd) %>%
  unique()

var <- rbind(sens_sy, sens_ww) 

var %>%
  ggplot(., aes(y = sd, x = parameter_varied, colour = held_constant)) + 
  geom_point(position = position_dodge(width = 0.1)) + 
  labs(x = "Sensitivity parameter varied", y = "Standard deviation of change in spectral exponent") +
  scale_x_discrete(labels = c("Start year", "Window width"))
