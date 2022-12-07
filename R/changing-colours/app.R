library(ggplot2)
library(shiny)
library(fractal)
library(tidyverse)
library(broom)

# start_colour = function(){return(0.9)}
# end_colour = function(){return(1)}
# num_steps = function(){return(110)}
# by = function() {
#   return((end_colour() - start_colour())/num_steps())
#   }
# each = function() {
#   return(ceiling(200000/num_steps()))
# }
# noise <- function() {
#   by = (end_colour() - start_colour())/num_steps()
#   each = ceiling(200000/num_steps())
# 
#   if(by == 0) {
#     alpha = rep(start_colour(), 200000)
#     delta = -alpha/-2
#   }
#   else {
#     ## simulate the time series
#     alpha <- seq(start_colour(), end_colour(), by = by)
#     alpha = rep(alpha, each = each)[1:200000]
#     delta = -alpha/-2
#   }
# 
#   ## set the innovations variance to unity
#   innovation <- rep(1, length(delta))
# 
#   ## simulate a time-varying FD process
#   noise <- FDSimulate(delta = delta, innovation = innovation)
# 
#   return(noise)
# }
# correct_noise <- function() {
#   noise_data = noise()
# 
#   data = data.frame(Noise = c(noise_data), "Time" = 1:200000)
# 
#   ## fix it so that each time step begins at last time step value:
#   breaks = seq(1, 220000, by = each())
#   new_noise = data$Noise
# 
#   i = 2
#   while(i < length(breaks)) {
#     to_add = new_noise[breaks[i] - 1] - new_noise[breaks[i]]
# 
#     new_noise[breaks[i]:(breaks[i+1] - 1)] = new_noise[breaks[i]:(breaks[i+1] - 1)] + to_add
#     i = i + 1
#   }
# 
#   return(new_noise)
# }
# correct_noise_conditional <- function() {
#   noise_data = noise()
# 
#   alpha <- seq(start_colour(), end_colour(), by = by())
# 
#   if(any(alpha < 1) & !all(alpha < 1)) {
#     inc = start_colour() < end_colour()
# 
#     ## define point at which alpha switches to greater/less than 1
#     if(inc == TRUE) {
#       point = last(which(alpha < 1))
#     }
#     else if (inc == FALSE) {
#       point = first(which(alpha < 1))
#     }
# 
#     data = data.frame(Noise = c(noise_data), "Time" = 1:200000)
# 
#     ## fix it so that each time step begins at last time step value IF alpha < 1:
#     breaks = seq(1, 220000, by = each())
#     new_noise = data$Noise
# 
#     i = 2
#     while(i < length(breaks)) {
#       if((i > point & inc == TRUE) | (i < point & inc == FALSE)) {
#         to_add = new_noise[breaks[i] - 1] - new_noise[breaks[i]]
# 
#         new_noise[breaks[i]:(breaks[i+1] - 1)] = new_noise[breaks[i]:(breaks[i+1] - 1)] + to_add
# 
#         i = i + 1
#       }
#       else {
#         i = i + 1
#       }
#     }
#   }
# 
#   else if (any(alpha > 1)) {
#     data = data.frame(Noise = c(noise_data), "Time" = 1:200000)
# 
#     ## fix it so that each time step begins at last time step value:
#     breaks = seq(1, 220000, by = each())
#     new_noise = data$Noise
# 
#     i = 2
#     while(i < length(breaks)) {
#       to_add = new_noise[breaks[i] - 1] - new_noise[breaks[i]]
# 
#       new_noise[breaks[i]:(breaks[i+1] - 1)] = new_noise[breaks[i]:(breaks[i+1] - 1)] + to_add
# 
#       i = i + 1
#     }
#   }
#   else {
#     new_noise = c(noise_data)
#   }
#   return(new_noise)
# }
# 
# plot_spectral_change <- function() {
# 
#   ## get simulated noise
#   noise_data = noise()
#   corr_noise = correct_noise()
#   corr_noise_cond_data = correct_noise_conditional()
# 
#   ## make dataframe
#   data = data.frame(noise_data = c(noise_data)[1:200000],
#                     corr_noise = corr_noise[1:200000],
#                     corr_noise_cond_data = corr_noise_cond_data[1:200000],
#                     time = 1:200000)
# 
#   ## change NA to 0
#   data[which(is.na(data))] <- 0
#   element = 1
#   spec_exp_list <- list()
# 
#   n = 5
#   while (n < 11) {
#     year_start <- 1
#     year_stop <- 1 + n*365
# 
#     while (year_start <= (150000 - n*365)) {
#       ## extract temps within time window
#       ts_chunk <- data[year_start:year_stop,c(1,2,3)]
# 
#       ## get length
#       L = nrow(ts_chunk)
# 
#       ## preprocess the time series:
#       ## a. subtracting mean
#       ts_n <- ts_chunk$noise_data - mean(ts_chunk$noise_data)
#       ts_cn <- ts_chunk$corr_noise - mean(ts_chunk$corr_noise)
#       ts_cnc <- ts_chunk$corr_noise_cond_data - mean(ts_chunk$corr_noise_cond_data)
# 
#       ## b. windowing - multiply by a parabolic window
#       window <- parabolic_window(series = ts_n, N = L)
#       ts_n <- ts_n*window
#       window <- parabolic_window(series = ts_cn, N = L)
#       ts_cn <- ts_cn*window
#       window <- parabolic_window(series = ts_cnc, N = L)
#       ts_cnc <- ts_cnc*window
# 
#       ## c. bridge detrending (endmatching)
#       ## ie. subtracting from the data the line connecting the first and last points of the series
#       ts_n <- bridge_detrender(windowed_series = ts_n, N = L)
#       ts_cn <- bridge_detrender(windowed_series = ts_cn, N = L)
#       ts_cnc <- bridge_detrender(windowed_series = ts_cnc, N = L)
# 
#       ## calculate spectral exponent in window using PSD and AWC methods
#       exp_PSD_n <- spectral_exponent_calculator_PSD(ts_n, l = L)
#       exp_PSD_cn <- spectral_exponent_calculator_PSD(ts_cn, l = L)
#       exp_PSD_cnc <- spectral_exponent_calculator_PSD(ts_cnc, l = L)
# 
#       exp_PSD_low_n <- exp_PSD_n[[1]]
#       exp_PSD_all_n <- exp_PSD_n[[3]]
#       exp_PSD_low_cn <- exp_PSD_cn[[1]]
#       exp_PSD_all_cn <- exp_PSD_cn[[3]]
#       exp_PSD_low_cnc <- exp_PSD_cnc[[1]]
#       exp_PSD_all_cnc <- exp_PSD_cnc[[3]]
# 
#       spec_exp_list[[element]] <- c(exp_PSD_low_n, exp_PSD_all_n,
#                                     exp_PSD_low_cn, exp_PSD_all_cn,
#                                     exp_PSD_low_cnc, exp_PSD_all_cnc,
#                                     year_start, year_stop,
#                                     paste(n, "years"))
# 
# 
#       ## move to next window
#       year_start = year_stop + 1
#       year_stop = year_stop + n*365
# 
#       element = element + 1
#     }
# 
#     ## move to next window width
#     n = n + 1
#   }
# 
#   ## bind rows in list into data frame
#   spec_exp_df <- data.frame(do.call(rbind, spec_exp_list), stringsAsFactors = FALSE)
#   colnames(spec_exp_df) <- c("exp_PSD_low_n", "exp_PSD_all_n",
#                              "exp_PSD_low_cn", "exp_PSD_all_cn",
#                              "exp_PSD_low_cnc", "exp_PSD_all_cnc",
#                              "window_start_year","window_stop_year", "time_window_width")
#   spec_exp_df <- spec_exp_df %>%
#     gather(key = "Time series type", value = "Measured spectral exponent",
#            c(exp_PSD_low_n, exp_PSD_all_n,
#              exp_PSD_low_cn, exp_PSD_all_cn,
#              exp_PSD_low_cnc, exp_PSD_all_cnc)) %>%
#     mutate(ts_type = ifelse(`Time series type` %in% c("exp_PSD_low_n", "exp_PSD_all_n"),
#                             "original simulation",
#                             ifelse(`Time series type` %in% c("exp_PSD_low_cn", "exp_PSD_all_cn"),
#                                    "corrected simulation",
#                                    "conditionally-corrected simulation")),
#            se_type = ifelse(`Time series type` %in% c("exp_PSD_low_n", "exp_PSD_low_cn",
#                                                       "exp_PSD_low_cnc"),
#                             "low", "all"))
# 
#   ## convert numbers to numeric
#   spec_exp_df[,c(1,2,5)] <- sapply(spec_exp_df[,c(1,2,5)], as.numeric)
# 
#   return(spec_exp_df)
# }


 
# Define UI for application that draws a histogram
ui <- fluidPage(
  
    # Application title
    titlePanel("Simulating noise with changing colour"),
    
    br(),
    
    fluidRow(
      column(10, align="left", offset = 0,
             textOutput("ask_input")
      )), 
    
    br(),
    
    # take input
    fluidRow(
      column(6, numericInput("baseline_colour", 
                             label = "Starting spectral exponent:", 
                             min = 0.001, 
                             max = 3.0, 
                             step = 0.001,
                             value = 1,
                             width = "auto"),
             numericInput("ending_colour", 
                          label = "Ending spectral exponent:", 
                          min = 0.001, 
                          max = 3.0, 
                          step = 0.001,
                          value = 2,
                          width = "auto")),
      column(6, numericInput("num_steps", 
                             label = "Number of steps", 
                             min = 10000, 
                             max = 1, 
                             step = 1,
                             value = 10,
                             width = "auto"),
             textOutput("window_width"))
    ),

    br(),
    br(),
    
    fluidRow(
      column(12, align="center", offset = 0,
             actionButton("vis", "Click to visualize parameters")
      )),
    
    # plots
    plotOutput("mainPlot", height = "200px"),
    
    br(),
    br(),
    
    fluidRow(
      column(12, align="center", offset = 0,
             actionButton("sim", "Click to simulate noise")
      )),
   
    plotOutput("secondPlot", height = "200px"),
    plotOutput("thirdPlot",  height = "200px"),
    plotOutput("fifthPlot", height = "200px"),
    
    br(),
    br(),
    
    fluidRow(
      column(12, align="center", offset = 0,
             actionButton("spec", "Click to measure spectral change")
      )),
    plotOutput("spec_plot", height = "600px"),
    plotOutput("spec_lms", height = "400px"),
    textOutput("plot_desc"),
    
    br(),
    br()
    
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$ask_input = renderText({
    "Select a starting spectral exponent, an ending spectral exponent, and the number of steps it will take to go from one to the other."
  })
  
  ## print out time step width
  output$window_width <- renderText({
    paste("Time step width: ~", round(ceiling(200000/input$num_steps)/365, 2), " years", sep = "")
  })
  
  ## get input parameters
  start_colour = reactive({
    return(input$baseline_colour)
  })
  end_colour = reactive({
    return(input$ending_colour)
  })
  num_steps = reactive({
    return(input$num_steps)
  })
  by = reactive({
    return((end_colour() - start_colour())/num_steps())
  })
  each = reactive({
    return(ceiling(200000/num_steps()))
  })
  
  ## simulate noise based on input
  noise <- reactive({
    by = (end_colour() - start_colour())/num_steps()
    each = ceiling(200000/num_steps())
    
    if(by == 0) {
      alpha = rep(start_colour(), 200000)
      delta = -alpha/-2
    }
    else {
      ## simulate the time series
      alpha <- seq(start_colour(), end_colour(), by = by)
      alpha = rep(alpha, each = each)[1:200000]
      delta = -alpha/-2
    }
    
    ## set the innovations variance to unity
    innovation <- rep(1, length(delta))
    
    ## simulate a time-varying FD process
    noise <- FDSimulate(delta = delta, innovation = innovation)
    
    return(noise)
  })
 
  ## correct the noise to make less jumpy
  correct_noise <- reactive({
    noise_data = noise()
    
    data = data.frame(Noise = c(noise_data), "Time" = 1:200000)
    
    ## fix it so that each time step begins at last time step value:
    breaks = seq(1, 220000, by = each())
    new_noise = data$Noise
    
    i = 2
    while(i < length(breaks)) {
      to_add = new_noise[breaks[i] - 1] - data$Noise[breaks[i]]
      
      new_noise[breaks[i]:(breaks[i+1] - 1)] = new_noise[breaks[i]:(breaks[i+1] - 1)] + to_add
      i = i + 1
    }
    
    return(new_noise)
  })
  
  correct_noise_conditional <- reactive ({
    noise_data = noise()
    
    alpha <- seq(start_colour(), end_colour(), by = by())
    
    if(any(alpha < 1) & !all(alpha < 1)) {
      inc = start_colour() < end_colour()
      
      ## define point at which alpha switches to greater/less than 1
      if(inc == TRUE) {
        point = last(which(alpha < 1))
      }
      else if (inc == FALSE) {
        point = first(which(alpha < 1))
      }
      
      data = data.frame(Noise = c(noise_data), "Time" = 1:200000)
      
      ## fix it so that each time step begins at last time step value IF alpha < 1:
      breaks = seq(1, 220000, by = each())
      new_noise = data$Noise
      
      i = 2
      while(i < length(breaks)) {
        if((i > point & inc == TRUE) | (i < point & inc == FALSE)) {
          to_add = new_noise[breaks[i] - 1] - new_noise[breaks[i]]
          
          new_noise[breaks[i]:(breaks[i+1] - 1)] = new_noise[breaks[i]:(breaks[i+1] - 1)] + to_add
          
          i = i + 1
        }
        else {
          i = i + 1
        }
      }
    }
    
    else if (any(alpha > 1)) {
      data = data.frame(Noise = c(noise_data), "Time" = 1:200000)
      
      ## fix it so that each time step begins at last time step value:
      breaks = seq(1, 220000, by = each())
      new_noise = data$Noise
      
      i = 2
      while(i < length(breaks)) {
        to_add = new_noise[breaks[i] - 1] - new_noise[breaks[i]]
        
        new_noise[breaks[i]:(breaks[i+1] - 1)] = new_noise[breaks[i]:(breaks[i+1] - 1)] + to_add
        
        i = i + 1
      }
    }
    else {
      new_noise = c(noise_data)
    }
    return(new_noise)
  })
  
  spectral_change <- reactive({
    ## get simulated noise
    noise_data = noise()
    corr_noise = correct_noise()
    corr_noise_cond_data = correct_noise_conditional()
    
    ## make dataframe 
    data = data.frame(noise_data = c(noise_data)[1:200000], 
                      corr_noise = corr_noise[1:200000],
                      corr_noise_cond_data = corr_noise_cond_data[1:200000],
                      time = 1:200000)
    
    ## change NA to 0
    data[which(is.na(data))] <- 0
    element = 1
    spec_exp_list <- list()
    
    n = 5
    while (n < 11) {
      year_start <- 1
      year_stop <- 1 + n*365
      
      while (year_start <= (150000 - n*365)) {
        ## extract temps within time window
        ts_chunk <- data[year_start:year_stop,c(1,2,3)]
        
        ## get length
        L = nrow(ts_chunk)
        
        ## preprocess the time series:
        ## a. subtracting mean
        ts_n <- ts_chunk$noise_data - mean(ts_chunk$noise_data)
        ts_cn <- ts_chunk$corr_noise - mean(ts_chunk$corr_noise)
        ts_cnc <- ts_chunk$corr_noise_cond_data - mean(ts_chunk$corr_noise_cond_data)
        
        ## b. windowing - multiply by a parabolic window 
        window <- parabolic_window(series = ts_n, N = L)
        ts_n <- ts_n*window
        window <- parabolic_window(series = ts_cn, N = L)
        ts_cn <- ts_cn*window
        window <- parabolic_window(series = ts_cnc, N = L)
        ts_cnc <- ts_cnc*window
        
        ## c. bridge detrending (endmatching)
        ## ie. subtracting from the data the line connecting the first and last points of the series
        ts_n <- bridge_detrender(windowed_series = ts_n, N = L)
        ts_cn <- bridge_detrender(windowed_series = ts_cn, N = L)
        ts_cnc <- bridge_detrender(windowed_series = ts_cnc, N = L)
        
        ## calculate spectral exponent in window using PSD and AWC methods
        exp_PSD_n <- spectral_exponent_calculator_PSD(ts_n, l = L)
        exp_PSD_cn <- spectral_exponent_calculator_PSD(ts_cn, l = L)
        exp_PSD_cnc <- spectral_exponent_calculator_PSD(ts_cnc, l = L)
        
        exp_PSD_low_n <- exp_PSD_n[[1]]
        exp_PSD_all_n <- exp_PSD_n[[3]]
        exp_PSD_low_cn <- exp_PSD_cn[[1]]
        exp_PSD_all_cn <- exp_PSD_cn[[3]]
        exp_PSD_low_cnc <- exp_PSD_cnc[[1]]
        exp_PSD_all_cnc <- exp_PSD_cnc[[3]]
        
        spec_exp_list[[element]] <- c(exp_PSD_low_n, exp_PSD_all_n,
                                      exp_PSD_low_cn, exp_PSD_all_cn,
                                      exp_PSD_low_cnc, exp_PSD_all_cnc,
                                      year_start, year_stop, 
                                      paste(n, "years"))
        
        
        ## move to next window
        year_start = year_stop + 1
        year_stop = year_stop + n*365 
        
        element = element + 1
      }
      
      ## move to next window width
      n = n + 1
    }
    
    ## bind rows in list into data frame
    spec_exp_df <- data.frame(do.call(rbind, spec_exp_list), stringsAsFactors = FALSE)
    colnames(spec_exp_df) <- c("exp_PSD_low_n", "exp_PSD_all_n",
                               "exp_PSD_low_cn", "exp_PSD_all_cn",
                               "exp_PSD_low_cnc", "exp_PSD_all_cnc",
                               "window_start_year","window_stop_year", "time_window_width")
    spec_exp_df <- spec_exp_df %>%
      gather(key = "Time series type", value = "Measured spectral exponent",
             c(exp_PSD_low_n, exp_PSD_all_n,
               exp_PSD_low_cn, exp_PSD_all_cn,
               exp_PSD_low_cnc, exp_PSD_all_cnc)) %>%
      mutate(ts_type = ifelse(`Time series type` %in% c("exp_PSD_low_n", "exp_PSD_all_n"),
                              "original simulation",
                              ifelse(`Time series type` %in% c("exp_PSD_low_cn", "exp_PSD_all_cn"),
                                     "corrected simulation",
                                     "conditionally-corrected simulation")),
             se_type = ifelse(`Time series type` %in% c("exp_PSD_low_n", "exp_PSD_low_cn",
                                                        "exp_PSD_low_cnc"),
                              "low", "all"))
    
    ## convert numbers to numeric
    spec_exp_df[,c(1,2,5)] <- sapply(spec_exp_df[,c(1,2,5)], as.numeric)
    
    return(spec_exp_df)
  })
  
  observeEvent(input$vis, {
    ## plot spectral exponent over time according to input parameters
    
    param_plot = eventReactive(input$vis, {
      by = (end_colour() - start_colour())/num_steps()
      each = each()
      
      if(by == 0) {
        alpha = rep(start_colour(), 200000)
        delta = -alpha/-2
      }
      else {
        ## simulate the time series
        alpha <- seq(start_colour(), end_colour(), by = by)
        alpha = rep(alpha, each = each)[1:200000]
        delta = -alpha/-2
      }
      
      data = data.frame("Spectral exponent" = alpha, "Time" = 1:200000)
      
      ggplot(data = data, aes(x = Time, y = Spectral.exponent)) + geom_line() + theme_bw() +
        labs(y = "Spectral exponent", title = "Change in spectral exponent:") + 
        scale_y_continuous(limits = c(0, 3)) +
        theme(panel.grid = element_blank())
    })
      
    output$mainPlot <- renderPlot({
      param_plot()
    })
  })
  
    observeEvent(input$sim, {
      
      output$corrected_noise = renderText({
        "Simulated noise after correcting so that there are no huge jumps between steps:"
      })
      
      output$mv_noise = renderText({
        "Corrected noise after removing mean and standardizing variance:"
      })
      
      second_plot = eventReactive(input$sim, {
        ## plot simulated noise 
        noise_data = noise()
          
        data = data.frame(Noise = c(noise_data), "Time" = 1:200000)
          
        ggplot(data = data, aes(x = Time, y = Noise)) + geom_line(linewidth = 0.1) + theme_bw() +
            theme(panel.grid = element_blank()) +
          labs(title = "Original simulation from FDSimulate()")
      })
      
      third_plot = eventReactive(input$sim, {
        ## plot corrected noise
        corrected_noise = correct_noise()
          
        new_data = data.frame(corrected_noise = corrected_noise[1:200000], "Time" = 1:200000)
          
        ggplot(data = new_data, aes(x = Time, y = corrected_noise)) + geom_line(linewidth = 0.1) + theme_bw() +
            theme(panel.grid = element_blank()) + 
            labs(y = "Corrected noise", title = "Corrected simulation")
      })
      
      fifth_plot = eventReactive(input$sim, {
        ## plot corrected noise
        new_correction = correct_noise_conditional()
        
        ## remove mean and change variance to 4140:
        new_noise <- new_correction[1:200000]*1/sqrt(var(new_correction[1:200000]))*sqrt(4140)
        new_noise <- new_noise - mean(new_noise)
        
        new_data = data.frame(new_noise = new_noise, "Time" = 1:200000)
        
        ggplot(data = new_data, aes(x = Time, y = new_noise)) + geom_line(linewidth = 0.1) + theme_bw() +
          theme(panel.grid = element_blank()) + 
          labs(y = "Conditionally-corrected noise", title = "Conditionally-corrected simulation")
      })
      
      output$secondPlot <- renderPlot({
        second_plot()
      })
      output$thirdPlot <- renderPlot({
        third_plot()
      })
      output$fifthPlot <- renderPlot({
        fifth_plot()
      })
      
    })
    
    observeEvent(input$spec, {
      spec_plot = eventReactive(input$spec, {
        ## calculate spectral change
        spec_exp_df = spectral_change()
        
        exp_slope = (end_colour() - start_colour())/200000
        
        ## plot
       spec_exp_df %>%
          mutate(time_window_width = factor(.$time_window_width, 
                                            levels = c("5 years", "6 years", "7 years", 
                                                       "8 years", "9 years", "10 years"), 
                                            ordered = TRUE)) %>%
          mutate(ts_type = factor(.$ts_type,  levels = c("original simulation", "corrected simulation",
                                                       "conditionally-corrected simulation"), 
                                            ordered = TRUE)) %>%
          filter(se_type != "all") %>%
          ggplot(aes(x = window_start_year, y = `Measured spectral exponent`, colour = ts_type)) +
          geom_abline(slope = exp_slope, intercept = start_colour(), colour = "black") +
          geom_point() +
          labs(x = "Time", y = "Measured spectral exponent", colour = "") +
          facet_wrap(~time_window_width) +
          geom_smooth(method = "lm", se = FALSE) + theme_bw() +
          theme(panel.grid = element_blank(), legend.position = "top") 
        
      })
      
      spec_lms = eventReactive(input$spec, {
        ## calculate spectral change
        spec_exp_df = spectral_change()
        
        exp_slope = (end_colour() - start_colour())/200000
        
        ## fit linear models and plot different slope estimates beside expected slope
        spec_exp_df %>%
          mutate(ts_type = factor(.$ts_type,  levels = c("original simulation", "corrected simulation",
                                                         "conditionally-corrected simulation"), 
                                  ordered = TRUE)) %>%
          mutate(time_window_width = factor(.$time_window_width, 
                                            levels = c("5 years", "6 years", "7 years", 
                                                       "8 years", "9 years", "10 years"), 
                                            ordered = TRUE)) %>%
          filter(se_type != "all") %>%
          group_by(time_window_width, ts_type) %>%
          do(tidy(lm(`Measured spectral exponent` ~ window_start_year, data = .), conf.int = TRUE)) %>% 
          filter(term == "window_start_year") %>% 
          ungroup() %>%
          ggplot(aes(x = time_window_width, colour = ts_type, y = estimate)) + geom_point() +
          geom_point(y = exp_slope, colour = "black", shape = 3) + theme_bw() +
          theme(panel.grid = element_blank(), legend.position = "top") + 
          scale_y_continuous(limits = c(exp_slope - 0.000005, exp_slope + 0.000005)) +
          labs(x = "Time window width", colour = "", y = "Estimated change in spectral exponent")
      })
      
      output$spec_plot <- renderPlot({
        spec_plot()
      })
      
      output$spec_lms <- renderPlot({
        spec_lms()
      })
      
      output$plot_desc <- renderText({
        "*Black line and black plus represent the theoretical change in spectral exponent based on the simulation input parameters."
      })
      
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
