### shiny app to visualize a&s model
library(shiny)
library(tidyverse)

fitness_plot <- function(temp_min, temp_max, s, d_j, d_a, Ad_j, Ad_a, Aa, alpha, b_tr, Topt, Tr) {
  ## call function over range of temperatures
  range <- seq(from = temp_min, to = temp_max, by = 0.001)
  temps <- c(273.15 + range)
  
  fitness <- function(temp) {
    TD = ((1/Tr) - (1/temp)) ## 
    bT = b_tr*exp(-((temp-Topt)^2/(2*s^2))) 
    inv_alphaT = alpha*exp(Aa*TD)
    term1 = d_a*exp(Ad_a*TD)
    term2 = (1/(alpha*exp(Aa*TD)))*(LambertW::W(b_tr*alpha*exp(Aa*TD - ((temp-Topt)^2/(2*s^2)) + alpha*exp(Aa*TD)*(d_a*exp(Ad_a*TD) - d_j*exp(Ad_j*TD)))))
    rT = -term1 + term2
    
    return(c(term1, term2, rT, bT, inv_alphaT))
  }
  
  data <- sapply(FUN = fitness, temps)
  
  ## plot results
  data = data.frame(term1 = data[1,], term2 = data[2,], rT = data[3,], temps = temps,
                    bT = data[4,], inv_alphaT = data[5,]) %>%
    select(-bT, inv_alphaT) %>%
    gather(key = "term", value = "value", c(term1, term2, rT)) %>%
    mutate(facet = ifelse(term == "rT", "rT", "terms"))
  
  ymax = max(data$value)
  
   plot <- data %>%
    ggplot(., aes(x = temps, y = value, colour = term)) + geom_point(size = 0.1) + theme_bw() +
    labs(x = "Temperature (K)", y = "rm(T) or components of rm(T)", colour = "") +
    scale_y_continuous(limits = c(0, ymax)) +
    facet_wrap(~facet, nrow = 2)  + 
    #geom_vline(xintercept = 223.15) + geom_vline(xintercept = 313.15) + 
    scale_color_discrete(labels = c("rm(T)", "Term 1", "Term 2")) +
    theme(strip.text.x = element_blank())

  return(plot)
}

K_plot <- function(M, Ea) {
  
  ## vary activation energy of metabolism
  temp_max = 50
  temp_min = -50
  k = 8.617e-5 ## unts: eV K-1
  T=NA
  Kt = M^(-3/4)*exp(Ea/(k*T))
  
  range <- seq(from = temp_min, to = temp_max, by = 0.001)
  temps <- c(273.15 + range)
  

  Kt = M^(-3/4)*exp(Ea/(k*temps))

  Eas <- rep(Ea, each = length(temps))
  df <- data.frame(Kt = Kt, Eas = Ea, temps = temps)
  
  plot <- df %>%
    ggplot(., aes(x = temps, y = Kt)) + geom_point() + theme_bw() +
    labs(x = "Temperature (K)", y = "Carrying capacity")
  
  return(plot)
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Temperature dependence of population growth rate"),
    
    br(),
    
    p("Model by Amarasekare and Savage (2012) - ", em("A framework for elucidating the temperature dependence of fitness.")),
    
    p("Population growth rate varies as a function of temperature because of the temperature dependence of mortality, development, and fecundity."), 
    
    br(),
    
    plotOutput("rPlot", width = "75%"),
   
    br(),
    
    fluidRow(
      column(3,
             h4("Mortality"),
             sliderInput("d_j",
                         "Juvenile mortality at reference temperature (d_j)",
                         min = 0.1,
                         max = 0.2,
                         value = 0.03, step = 0.01),
             sliderInput("d_a",
                         "Adult mortality at reference temperature (d_a)",
                         min = 0.1,
                         max = 0.2,
                         value = 0.05, step = 0.01),
             sliderInput("Ad_j",
                         "Arrhenius constant for juvenile mortality (Ad_j)",
                         min = 5000,
                         max = 20000,
                         value = 7500, step = 10),
             sliderInput("Ad_a",
                         "Arrhenius constant for adult mortality (Ad_a)",
                         min = 5000,
                         max = 20000,
                         value = 10000, step = 10)),
      column(3,
             h4("Development"),
             sliderInput("Aa",
                         "Arrhenius constant for development (Aa)",
                         min = -20000,
                         max = 0,
                         value = -4000, step = 10),
             sliderInput("alpha",
                         "Age at maturity (alpha)",
                         min = 0.5,
                         max = 200,
                         value = 1, step = 0.5)),
      column(3,
             h4("Fecundity"),
             sliderInput("Topt",
                         "Temperature at which avg. fecundity is maximized (Topt)",
                         min = 250,
                         max = 330,
                         value = 273.15, step = 1),
             sliderInput("s",
                         "Variability around maximum fecundity (s)",
                         min = 0.5,
                         max = 20,
                         value = 4.8, step = 0.1),
             sliderInput("b_tr",
                         "Average per capita fecundity at reference temperature (b_tr)",
                         min = 0.5,
                         max = 200,
                         value = 50, step = 0.5),
             sliderInput("Tr",
                         "Reference temperature (Tr)",
                         min = 250,
                         max = 330,
                         value = 298, step = 1))),
    # Application title
    titlePanel("Temperature dependence of carrying capacity"),
    
    br(),
    
    p("Model by Savage et al. (2004) - ", em("Effects of body size and temperature on population growth.")),
    
    p("Carrying capacity varies as a function of temperature and mass according to the metabolic theory of ecology (MTE)."), 
    
    br(),
    
    plotOutput("KPlot", width = "75%"),
    
    fluidRow(
      column(3,
             sliderInput("M",
                         "Average mass of individual (M)",
                         min = 1,
                         max = 200,
                         value = 50, step = 1)
             
      ),
      column(3,
             sliderInput("Ea",
                         "Activation energy of metabolism (Ea)",
                         min = 0.1,
                         max = 0.7,
                         value = 0.2, step = 0.1)
             )
      
    )
      
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$rPlot <- renderPlot({
      
      temp_min = -50
      temp_max = 50
      s = input$s
      d_j = input$d_j # juvenile mortality at reference temperature 
      d_a = input$d_a # adult mortality rate at reference temperature 
      Ad_j = input$Ad_j # Arrhenius constant for juvenile mortality  
      Ad_a = input$Ad_a # Arrhenius constant for adult mortality  
      Aa = input$Aa # Arrhenius constant for development
      alpha = input$alpha # age at maturity
      b_tr = input$b_tr # average per capita fecundity at reference temperature, height and shape of term 2
      Topt = input$Topt # temperature at which avg fecundity is maximal, shifts term 2 
      Tr = input$Tr # reference temperature, changes height of term curves
      
      ## call plot function with parameters
      fitness_plot(temp_min = temp_min,
                   temp_max = temp_max, s = s, d_a = d_a, Ad_j = Ad_j, 
                   Ad_a = Ad_a, Aa = Aa, alpha = alpha, b_tr = b_tr, 
                   Topt = Topt, Tr = Tr, d_j = d_j)
        
    })
    
    output$KPlot <- renderPlot({ 
     
      M = input$M
      Ea = input$Ea ## units: eV
      
      K_plot(M = M, Ea = Ea)
      })
}

# Run the application 
shinyApp(ui = ui, server = server)


