library(shiny)
library(tidyverse)
library(raster)
library(shinyWidgets)
library(shinythemes)

## read in raster stack of temps on may 28:
l_maps = readRDS("l_mosaic-02_CMCC-CM.rds")
s_maps = readRDS("s_mosaic-02_CMCC-CM.rds")

# Define UI for application that draws a histogram
ui <- fluidPage(
    

    # Application title
    titlePanel("Visualizing global climate models"),
   
    theme = shinytheme("paper"), 
    
    fluidRow(column(6, offset = 4,
                    sliderTextInput("year","Year:" , 
                                    choices =  as.character(1871:2100)))),
    fluidRow(plotOutput("lPlot")),
    fluidRow(plotOutput("sPlot"))
        
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$lPlot <- renderPlot({
        # grab map data based on year input from ui.R
        m <- l_maps[[input$year - 1870]]
        
        xyz <- data.frame(rasterToPoints(m))
        colnames(xyz) <- c('x', 'y', 'z')
        
        # draw the map
        ggplot(xyz, aes(x = x, y = y, fill = z)) + geom_raster() + coord_fixed() + 
            theme_minimal() + 
            theme(panel.grid = element_blank()) +
            labs(x = "Longitude", y = "Latitude", fill = "Temperature (C)",
                 title = "Linearly detrended air surface temperate on May 28") +
            scale_fill_gradientn(limits = c(-40, 40),
                                 colours=c("darkblue", "white", "darkred"))
    })
    
    output$sPlot <- renderPlot({
        # grab map data based on year input from ui.R
        m <- s_maps[[input$year - 1870]]
        
        xyz <- data.frame(rasterToPoints(m))
        colnames(xyz) <- c('x', 'y', 'z')
        
        # draw the map
        ggplot(xyz, aes(x = x, y = y, fill = z)) + geom_raster() + coord_fixed() + 
            theme_minimal() + 
            theme(panel.grid = element_blank()) +
            labs(x = "Longitude", y = "Latitude", fill = "Temperature (C)",
                 title = "Seasonally detrended air surface temperate on May 28") +
            scale_fill_gradientn(limits = c(-40, 60),
                                 colours=c("darkblue", "white", "darkred"))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)


