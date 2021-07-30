library(shiny)
library(tidyverse)
library(raster)

setwd("/Users/nikkimoore/Documents/velocities-of-variability")

## read in raster stack of temps on may 28:
maps = stack("data-processed/mosaic.grd")
l_maps = stack("data-processed/l_mosaic.grd")
s_maps = stack("data-processed/s_mosaic.grd")

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Visualizing global climate models"),
    
    # Sidebar with a slider input for which year 
    sidebarLayout(
        
        # Show a map of air surface temperature on May 28
        mainPanel(
            plotOutput("normalPlot"),
            plotOutput("lPlot"),
            plotOutput("sPlot")
        ),
        sidebarPanel(
            sliderInput("year",
                        "Year:",
                        min = 1871,
                        max = 2100,
                        value = 1)
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$normalPlot <- renderPlot({
        # grab map data based on year input from ui.R
        m <- maps[[input$year - 1870]]
        
        xyz <- data.frame(rasterToPoints(m))
        colnames(xyz) <- c('x', 'y', 'z')
        
        # draw the map
        ggplot(xyz, aes(x = x, y = y, fill = z - 273.15)) + geom_raster() + coord_fixed() + 
            theme_void() +
            labs(x = "Longitude", y = "Latitude", fill = "Temperature (C)",
                 main = "Air surface temperate on May 28") +
            scale_fill_gradientn(limits = c(-80, 60),
                                 colours=c("darkblue", "white", "darkred"))
    })
    
    output$lPlot <- renderPlot({
        # grab map data based on year input from ui.R
        m <- l_maps[[input$year - 1870]]
        
        xyz <- data.frame(rasterToPoints(m))
        colnames(xyz) <- c('x', 'y', 'z')
        
        # draw the map
        ggplot(xyz, aes(x = x, y = y, fill = z)) + geom_raster() + coord_fixed() + 
            theme_void() +
            labs(x = "Longitude", y = "Latitude", fill = "Temperature (C)",
                 main = "Linearly detrended air surface temperate on May 28") +
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
            theme_void() +
            labs(x = "Longitude", y = "Latitude", fill = "Temperature (C)",
                 main = "Seasonally detrended air surface temperate on May 28") +
            scale_fill_gradientn(limits = c(-40, 60),
                                 colours=c("darkblue", "white", "darkred"))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)


