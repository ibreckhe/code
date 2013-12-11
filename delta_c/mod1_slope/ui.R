
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel(
    HTML(
      '<div id="stats_header">
          <img id="stats_logo" align="right" 
              alt="deltaC: Statistics for environmental change" 
              src="https://dl.dropboxusercontent.com/u/596355/deltac_logo_small.png" />
          <h2>Module 1: When is a trend actually trending?</h2>
          <h3>Autocorrelation and Slope</h3>
        </div>'
    ),
    windowTitle="deltaC - Autocorrelation and Trend Uncertainty"
  ), 
  # Sidebar with instructions
  sidebarPanel(
    h4("This time-series does not have autocorrelation."),
    sliderInput("nomem_center", 
                "Slide or press play to change the regression window", 
                min = 20, 
                max = 134, 
                value = 20,
                step=1,
                animate=animationOptions(interval=200, loop=T)),
    HTML(
      '<br><br><br>'
    ),
    numericInput("window","Width of Regression Window (time steps)",
                 min=4,
                 max=19,
                 value=15),
    HTML(
      '<br><br><br><br><br><br>'
    ),
    h4("This time-series has autocorrelation."),
    sliderInput("mem_center", 
                "Slide or press play to change the regression window", 
                min = 20, 
                max = 134, 
                value = 20,
                step=1,
                animate=animationOptions(interval=200, loop=T))
    ),
  
  # Show a plot of the generated distribution
  mainPanel(
                plotOutput("nomemPlot",height=300),
                plotOutput("memPlot",height=300)
        )
    
  )
)
