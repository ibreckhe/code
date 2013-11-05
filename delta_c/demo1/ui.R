library(shiny)

# Define UI for application that plots google searches for cats over time.
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel(
    HTML(
      '<div id="stats_header">
          <img id="stats_logo" align="right" 
              alt="deltaC: Statistics for environmental change" 
              src="https://dl.dropboxusercontent.com/u/596355/deltac_logo_small.png" />
          <h2>Tech Demo #1</h2>
          <h3>When is a trend actually trending?</h3>
        </div>'
    ),
    windowTitle="deltaC Tech Demo #1"
  ),
    
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    
    #Search Term
    textInput("term", "Google Ngrams Search Term(s):", value = "climate change"),
    
    # Specification of range within an interval
    numericInput("start","From Year:",value=1980),
    numericInput("end","To Year:",value=2012),
    submitButton("Update View")
  ),
  
  # Show a plot of the generated timeseries.
  mainPanel(
    plotOutput("timePlot")
  )
))