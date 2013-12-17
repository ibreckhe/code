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
          <h2>Module 1: When is a trend actually trending?</h2>
          <h3>Autocorrelation in Regression</h3>
        </div>'
    ),
    windowTitle="deltaC Module #1: When is a trend actually trending?"
  ),
    
  # Sidebar with a slider input for number of observations
  sidebarPanel(
      #Help text
      helpText("Explore different methods for dealing with temporal autocorrelation. Enter a search term below to get started."),
      
      #Search Term
      textInput("term", "Google Ngrams Search Term(s):", value = "climate"),
      
      # Specification of range within an interval
      numericInput("start","From Year:",value=1950),
      numericInput("end","To Year:",value=2010),
      
      # Options to display trend and diagnostic plots.
      checkboxInput("trend","Show Linear Trend",value=F),
      checkboxInput("diag","Show Diagnostic Plots",value=F),
      checkboxInput("mod","Show Model Output",value=F),
      
      # Displays model output in sidebar
      conditionalPanel(condition="input.conditionedPanels=='OLS'",
                       tableOutput("ols_modTable")),
      conditionalPanel(condition="input.conditionedPanels=='Subsampling'",
                       tableOutput("sub_modTable")),
      conditionalPanel(condition="input.conditionedPanels=='Lagging'",
                       tableOutput("lag_modTable")),
    
    # Displays options for subsampling.
    conditionalPanel(condition="input.conditionedPanels=='Subsampling'",
           wellPanel(
             helpText("Subsampling Options:"),
             sliderInput("sub","Subsampling Interval (Years)",
                         min=0,
                         max=5,
                         value=1))
           ),
    conditionalPanel(condition="input.conditionedPanels=='Lagging'",
                     wellPanel(
                       helpText("Lagging Options:"),
                       sliderInput("lag","Interval of Lag (Years)",
                                   min=0,
                                   max=5,
                                   value=1))
    )
  ),
  
  # Draw plots of the generated timeseries.
  mainPanel(
    tabsetPanel(id = "conditionedPanels",
      tabPanel("OLS",    
               plotOutput("ols_timePlot",height=300),
               plotOutput("ols_diagPlot",height=300)),
      tabPanel("Subsampling",
               plotOutput("sub_timePlot",height=300),
               plotOutput("sub_diagPlot",height=300)),
      tabPanel("Lagging",
               plotOutput("lag_timePlot",height=300),
               plotOutput("lag_diagPlot",height=300)),
      tabPanel("GLS",
               helpText("Nothing Here Yet"))
    )
  )
))