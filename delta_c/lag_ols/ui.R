
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Lagged Regression"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    dateRangeInput("date_range", "Date Range of Analysis:",start="2004-01-01",end="2013-7-28",
                   min="2004-01-01",max="2013-7-28",startview="month"),
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    plotOutput("timePlot")
  )
))
