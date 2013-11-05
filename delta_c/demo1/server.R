library(shiny)
library(ngramr)

# Define server logic required to generate and trend data from google trends
shinyServer(function(input, output) {
  ##Gets data from Google
  data <- reactive({ngram(as.character(input$term),corpus="eng_2012",year_end=2012,smoothing=0)})
  
  ##Takes inputs from the date range input.
  year_start <- reactive({input$start})
  year_end <- reactive({input$end})
  
  output$timePlot <- renderPlot({
    ##Subsets the data.
    data2 <- data()[data()$Year>=year_start() & data()$Year<=year_end(),]
    
    ##Plots the data.
    plot(data2$Year,data2$Frequency,type='b',col="slateblue",xlab="Time",ylab="Ngram Frequency",main="")
    
    ##Adds a trend line to the plot.
    model <- lm(data2$Frequency~data2$Year)
    abline(model,lty=2)
    legend("topleft",legend=paste("Trend: ",format(model$coefficients[2],digits=4),
                                  ", OLS p-value: ",format(summary(model)$coefficients[2,4], digits=4)),bty="n")
    
  })
})