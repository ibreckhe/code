
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

##Generates two time-series: one with memory, one without.
set.seed(42)
int <- 0.5
slope <- 0.03
n <- 150
mem_sim <- int+slope*seq(1:n)+arima.sim(n = n, 
                    model=list(ar = c(1, -0.2), ma = c(-0.2279, 0.2488)),
                    rand.gen = function(n, ...) sqrt(0.1796) * rt(n, df = 5))
nomem_sim <- int+slope*seq(1:n)+arima.sim(n=n,model=list())

##Function to subset the data,fits a regression to the data, and plot it.

plotfun <- function(series,window,center){ 
  ##subsets and computes the regression.
  start <- center-window
  end <- center+window
  sub <- series[start:end]
  sub_lm <- lm(sub~c(start:end))
  ##creates the plot
  ymin <- min(sub) - 20
  ymax <- max(sub) + 20
  ts.plot(series,col=rgb(0.3,0.3,0.8,1),lwd=1.5)
  polygon(x=c(start,end,rev(c(start,end))),
          y=c(rep(ymin,2),rep(ymax,2)),
          col=rgb(0.5,0.5,0.5,0.5),border=NA)
  abline(sub_lm,lty=2,lwd=1.5)
  points(c(start:end),sub,pch=21,bg="white")
}

shinyServer(function(input, output) {
   
  output$nomemPlot <- renderPlot({
    ##Plots the no-memory time series.
    plotfun(nomem_sim,input$window,input$nomem_center)
  })

  output$memPlot <- renderPlot({
    ##Plots the no-memory time series.
    plotfun(mem_sim,input$window,input$mem_center)
  })
  
})
