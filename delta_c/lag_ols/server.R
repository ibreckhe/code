
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

##Shiny application to explore lagged regression based on ordinary least squares.
library(shiny)

##Reads in the data and converts it to a date format.
data <- read.csv('report.csv',skip=4,nrows=500)
data$Week <- as.character(data$Week)
weeks <- matrix(unlist(strsplit(data$Week,split=" - "),recursive=F),byrow=T,ncol=2)
data <- data.frame(data,start=weeks[,1],end=weeks[,2])
data$start <- as.Date(as.character(data$start))
data$end <- as.Date(as.character(data$end))


shinyServer(function(input, output) {
   
  output$distPlot <- renderPlot({
     
    ##Subsets the data.
    all_dat <- data[data$start>=input$date_range[1]&data$start<=input$date_range[2],]
    
    ##Fits the lagged regression
    model <- lm(cats~start,data=all_dat)
    
    ##Creates a plot
    plot(cats~start,data=all_dat,type='n',xlim=range(all_dat$start),
         xlab="Time",ylab="Search Freqency",ylim=c(min(prd[,2]),max(prd[,3])))
    
    polygon(c(all_dat$start,rev(all_dat$start)),c(prd[,2],rev(prd[,3])), col="grey80",border=NULL,lty="blank")
    polygon(c(all_dat$start,rev(all_dat$start)),c(cnf[,2],rev(cnf[,3])), col="grey60",border=NULL,lty="blank")
    lines(all_dat$start,all_dat$cats,type='l',col=rgb(0,0,1,0.8),lwd=1.5)
    abline(v=input$mid_date)
    lines(all_dat$start,prd[,1],lty=2)
    
  })
  
})
