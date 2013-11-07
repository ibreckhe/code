##Shiny application to explore regression predictions based on ordinary least squares.
library(shiny)

##Reads in the data and converts it to a date format.
data <- read.csv('report.csv',skip=4,nrows=500)
data$Week <- as.character(data$Week)
weeks <- matrix(unlist(strsplit(data$Week,split=" - "),recursive=F),byrow=T,ncol=2)
data <- data.frame(data,start=weeks[,1],end=weeks[,2])
data$start <- as.Date(as.character(data$start))
data$end <- as.Date(as.character(data$end))

shinyServer(function(input, output) {
  
#   observe({
#     # We'll use the input$controller variable multiple times, so save it as x
#     # for convenience.
#     x <- input$controller
#     
#     updateDateRangeInput(session, "inDateRange",
#                          label = paste("Date range label", x),
#                          start = paste("2013-01-", x, sep=""),
#                          end = paste("2013-12-", x, sep=""))
#     }),
  
  output$timePlot <- renderPlot({   
    ##Subsets the data.
    start_dat <- data[data$start>=input$date_range[1]&data$start<=input$mid_date,]
    end_dat <- data[data$start>=input$mid_date&data$start<=input$date_range[2],]
    all_dat <- data[data$start>=input$date_range[1]&data$start<=input$date_range[2],]
    
    ##Adds a trend line to the plot.
    model <- lm(cats~start,data=start_dat)
    
    ##Adds the confidence and predictive interval to the dataset.
    prd <- predict(model,newdata=all_dat,interval="prediction",level=0.95,type="response")
    cnf <- predict(model,newdata=all_dat,interval="confidence",level=0.95,type="response")
    
    ##Creates a plot
    plot(cats~start,data=start_dat,type='n',xlim=range(all_dat$start),
         xlab="Time",ylab="Search Freqency",ylim=c(min(prd[,2]),max(prd[,3])))
    
    polygon(c(all_dat$start,rev(all_dat$start)),c(prd[,2],rev(prd[,3])), col="grey80",border=NULL,lty="blank")
    polygon(c(all_dat$start,rev(all_dat$start)),c(cnf[,2],rev(cnf[,3])), col="grey60",border=NULL,lty="blank")
    lines(start_dat$start,start_dat$cats,type='l',col=rgb(0,0,1,0.8),lwd=1.5)
    abline(v=input$mid_date)
    lines(end_dat$start,end_dat$cats,type="l",lty=2,col=rgb(0,0,1,0.8),pch=20)
    lines(all_dat$start,prd[,1],lty=2)
    
  })
  
})
