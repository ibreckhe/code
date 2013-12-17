library(shiny)
library(ngramr)
library(xtable)

# Define server logic required to generate and trend data from google trends
shinyServer(function(input, output) {
  ##Gets data from Google
  data <- reactive({ngram(as.character(input$term),corpus="eng_2012",year_end=2012,smoothing=0)})
  
  ##Takes inputs from the date range input.
  year_start <- reactive({input$start})
  year_end <- reactive({input$end})
  
  ##Takes input from the subsampling control.
  sub_factor <- reactive({input$sub})
  
  ##Takes input from the lagging control
  lag_factor <- reactive({input$lag})
  
  output$ols_timePlot <- renderPlot({
    ##Subsets the data.
    data2 <- data()[data()$Year>=year_start() & data()$Year<=year_end(),]
    
    ##Plots the data.
    plot(data2$Year,data2$Frequency,type='l',col="slateblue",xlab="Year",ylab="Ngram Frequency",main="Trend over Time")
    points(data2$Year,data2$Frequency,pch=21,bg="white")
    
    ##Adds a trend line to the plot.
    model <- lm(data2$Frequency~data2$Year)
    if(input$trend==T){
      abline(model,lty=2)
     }
  })
  
  output$ols_diagPlot <- renderPlot({
    if(input$diag==T){
      data2 <- data()[data()$Year>=year_start() & data()$Year<=year_end(),]
      model <- lm(data2$Frequency~data2$Year)
      resid <- residuals(model,type="response")
      acf <- acf(resid)
      par(mfrow=c(1,2))
      plot(data2$Year,resid,type="p",main="Residuals",ylab="Residuals",xlab="Year")
      abline(h=0)
      plot(acf,main="Auto-correlation Function",xlab="Lag (Years)")
    }
    
  output$ols_modTable <- renderTable({
    if(input$mod==T){
      data2 <- data()[data()$Year>=year_start() & data()$Year<=year_end(),]
      model <- lm(data2$Frequency~data2$Year)
      coefs <- summary(model)$coefficients
      rownames(coefs) <- c("Int.","Year")
      return(coefs)
    }
  },digits=-2)
    
    output$sub_timePlot <- renderPlot({
      ##Subsets the data.
      data2 <- data()[data()$Year>=year_start() & data()$Year<=year_end(),]
      
      ##Plots the data.
      plot(data2$Year,data2$Frequency,type='l',col="slateblue",xlab="Year",ylab="Ngram Frequency",main="Trend over Time")
      
      ##Subsets the data.
      data3 <- data2[c(rep(FALSE,sub_factor()),TRUE),]
      points(data3$Year,data3$Frequency,pch=21,bg="white")
      
      ##Adds a trend line to the plot.
      model <- lm(data3$Frequency~data3$Year)
      if(input$trend==T){
        abline(model,lty=2)
      }
    })
    
    output$sub_diagPlot <- renderPlot({
      if(input$diag==T){
        data2 <- data()[data()$Year>=year_start() & data()$Year<=year_end(),]
        ##Subsets the data.
        data3 <- data2[c(rep(FALSE,sub_factor()),TRUE),]
        
        model <- lm(data3$Frequency~data3$Year)
        resid <- residuals(model,type="response")
        acf <- acf(resid)
        par(mfrow=c(1,2))
        plot(data3$Year,resid,type="p",main="Residuals",ylab="Residuals",xlab="Year")
        abline(h=0)
        plot(acf,main="Auto-correlation Function",xlab="Lag (Years)")
      }
    })
      
      output$sub_modTable <- renderTable({
        if(input$mod==T){
          data2 <- data()[data()$Year>=year_start() & data()$Year<=year_end(),]
          data3 <- data2[c(rep(FALSE,sub_factor()),TRUE),]
          model <- lm(data3$Frequency~data3$Year)
          coefs <- summary(model)$coefficients
          rownames(coefs) <- c("Int.","Year")
          return(coefs)
        }
      },digits=-2)
    
    output$lag_timePlot <- renderPlot({
      ##Subsets the data.
      data2 <- data()[data()$Year>=year_start() & data()$Year<=year_end(),]
      
      ##Plots the data.
      plot(data2$Year,data2$Frequency,type='l',col="slateblue",xlab="Year",ylab="Ngram Frequency",main="Trend over Time")
      points(data2$Year,data2$Frequency,pch=21,bg="white")
      
      ##Creates a lag variable and binds it to the data.
      Frequency_ts <- as.ts(data2$Frequency,frequency=(1/365))
      Frequency_lag <- lag(Frequency_ts,k=lag_factor())
      lag_data <- cbind(Frequency_ts,Frequency_lag)
      data2$Frequency_lag <- lag_data[,2][(lag_factor()+1):length(lag_data[,2])]
      
      ##Adds to the plot
      lines(data2$Year,data2$Frequency_lag,lty=2,col='slateblue')
      
      ##Creates a lagged regression model.
      model <- lm(data2$Frequency~data2$Year+data2$Frequency_lag)
      
      if(input$trend==T){
        data3 <- data.frame(data2$Frequency,data2$Year,Frequency_lag=mean(data2$Frequency_lag,na.rm=T))
        lag_pred <- predict(model,newdata=data3,type='response')
        lines(data2$Year[1:(length(data2$Year)-lag_factor())],lag_pred,lty=2)
      }
    })
    
    output$lag_diagPlot <- renderPlot({
      if(input$diag==T){
        data2 <- data()[data()$Year>=year_start() & data()$Year<=year_end(),]
        
        ##Creates a lag variable and binds it to the data.
        Frequency_ts <- as.ts(data2$Frequency,frequency=(1/365))
        Frequency_lag <- lag(Frequency_ts,k=lag_factor())
        lag_data <- cbind(Frequency_ts,Frequency_lag)
        data2$Frequency_lag <- lag_data[,2][(lag_factor()+1):length(lag_data[,2])]
        
        ##Creates a lagged regression model.
        model <- lm(data2$Frequency~data2$Year+data2$Frequency_lag)
        resid <- residuals(model,type="response")
        acf <- acf(resid)
        par(mfrow=c(1,2))
        plot(data2$Year[1:(length(data2$Year)-lag_factor())],resid,type="p",main="Residuals",ylab="Residuals",xlab="Year")
        abline(h=0)
        plot(acf,main="Auto-correlation Function",xlab="Lag (Years)")
      }
    })
    
    output$lag_modTable <- renderTable({
      if(input$mod==T){
        data2 <- data()[data()$Year>=year_start() & data()$Year<=year_end(),]
        
        ##Creates a lag variable and binds it to the data.
        Frequency_ts <- as.ts(data2$Frequency,frequency=(1/365))
        Frequency_lag <- lag(Frequency_ts,k=lag_factor())
        lag_data <- cbind(Frequency_ts,Frequency_lag)
        data2$Frequency_lag <- lag_data[,2][(lag_factor()+1):length(lag_data[,2])]
        
        ##Creates a lagged regression model.
        model <- lm(data2$Frequency~data2$Year+data2$Frequency_lag)
        coefs <- summary(model)$coefficients
        rownames(coefs) <- c("Int.","Year","Lagged Response")
        return(coefs)
      }
    },digits=-2)
    
  })
})