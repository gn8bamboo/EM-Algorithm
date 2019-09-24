library(shiny)


DealData <- function(data, nModes)
{
  
  premu<-seq(from=min(data),to=max(data),by=(max(data)-min(data))/(nModes-1))
  prepi<-rep(1/nModes,times=nModes)
  presigma<-rep(1,times=nModes)
  data.m <- matrix(NA, nrow = nModes, ncol = 3)
  data.m[,1] <- prepi
  data.m[,2] <- premu
  data.m[,3] <- presigma
  
  list("prepi" = data.m[,1],
       "premu" = data.m[,2],
       "presigma" = data.m[,3])
}

E_Step <- function(data,nModes,pi,mu,sigma)
{
  post.pro <- matrix(0, nrow = length(data), ncol = nModes)
  post.pro.m<-matrix(0, nrow = length(data), ncol = nModes)
  sum.of.prob <- 0
  for (i in 1: nModes) {
    post.pro[,i] <- pi[i]*dnorm(data,mu[i],sigma[i])
    sum.of.prob <- sum.of.prob + post.pro[,i]
  } 
  
  for (i in 1: nModes) {
    post.pro.m[,i] <- post.pro[,i]/sum.of.prob
  } 
  
  sum.of.comps.ln <- log(sum.of.prob, base = exp(1))
  sum.of.comps.ln.sum <- sum(sum.of.comps.ln)
  loglik <- sum.of.comps.ln.sum
  list("loglik" = loglik,
       "posterior" = post.pro.m)
}

M_Step <- function(data,nModes,posterior)
{
  ls <- matrix(0, nrow = nModes, ncol = 3)
  colnames(ls) <- paste(c('pi','mu','sigma'))
  mu <- 0
  pi<- 0
  sigma <- 0
  temp <- 0
  
  for (i in 1: nModes) {
    
    temp[i] <- sum(posterior[, i])
    pi[i] <- temp[i] / length(data)
    mu[i] <- 1/temp[i] * sum(posterior[, i] * data)
    sigma[i] <- sqrt(sum(posterior[, i] * (data - mu[i])^2)/temp[i])
    
    ls[i,] <- c(pi[i],mu[i],sigma[i])
  } 
  
  list("result" = ls)
}


Iteration <- function(x,nModes)
{
  
  data.ready <- DealData(x,nModes)
  pi <- data.ready$prepi
  mu <- data.ready$premu
  sigma <- data.ready$presigma
  
  
  for (i in 1:5000) {
    if (i == 1) {
      # Initialization
      e.step <- E_Step(x,nModes,pi,mu,sigma)
      m.step <- M_Step(x,nModes,e.step[["posterior"]])
      cur.loglik <- e.step[["loglik"]]
      loglik.vector <- e.step[["loglik"]]
    } else {
      # Repeat E and M steps till convergence
      e.step <- E_Step(x,nModes,m.step$result[,1],m.step$result[,2],m.step$result[,3])
      m.step <- M_Step(x,nModes,e.step[["posterior"]])
      loglik.vector <- c(loglik.vector, e.step[["loglik"]])
      loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
      if(loglik.diff < 1e-6) {
        break
      } else {
        cur.loglik <- e.step[["loglik"]]
      }
    }
  }
  list("result" = m.step,
       "logliklihood" = loglik.vector)
}

Stepbystep <- function(x,nModes,times)
{
  
  data.ready <- DealData(x,nModes)
  pi <- data.ready$prepi
  mu <- data.ready$premu
  sigma <- data.ready$presigma
  
  
  for (i in 1:times) {
    if (i == 1) {
      # Initialization
      e.step <- E_Step(x,nModes,pi,mu,sigma)
      m.step <- M_Step(x,nModes,e.step[["posterior"]])
      cur.loglik <- e.step[["loglik"]]
      loglik.vector <- e.step[["loglik"]]
    } else {
      # Repeat E and M steps till convergence
      e.step <- E_Step(x,nModes,m.step$result[,1],m.step$result[,2],m.step$result[,3])
      m.step <- M_Step(x,nModes,e.step[["posterior"]])
      loglik.vector <- c(loglik.vector, e.step[["loglik"]])
      loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
      if(loglik.diff < 1e-6) {
        break
      } else {
        cur.loglik <- e.step[["loglik"]]
      }
    }
  }
  list("result" = m.step,
       "logliklihood" = loglik.vector)
}


function(input,output){
  
  data <- reactive({
    file1 <- input$file
    if(is.null(file1)){return()} 
    read.csv(file=file1$datapath, sep=",", header = input$header)
  })
  
  

  output$obs <- renderTable({
    if(is.null(data()))
    {return()}
    data()
  })
  
  output$summary <- renderPrint({
    if(is.null(data()))
    {return()}
    summary(data())
  })
  
  output$stepbystep <- renderPrint({
    if(is.null(data()))
    {return()}
    x <- data()[,input$pickcolumn]
    nModes = input$numofmodes
    times = input$times
    
    s = Stepbystep(x,nModes,times)
    print(s$result)
    print(s$logliklihood)
  
  })
  
  output$convergency <- renderPrint({
    if(is.null(data()))
    {return()}
    x <- data()[,input$pickcolumn]
    nModes = input$numofmodes
    
    i = Iteration(x,nModes)
    i$logliklihood
    noi <- length(i$logliklihood)
    
    cat("The number of iterations:",noi, "\n" )
    print(i$logliklihood)
  })

  output$emr <- renderPrint({
    
    if(is.null(data()))
    {return()}
    
    x <- data()[,input$pickcolumn]
    nModes = input$numofmodes
    
    i = Iteration(x,nModes)
    i$result
   
  })
  
  output$modelsel <- renderPrint({
    if(is.null(data()))
    {return()}
   
    x <- data()[,input$pickcolumn]
    nModes = input$numofmodes
    
    LogL <- Iteration(x,nModes)$logliklihood[length(Iteration(x,nModes)$logliklihood)]
    k <- 3*nModes-1
    n <- length(x)
    
    AICc <- -2*LogL + 2*k + 2*k*(k+1)/(n-k-1)
    BIC <- -2*LogL + k*log(n,base = exp(1))
    
    cat("The AICc Value is:", AICc,"\n") 
    cat("The BIC Value is:", BIC)
    
  })
  
  output$aic <- renderPlot({
    if(is.null(data()))
    {return()}
    
    x <- data()[,input$pickcolumn]
    nModes = input$numofmodes

    n <- length(x)
    
    LogL <- c(Iteration(x,2)$logliklihood[length(Iteration(x,2)$logliklihood)],
               Iteration(x,3)$logliklihood[length(Iteration(x,3)$logliklihood)],
               Iteration(x,4)$logliklihood[length(Iteration(x,4)$logliklihood)],
              Iteration(x,5)$logliklihood[length(Iteration(x,5)$logliklihood)]
              )
  
    nummode <- c(2,3,4,5)
    AICc <- -2*LogL + 2*(3*nummode -1) + 2*(3*nummode -1)*((3*nummode -1)+1)/(n-(3*nummode -1)-1)
   
    plot(nummode,AICc)
    lines(nummode,AICc)
   
  })
  
  output$bic <- renderPlot({
    if(is.null(data()))
    {return()}
    
    x <- data()[,input$pickcolumn]
    nModes = input$numofmodes
    
    n <- length(x)
    
    
    LogL <- c(Iteration(x,2)$logliklihood[length(Iteration(x,2)$logliklihood)],
              Iteration(x,3)$logliklihood[length(Iteration(x,3)$logliklihood)],
              Iteration(x,4)$logliklihood[length(Iteration(x,4)$logliklihood)],
              Iteration(x,5)$logliklihood[length(Iteration(x,5)$logliklihood)]
    )
    
    
    nummode <- c(2,3,4,5)
    BIC <- -2*LogL + (3*nummode -1)*log(n,base = exp(1))
  
    plot(nummode,BIC)
    lines(nummode,BIC)
  })
  
  output$visual <- renderPlot({
    if(is.null(data()))
    {return()}
    x <- data()[,input$pickcolumn]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    hist(x, breaks = bins,col = 'darkgray', border = 'white',main = "Histograms of the Data", xlab = "Data",prob=TRUE)
    xfit<-seq(min(x),max(x),length=40)
    lines(density(x), lty=2, lwd=2)
   
  })
  
  
  output$content <- renderUI({
    if(is.null(data()))
      h5("No dataset loaded!")
   
  })
  
}



