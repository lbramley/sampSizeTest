library(shiny)

server <- function(input, output) 
{
  output$coveredSlider <- renderUI({
    maxCovered <- ifelse(input$likeN < input$likeX, input$likeN, input$likeX)
    sliderInput("likeX", "Number covered", min = 1, max = input$likeN, maxCovered, step = 1)
  }
  )
  
  observe(
    {
      r <- input$priorA - 1.0
      n <- input$priorA + input$priorB - 2.0
      #
      # Some special cases (i.e. for a Beta(1,1) prior and for extreme proportions)
      #
      if(n == 0) 
      {
        p <- 0.5
      } else 
      {
        p <- r/n
      }
      if(p < 0.1) p <- 0.1
      if(p > 0.9) p <- 0.9
      #
      # Calculate sample size (normal approximation)
      #  
      e  <- input$surveyPrecision/100.0
      s <- (p*(1 - p))/(e/1.96)^2 - n
      #
      # Apply finite population correction (FPC)
      #
      # NOTE : May want to move the sample size suggestion to a separate function to avoid recalculation
      #        when the required precsion is changed.
      #
      #        Also ... may need to revise the sample size suggestion.
      #
      s <- (s*input$surveyPopulation)/(s + (input$surveyPopulation - 1))
      s <- ceiling(s)
      #
      # Apply minimum sample size rule
      #
      if(s < n) s <- min(c(ceiling(n), input$surveyPopulation))
      output$sampleSize <- renderText(paste("Suggested sample size :", s))
      #
      # Calculate parameters to the likelihood distribution (from likeN and likeX)
      #
      likeA <- input$likeX + 1.0
      likeB <- input$likeN - input$likeX + 1.0
      #
      # Calculate parameters for posterior distribution
      #
      postA <- input$priorA + likeA - 1.0
      postB <- input$priorB + likeB - 1.0
      #
      # Calculate data to be plotted
      #
      x <- seq(0, 1, 0.01)
      priorPoints <- dbeta(x, input$priorA, input$priorB)
      likePoints <- dbeta(x, likeA, likeB)
      postPoints <- dbeta(x, postA, postB)
      xPerCent <- x * 100
      
      ymax <- max(priorPoints, likePoints, postPoints)
      #
      # Plot prior, likelihood, and posterior 
      #
      output$bayesPlot <- renderPlot({
        par(cex = 1.3)
        plot(priorPoints ~ xPerCent, type = "l", col = "#425cf4", ylim = c(0, ymax), 
             yaxt = "n", xlab = "Proportion (% Coverage)", ylab = "", lwd = 2)
        lines(likePoints ~ xPerCent, col = "#31d837", lwd = 2)
        lines(postPoints ~ xPerCent, col = "#ff0000", lwd = 2)
      })
      #
      # Estimate posterior mode
      #
      est <- (postA - 1)/(postA + postB - 2)
      #
      # Calculate estimate and 95% CI by bootstrap aggregation (bagging) method
      #
      lciList <- sapply(1:100, function(x) quantile(rbeta(100, postA, postB), na.rm = T, probs = .025))
      uciList <- sapply(1:100, function(x) quantile(rbeta(100, postA, postB), na.rm = T, probs = .975))
      
      lci <- mean(lciList)
      uci <- mean(uciList)
      #
      # Convert to percentages and round to 1 d.p.
      #
      est <- round(est * 100, 1)
      lci <- round(lci * 100, 1)
      uci <- round(uci * 100, 1)
      #
      # Output string
      #
      estimateText <- paste(est, " (", lci, " - ", uci, ")", sep = "")
      output$outputText <- renderText(paste("Estimate: ", estimateText, sep = ""))
      #
      # z-test
      #    
      r1 <- input$priorA - 1.0
      n1 <- input$priorA + input$priorB - 2.0
      #
      # Check for the non-informative prior (avoid division by zero)
      #
      if(n1 == 0)
      {
        zTest <- "** Not Applicable **"
      } else { 
        r2 <- input$likeX * 1.0
        n2 <- input$likeN * 1.0
        p0 <- (r1 + r2)/(n1 + n2)
        p1 <- r1/n1
        p2 <- r2/n2
        z <- (p1 - p2)/sqrt((p0*(1 - p0))*(1/n1 + 1/n2))
        absZ <- abs(z)
        p <- pnorm(absZ, 0, 1)
        p <- 2.0*(1.0 - p) 
        z <- round(z, 2)
        p <- round(p, 4)
        zTest <- paste("z = ", z, ", ", "p = ", p, sep = "")
      }
      output$zTestText <- renderText(paste("z-test:", zTest))
    })
}

ui <- fluidPage(sidebarLayout(
  sidebarPanel(
    sliderInput("priorA", "Prior alpha", min = 1.0, max = 64.0, value = 5.0, step = .1),
    sliderInput("priorB", "Prior beta", min = 1.0, max = 64.0, value = 5.0, step = .1),
    sliderInput("surveyPrecision", "Precision", min = 1.0, max = 20.0, value = 10.0, step = 1.0),
    sliderInput("surveyPopulation", "Population", min = 50, max = 1000, value = 500, step = 50),
    textOutput("sampleSize"),
    #checkboxInput("surveyData", "Use survey data"),
    sliderInput("likeN", "Achieved sample size", min = 10, max = 500, value = 100, step = 1),
    uiOutput("coveredSlider"),
    textOutput("outputText"),
    textOutput("zTestText")
  ),
  mainPanel(
    plotOutput("bayesPlot", height = 550),
    tags$div(class = "header",
             list(
               tags$div(id = "legend",
                 h3(span("Prior", class = "blueline"),
                 span("Posterior", class = "greenline"),
                 span("Likelihood", class = "redline")
               ))
             ))
  )
), theme = "bayesTheme.css"
)

shinyApp(ui = ui, server = server)