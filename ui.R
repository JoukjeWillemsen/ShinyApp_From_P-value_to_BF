ui <- fluidPage(
  titlePanel("From critical P-value to critical BF for binomial data"),
   fluidRow(
     column(4, 
           "Frequentist test",
           numericInput("N", "samplesize", value = 20, min = 1, max = 10000, step = 1),
           numericInput("critp", "critical P", value = 0.05, min = 0.0001, max = 0.1, step = 0.01),
           numericInput("nullvalue", "proportion corresponding to H0", value = 0.5, min = 0, max = 1, step = 0.1)
    ),
    column(4, 
           "Bayesian test",
           numericInput("aprior", "a", value = 1, min = 0, max = 1000),
           numericInput("bprior", "b", value = 1, min = 0, max = 1000),
           
           selectInput("ratio", label = "Ratio", choices = c("BF10", "BF01"), selected = "BF01")
    ),
    column(4, 
           "y-axis BF graph",
           sliderInput("ymax", "max",
                       min = 0.1, max = 100,
                       value = 5)

    )

    ),
  
  fluidRow(
    column(9, plotOutput("hist1")),
    column(3,
           h4("The Bayes Factor corresponding is:"),
           textOutput("summary"))
  )
)

