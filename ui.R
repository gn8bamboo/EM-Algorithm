library(shiny)

fluidPage(
  titlePanel("EM Algorithm in R -- Probability Project"),
  sidebarLayout(
    sidebarPanel(
      fileInput('file', 'Choose CSV File'),
      checkboxInput('header', 'Header', TRUE),
      tags$hr(),
      numericInput("numofmodes", "Number of Modes:", min = 2, max = 10, value = 2),
      numericInput("pickcolumn", "Pick the Column:", min = 1, value = 1),
      numericInput("times", "Number of iterations:",min = 1,max = 1000,value = 1),
      sliderInput("bins", "Number of bins:",min = 1,max = 50,value = 30),
      tags$hr(),
      submitButton("Update View")
    ),
    
    mainPanel(
      uiOutput('content'),
      tabsetPanel(tabPanel("Data",tableOutput("obs")),
                  tabPanel("Summary",verbatimTextOutput("summary")),
                  tabPanel("Convergency",verbatimTextOutput("stepbystep"),verbatimTextOutput("convergency")),
                  tabPanel("EM Result",verbatimTextOutput("emr")),
                  tabPanel("Model Selection",verbatimTextOutput("modelsel"),plotOutput("aic"),plotOutput("bic")),
                  tabPanel("Visualisation",plotOutput("visual")))
    )
  )
)