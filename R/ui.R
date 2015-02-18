require(shiny)

ui.edaseq <- function(files){
  
  ## shinyUI(navbarPage(
  navbarPage("EDASeq",
    tabPanel("Counts",
      sidebarLayout(
        sidebarPanel(
        checkboxGroupInput("counts", label = "Samples", 
                           choices = sapply(names(files),function(f){return(which.max(names(files)==f))}),
                           selected = 1:length(files))
        ),
        mainPanel(plotOutput("readcounts"))
        #mainPanel(verbatimTextOutput("readcounts"))
      )),
    tabPanel("Avg Quality",
      sidebarLayout(
        sidebarPanel(
          checkboxGroupInput("quality", label = "Samples", 
                             choices = sapply(names(files),function(f){return(which.max(names(files)==f))}),
                             selected = 1:length(files))
        ),
        mainPanel(plotOutput("readqualities"))
      )),
    tabPanel("Sample Quality",
             sidebarLayout(
               sidebarPanel(
                 radioButtons("quality1", label = "Sample", 
                                    choices = sapply(names(files),function(f){return(which.max(names(files)==f))}),
                                    selected = 1)
               ),
               mainPanel(plotOutput("readquality1"))
      ))
  )
}

#read counts, lengths, seqn analysis

