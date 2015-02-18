
server.edaseq <- function(files){
  function(input, output, session){
    require("shiny")
    require("ShortRead")
    require("Rsamtools")
    require("DESeq")
    require("IRanges")
    source("./AllClasses.R"); source("./AllGenerics.R"); source("./functions.R")
    source("./methods-BamFileList.R"); source("./methods-FastqFileList.R")
    source("./methods-Others.R"); source("./methods-SeqExpressionSet.R")
    
    output$readcounts <- renderPlot({
      barplot(files[as.integer(input$counts)],las=2)
    })

    output$readqualities <- renderPlot({
      plotQuality(files[as.integer(input$quality)],lty=1)
    })
    
    output$readquality1 <- renderPlot({
      plotQuality(files[[as.integer(input$quality1)]])
    })
    
    
    #output$readcounts <- renderPrint(({input$sampchecks}))
    
  }
}

