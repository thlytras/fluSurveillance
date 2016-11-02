source("include.R")

sentinelYears <- rev((resAll$yearweek[(resAll$yearweek %% 100)==40]-40)/100)

library(RColorBrewer)
library(shiny)

ui <- shinyUI(fluidPage(
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css"),
    tags$link(rel = "stylesheet", type = "text/css", href = "mine.css")
  ),
  
  titlePanel("Επιδημιολογική επιτήρηση της γρίπης"),
  tabsetPanel(
    tabPanel("Διάγραμμα sentinel έτους",
      sidebarLayout(
        sidebarPanel(
          selectInput("selSentYr", "Διαθέσιμα έτη", sentinelYears, selectize=TRUE,
                      selected=sentinelYears[1], multiple=TRUE),
          checkboxInput("sentCI", "Διάστημα εμπιστοσύνης"),
          img(src='keelpno.png', width=199, height=157, 
              style="display: block; margin-left: auto; margin-right: auto;")
        ), mainPanel(
          uiOutput("sentinelPlotPanel")
        )
      )
    ),
    tabPanel("Διαχρονικό διάγραμμα sentinel",
      sidebarLayout(
        sidebarPanel(
          dateRangeInput("diaxDateRange", "Εύρος ημερομηνιών", 
                         start=isoweekStart(diaxyear[1]), end=isoweekStart(diaxyear[2])+6, 
                         min=isoweekStart(diaxyear[1]), max=isoweekStart(diaxyear[2])+6, 
                         format="dd-mm-yyyy", lang="el", separator=" έως "),
          checkboxInput("sentCIdiax", "Διάστημα εμπιστοσύνης"),
          img(src='keelpno.png', width=199, height=157, 
              style="display: block; margin-left: auto; margin-right: auto;")
        ), mainPanel(
          plotOutput("diaxPlot")
        )
     )
    ),
    tabPanel("Sentinel κατά περιφέρεια",
      sidebarLayout(
        sidebarPanel(
          selectInput("selSentYrNuts", "Έτος", sentinelYears, selectize=TRUE,
                      selected=sentinelYears[1]),
          checkboxInput("sentCInuts", "Διάστημα εμπιστοσύνης"),
          img(src='keelpno.png', width=199, height=157, 
              style="display: block; margin-left: auto; margin-right: auto;")
        ), mainPanel(
          plotOutput("sentinelNutsPlot")
        )
      )
    ),
    tabPanel("Sentinel κατά αστικότητα",
      sidebarLayout(
        sidebarPanel(
          selectInput("selSentYrAsty", "Έτος", sentinelYears, selectize=TRUE,
                      selected=sentinelYears[1]),
          checkboxInput("sentCIasty", "Διάστημα εμπιστοσύνης"),
          img(src='keelpno.png', width=199, height=157, 
              style="display: block; margin-left: auto; margin-right: auto;")
        ), mainPanel(
          plotOutput("sentinelAstyPlot")
        )
      )
    )
  )
))


server <- shinyServer(function(input, output) {

  output$sentinelPlotPanel <- renderUI({
      plotOutput("sentinelPlot")
  })
  
  output$sentinelPlot <- renderPlot({
    if (is.null(input$selSentYr)) {
      selSentYr <- sentinelYears[1]
    } else {
      selSentYr <- as.integer(input$selSentYr)
    }
    colPal <- c(brewer.pal(9, "Set1")[-6], rainbow(8))
    if (sum(selSentYr>=2014)>0 & sum(selSentYr>=2014)<length(selSentYr)) {
      sentinel_graph(selSentYr, ci=input$sentCI, col=colPal, yaxis2=selSentYr[selSentYr<2014], mult=1/5, 
                     ylab="Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις\n(Νέο σύστημα επιτήρησης, από 2014-2015",
                     ylab2="Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις\n(Παλιό σύστημα επιτήρησης, έως 2013-2014")
    } else {
      sentinel_graph(selSentYr, ci=input$sentCI, col=colPal)
    }
  })

  output$diaxPlot <- renderPlot({
    dxy <- diaxyear
    if (!is.null(input$diaxDateRange)) {
      dxy <- isoweek(as.Date(input$diaxDateRange), "both_num")
    }
    diax_graph(dxy, ci=input$sentCIdiax, alpha=0.35)
  })
  
  output$sentinelNutsPlot <- renderPlot({
    if (is.null(input$selSentYrNuts)) {
      selSentYrNuts <- sentinelYears[1]
    } else {
      selSentYrNuts <- as.integer(input$selSentYrNuts)
    }
    sentinelGraphByGroup(resNutsAll, selSentYrNuts, ci=input$sentCInuts, alpha=0.1)
  })

  output$sentinelAstyPlot <- renderPlot({
    if (is.null(input$selSentYrAsty)) {
      selSentYrAsty <- sentinelYears[1]
    } else {
      selSentYrAsty <- as.integer(input$selSentYrAsty)
    }
    sentinelGraphByGroup(resAstyAll, selSentYrAsty, ci=input$sentCIasty, alpha=0.1)
  })
  
})

shinyApp(ui = ui, server = server)

