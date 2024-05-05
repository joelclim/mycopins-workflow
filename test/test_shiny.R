library(shiny)
ui <- fluidPage(
  "Hello, world!!!",
  textInput("name", "Enter a name:"),
  textOutput("q")
)

server <- function(input, output, session) {
  output$q <- renderText({
    paste("Do you prefer dogs or cats,", input$name, "?")
  })  
}

shinyApp(ui = ui, server = server)