library(shiny)

    ui <- mod_inputGroup_ui("help")

    server <- function(input, output, session) {
        mod_inputGroup_server("help",
            title = "Foo",
            content = "xxx"
        )
    }

    shinyApp(ui = ui, server = server)
