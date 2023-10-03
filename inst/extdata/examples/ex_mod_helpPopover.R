library(shiny)

    ui <- mod_helpPopover_ui("help")

    server <- function(input, output, session) {
        mod_helpPopover_server("help",
            title = "Foo",
            content = "xxx"
        )
    }

    shinyApp(ui = ui, server = server)
