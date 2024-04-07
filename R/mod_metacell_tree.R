#' @title Displays a correlation matrix of the quantitative data of a
#' numeric matrix.
#'
#' @description
#' xxxx
#' 
#' @param id xxx
#' @param type xxx
#' @param reset xxx
#'
#' @name metacell-tree
#' 
#' @return NA
#'
#' @examplesIf interactive()
#' data(Exp1_R25_pept, package = 'DaparToolshedData')
#' shiny::runApp(mod_metacell_tree(Exp1_R25_pept[[1]]))
#' 
#' data(Exp1_R25_prot, package = 'DaparToolshedData')
#' shiny::runApp(mod_metacell_tree(Exp1_R25_prot[[1]]))
#'
NULL



#' @import shinyBS
#' @import highcharter
#' @import shinyjs
#' 
#' @rdname metacell-tree
#' 
#' @export
#' 
mod_metacell_tree_ui <- function(id) {
    ns <- NS(id)
   require(shinyBS)
    fluidPage(
        shinyjs::useShinyjs(),
        fluidRow(
            column(width=6, 
                   actionButton(ns("openModalBtn"),
                   tagList(
                       #p('Select tags'),
                       tags$img(src = "images/metacelltree.png", height = "50px"))
                   ),
              shinyBS::bsTooltip(ns("openModalBtn"), "Click to open tags selection tool",
                             "right", options = list(container = "body"))),
            column(width=6, uiOutput(ns('selectedNodes'))
            )
            ),
        br(),
        uiOutput(ns('modaltree'))
        )

}



#' @import shinyBS
#' @import highcharter
#' @import shinyjs
#' 
#' @rdname metacell-tree
#' 
#' @export
#' 
mod_metacell_tree_server <- function(id, 
                                     obj = reactive({NULL}),
                                     reset = reactive({NULL})) {
   
  #pkgs.require(c("shinyBS", "shinyjs"))
    
    
    convertWidgetName <- function(name){
        # This function implements the transformations used to
        # create the names of the checkboxes
        ll <- lapply(name, function(x){
            tmp <- gsub('.', '', x, fixed = TRUE)
            tmp <- gsub(' ', '', tmp, fixed = TRUE)
            tmp <- tolower(tmp)
            tmp <- paste0(tmp, '_cb')
        })
        
        return(unlist(ll))
    }
    
    
    BuildMapping <- function(meta){
        mapping <- c()
        ll <- meta$node
        .ind <- which(ll == 'Any')
        ll <- ll[-.ind]
        colors <- setNames(meta$color[-.ind], nm = convertWidgetName(ll))
        widgets.names <- setNames(ll, nm = convertWidgetName(ll))
        return(list(names = widgets.names, colors = colors))
    }
    
    reverse.mapping <- function(mapping, target){
        req(mapping)
        return(names(mapping)[which(mapping == target)])
    }


    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        dataOut <- reactiveValues(
            trigger = NULL,
            values = ''
        )
        
        rv <- reactiveValues(
            tags = NULL,
            mapping = NULL,
            bg_colors = NULL,
            autoChanged = FALSE
        )
        
        
        Get_bg_color <- function(name){
            gsub("#bg-color#", rv$bg_colors[name], .style)
        }
        
        output$modaltree <- renderUI({
            tagList(
                shinyjs::inlineCSS(css),
                tags$script(paste0('$( document ).ready(function() {
                $("#', 
                ns('modalExample'), 
                '").on("hidden.bs.modal", function (event) {
                x = new Date().toLocaleString();
                Shiny.onInputChange("', ns('lastModalClose'), '",x);});})')),
            tags$head(tags$style(paste0(".modal-dialog { width: fit-content !important; z-index: 1000;}"))),
            #tags$head(tags$style("#modalExample{ display:none;")),
            
            shinyBS::bsModal(ns("modalExample"),
                 title = '',
                 # tagList(
                 #     p('Cells metadata tags'),
                 #     p('To get help about the organization of those tags, please refer to the FAQ')
                 # ),
                 trigger = ns("openModalBtn"),
                 size = "large",
                 #popover_for_help_ui(ns("metacellTag_help")),
                 div(
                     div(style = "align: center;display:inline-block; vertical-align: middle; margin: 5px; padding-right: 0px",
                         p('To get help about the organization of those tags, please refer to the FAQ'),
                         radioButtons(ns('checkbox_mode'), '',
                                      choices = c('Single selection' = 'single',
                                                  'Complete subtree' = 'subtree',
                                                  'Multiple selection' = 'multiple'),
                                      width = '150px')),
                     div(style = "align: center;display:inline-block; vertical-align: middle; margin: 5px; padding-right: 0px",
                         actionButton(ns('cleartree'), 'Clear')
                     )
                 ),
                 uiOutput(ns('tree'))
                 )
            )

            })

        
        output$selectedNodes <- renderUI({
            req(length(dataOut$values) > 0)
            p(paste0('Selected tags: ', paste0(dataOut$values, collapse=',')))
        })
     

        
        init_tree <- function(){
          req(omXplore::get_type(obj()))
          #print('------------ init_tree() ---------------')
          rv$meta <- DAPAR::metacell.def(omXplore::get_type(obj()))
          rv$mapping <- BuildMapping(rv$meta)$names
          rv$bg_colors <- BuildMapping(rv$meta)$colors
          
          tmp <- unname(rv$mapping[names(rv$mapping)])
          rv$tags <- setNames(rep(FALSE, length(tmp)), nm = gsub('_cb', '', tmp))
          rv$autoChanged <- TRUE
        }
        
        
        observeEvent(req(reset()), ignoreInit = TRUE, {
            #print('------------ observeEvent(req(reset()) ---------------')
            # init_tree()
            # update_CB()
            # updateRadioButtons(session, 'checkbox_mode', selected = 'single')
            # rv$autoChanged <- TRUE
            # dataOut$trigger <- as.numeric(Sys.time())
            # dataOut$values <- NULL
            
            if (!is.null(omXplore::get_type(obj())))
                init_tree()
            dataOut$trigger <- as.numeric(Sys.time())
            dataOut$values <- NULL
            }) 
        
        observeEvent(input$openModalBtn,{
            
            init_tree()
            update_CB()
            updateRadioButtons(session, 'checkbox_mode', selected = 'single')
            rv$autoChanged <- FALSE
            #dataOut$trigger <- as.numeric(Sys.time())
            #dataOut$values <- NULL
        })  
                     

# When OK button is pressed, attempt to load the data set. If successful,
# remove the modal. If not show another modal, but this time with a failure
# message.
observeEvent(input$lastModalClose,  ignoreInit = FALSE, ignoreNULL = TRUE, {
    #print('------------ input$lastModalClose ---------------')
    dataOut$trigger <- as.numeric(Sys.time())
    dataOut$values <- names(rv$tags)[which(rv$tags == TRUE)]
    #browser()
})




observeEvent(id, ignoreInit = FALSE, {
  if (!is.null(omXplore::get_type(obj())))
    init_tree()
  dataOut$trigger <- as.numeric(Sys.time())
  dataOut$values <- NULL
}, priority = 1000)



observeEvent(req(input$cleartree), ignoreInit = TRUE, {
    
  update_CB()
  updateRadioButtons(session, 'checkbox_mode', selected = 'single')
    
  rv$autoChanged <- TRUE
})


observeEvent(input$checkbox_mode, {
     
    update_CB()
    ind <- which(rv$meta$parent == 'Any')
    ll <- rv$meta$node[ind]
    ll.widgets <- switch(input$checkbox_mode,
                         single = {
                             for (l in names(rv$mapping))
                             shinyjs::toggleState(l, TRUE)
                         },
                         subtree = {
                             ll.widgets <- names(rv$mapping)[-match(ll, rv$mapping[names(rv$mapping)])]
                             for (l in ll.widgets)
                                 shinyjs::toggleState(l, FALSE)
                             },
                         multiple = {
                             for (l in names(rv$mapping))
                                 shinyjs::toggleState(l, TRUE)
                         }
                         )
    
    
    rv$autoChanged <- FALSE
})


output$tree <- renderUI({
    div(style = "overflow-y: auto;",
        uiOutput(ns(paste0('metacell_tree_', omXplore::get_type(obj()))))
    )
})

.style <- 'vertical-align: top; background: #bg-color#; color: white; padding: 5px;'


# Define tree for protein dataset
output$metacell_tree_protein <- renderUI({
    
  nb <- GetNbTags(obj())
  pourcentages <- round(100*nb/sum(nb), digits = 1)
  
  div(class='wtree',
    tags$ul(
      tags$li(
        checkboxInput(ns('quantified_cb'),
          tags$span(style = Get_bg_color('quantified_cb'), 
            paste0('Quantified (', nb['Quantified'], ' - ', pourcentages['Quantified'], '%)'))),
        tags$ul(
          tags$li(
            checkboxInput(ns('quantbydirectid_cb'),
              tags$span(style = Get_bg_color('quantbydirectid_cb'), 
                paste0('Quant. by direct id (', nb['Quant. by direct id'], ' - ', pourcentages['Quant. by direct id'], '%)')))),
          tags$li(
            checkboxInput(ns('quantbyrecovery_cb'),
              tags$span(style = Get_bg_color('quantbyrecovery_cb'), 
                paste0('Quant. by recovery (', nb['Quant. by recovery'], ' - ', pourcentages['Quant. by recovery'], '%)'))))
        )),
      
      
      tags$li(
        checkboxInput(ns('missing_cb'),
          tags$span(style = Get_bg_color('missing_cb'), 
            paste0('Missing (', nb['Missing'], ' - ', pourcentages['Missing'], '%)'))),
        tags$ul(
          tags$li(
            checkboxInput(ns('missingpov_cb'),
              tags$span(style = Get_bg_color('missingpov_cb'), 
                paste0('Missing POV (', nb['Missing POV'], ' - ', pourcentages['Missing POV'], '%)')))),
          tags$li(
            checkboxInput(ns('missingmec_cb'),
              tags$span(style = Get_bg_color('missingmec_cb'), 
                paste0('Missing MEC (', nb['Missing MEC'], ' - ', pourcentages['Missing MEC'], '%)'))))
        )),
      
      
      tags$li(
        checkboxInput(ns('imputed_cb'),
          tags$span(style = Get_bg_color('imputed_cb'), 
            paste0('Imputed (', nb['Imputed'], ' - ', pourcentages['Imputed'], '%)'))),
        tags$ul(
          tags$li(
            checkboxInput(ns('imputedpov_cb'),
              tags$span(style = Get_bg_color('imputedpov_cb'), 
                paste0('Imputed POV (', nb['Imputed POV'], ' - ', pourcentages['Imputed POV'], '%)')))),
          tags$li(
            checkboxInput(ns('imputedmec_cb'),
              tags$span(style = Get_bg_color('imputedmec_cb'), 
                paste0('Imputed MEC (', nb['Imputed MEC'], ' - ', pourcentages['Imputed MEC'], '%)'))))
        )
      ),
      
      tags$li(
        checkboxInput(ns('combinedtags_cb'),
          tags$span(style = Get_bg_color('combinedtags_cb'), 
            paste0('Combined tags (', nb['Combined tags'], ' - ', pourcentages['Combined tags'], '%)')))
        # tags$ul(
        #     tags$li(
        #         checkboxInput(ns('partiallyquantified_cb'),
        #                       tags$span(style = .style, 'Partially quantified')
        #         )
        #     )
        # )
      )
      
    )
  )
})




output$metacell_tree_peptide <- renderUI({
  nb <- GetNbTags(obj())
  pourcentages <- round(100*nb/sum(nb), digits = 1)
  
  div(class='wtree',
    tags$ul(
      tags$li(
        checkboxInput(ns('quantified_cb'),
          tags$span(style = Get_bg_color('quantified_cb'),
            paste0('Quantified (', nb['Quantified'], ' - ', pourcentages['Quantified'], '%)'))
        ),
        tags$ul(
          tags$li(
            checkboxInput(ns('quantbydirectid_cb'),
              tags$span(style = Get_bg_color('quantbydirectid_cb'),
                paste0('Quant. by direct id (', nb['Quant. by direct id'], ' - ', pourcentages['Quant. by direct id'], '%)'))
            )
          ),
          tags$li(
            checkboxInput(ns('quantbyrecovery_cb'),
              tags$span(style = Get_bg_color('quantbyrecovery_cb'),
                paste0('Quant. by recovery (', nb['Quant. by recovery'], ' - ', pourcentages['Quant. by recovery'], '%)'))
            )
          )
        )
      ),
      
      
      tags$li(
        checkboxInput(ns('missing_cb'),
          tags$span(style = Get_bg_color('missing_cb'),
            paste0('Missing (', nb['Missing'], ' - ', pourcentages['Missing'], '%)'))
        ),
        tags$ul(
          tags$li(
            checkboxInput(ns('missingpov_cb'),
              tags$span(style = Get_bg_color('missingpov_cb'),
                paste0('Missing POV (', nb['Missing POV'], ' - ', pourcentages['Missing POV'], '%)'))
            )
          ),
          tags$li(
            checkboxInput(ns('missingmec_cb'),
              tags$span(style = Get_bg_color('missingmec_cb'),
                paste0('Missing MEC (', nb['Missing MEC'], ' - ', pourcentages['Missing MEC'], '%)'))
            )
          )
        )
      ),
      
      
      tags$li(
        checkboxInput(ns('imputed_cb'),
          tags$span(style = Get_bg_color('imputed_cb'),
            paste0('Imputed (', nb['Imputed'], ' - ', pourcentages['Imputed'], '%)'))
        ),
        tags$ul(
          tags$li(
            checkboxInput(ns('imputedpov_cb'),
              tags$span(style = Get_bg_color('imputedpov_cb'),
                paste0('Imputed POV (', nb['Imputed POV'], ' - ', pourcentages['Imputed POV'], '%)'))
            )
          ),
          tags$li(
            checkboxInput(ns('imputedmec_cb'),
              tags$span(style = Get_bg_color('imputedmec_cb'),
                paste0('Imputed MEC (', nb['Imputed MEC'], ' - ', pourcentages['Imputed MEC'], '%)'))
            )
          )
        )
      )
      
    )
  )
})







update_CB <- function(nametokeep=NULL){
    widgets_to_disable <- names(rv$mapping)
    
    if(!is.null(nametokeep))
        widgets_to_disable <- names(rv$mapping)[-match(reverse.mapping(rv$mapping, nametokeep), names(rv$mapping))]
    
    lapply(widgets_to_disable, function(x){
        updateCheckboxInput(session, x, value = FALSE)
        rv$tags[rv$mapping[x]] <- FALSE
        }
    )
    
    rv$autoChanged <- TRUE
    
}



somethingChanged <- reactive({
  events <- unlist(lapply(names(rv$mapping), function(x) input[[x]]))
    compare <- rv$tags == events
    length(which(compare==FALSE)) > 0
})



# Catch a change in the selection of a node
observeEvent(somethingChanged(), ignoreInit = TRUE, {
    req(length(names(rv$mapping)) > 0)
    if (rv$autoChanged){
        rv$autoChanged <- FALSE
        return (NULL)
    }
    
    # Get the values of widgets corresponding to nodes in the tree
    events <- unlist(lapply(names(rv$mapping), function(x) input[[x]]))
    
    compare <- rv$tags == events
    
    # Deduce the new selected node
    newSelection <- names(rv$tags)[which(compare==FALSE)]
    #print(paste0('newSelection = ', paste0(newSelection, collapse = ', ')))
    # Update rv$tags vector with this new selection
    if (length(newSelection) > 0) {
        for (i in newSelection)
            rv$tags[i] <- input[[reverse.mapping(rv$mapping, i)]]
        
        
        switch(input$checkbox_mode,
               single = {
                   update_CB(newSelection)
                   },
               subtree = {
                   level <- omXplore::get_type(obj())
                   # As the leaves are disabled, this selection is a node
                   # by default, all its children must be also selected
                   for (i in newSelection){
                       if (i %in% metacell.def(level)$parent) {
                           #browser()
                           childrens <- Children(level, i)
                           if (!is.null(childrens) && length(childrens)>0){
                               lapply(childrens, function(x){
                                   updateCheckboxInput(session, 
                                               reverse.mapping(rv$mapping, x), 
                                               value = input[[reverse.mapping(rv$mapping, i)]])
                                   rv$tags[x] <- input[[reverse.mapping(rv$mapping, i)]]
                                   })
                               rv$autoChanged <- TRUE
                           }
                       }
                   }
               },
               multiple = {}
        )
        
    }
    
})




reactive({dataOut})

}
)


}




# This css is adapted from: https://codepen.io/willpower/pen/pJKdej
css <- "
* {
    margin: 0;
    padding: 0;
    box-sizing: border-box;
}

body {
    padding: 00px;
    font-family: helvetica, arial, sans-serif;
}

ul {
    margin: -10px 0px 30px 20px;
}

.wtree li {
    list-style-type: none;
    margin: 0px 0 -10px 10px;
    position: relative;
}
.wtree li:before {
    content: '';
    position: absolute;
    top: -15px;
    left: -25px;
    border-left: 1px solid #ddd;
    border-bottom: 1px solid #ddd;
    width: 20px;
    height: 30px;
}
.wtree li:after {
    position: absolute;
    content: '';
    top: 15px;
    left: -25px;
    border-left: 1px solid #ddd;
    border-top: 1px solid #ddd;
    width: 20px;
    height: 100%;
}
.wtree li:last-child:after {
    display: none;
    padding: 0px 0px 10px 0px;
}
.wtree li span {
    display: inline-block;
    border: 0px solid #ddd;
    border-radius: 10px;
    text-align: center;
    vertical-align: middle;
    padding: 0px 5px 0px 0px;
    color: #888;
    text-decoration: none;
    width: 150px;
}"




#' @export
#' @rdname metacell-tree
#' 
mod_metacell_tree <- function(obj){
  ui <- fluidPage(
    mod_metacell_tree_ui('tree')
    )
  
  server <- function(input, output) {
    
    tags <- mod_metacell_tree_server('tree', 
      obj = reactive({obj}))
    
    observeEvent(req(tags()$trigger), {
      print(tags()$values)
    })
  }
  
  app <- shinyApp(ui, server)
  
}