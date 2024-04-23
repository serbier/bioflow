#' qaGenoApp UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_qaGenoAdeApp_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::mainPanel(
      width = 12,
      tabsetPanel(
        id = ns("tabsMain"),
        type = "tabs",
        tabPanel(div(icon("book"), "Information-QA-Geno"),
                 br(),
                 shinydashboard::box(status = "success", width = 12,
                                     solidHeader = TRUE,
                                     column(width = 12, style = "height:580px; overflow-y: scroll;overflow-x: scroll;",
                                            tags$body(
                                              column(width = 6,
                                                     h1(strong(span("QA for genetic markers", style="color:green"))),
                                                     h2(strong("Status:")),
                                                     uiOutput(ns("warningMessage")),
                                                     img(src = "www/qaGeno.png", height = 300, width = 650)),
                                              column(width = 6, shiny::plotOutput(ns("plotDataDependencies")), ),
                                              column(width = 12,h2(strong("Details")))
                                       )
                                      )
                                     )
                 ),
        tabPanel(div(icon("arrow-right-to-bracket"), "Input"),
                 fluidRow(
                   column(width = 8,
                          titlePanel("Summary"),
                          ),
                   column(width = 4,
                          titlePanel("Filtering Parameters"),
                                     tabsetPanel(id = "filt_params",type = "tabs",
                                       tabPanel("By locus",
                                                uiOutput(ns("by_locus_filter_params"))
                                                ),
                                       tabPanel("By individual",
                                                uiOutput(ns("by_ind_filter_params"))
                                                )),
                          actionButton(inputId = ns('filter_btn'), "Filter")

                          ),
                 fluidRow(
                   column(width = 8,
                          titlePanel("Plots"),
                          uiOutput(ns("dist_params"))),
                   column(width = 4,
                          titlePanel("Preview"),
                          verbatimTextOutput(ns('log_filt')))
                 ),
              ),
        ),
        tabPanel(div(icon("arrow-right-from-bracket"), "Output"))
      )
    )
  )
}

#' qaGenoApp Server Functions
#'
#' @noRd
mod_qaGenoAdeApp_server <- function(id, data){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    output$plotDataDependencies <- shiny::renderPlot({ dependencyPlot() })
    ############################################################################ clear the console
    hideAll <- reactiveValues(clearAll = TRUE)
    observeEvent(data(), {
      hideAll$clearAll <- TRUE
    })
    ############################################################################
    # warning message
    output$warningMessage <- renderUI(
      if(is.null(data())){
        HTML( as.character(div(style="color: red; font-size: 20px;", "Please retrieve or load your data using the 'Data' tab.")) )
      }else{ # data is there
        if(!is.null(data()$data$geno)){
          HTML( as.character(div(style="color: green; font-size: 20px;", "Data is complete, please proceed to perform the marker QA specifying your input parameters under the Input tabs.")) )
        }else{HTML( as.character(div(style="color: red; font-size: 20px;", "Please retrieve or load your genotype data using the 'Data' tab. ")) )}
      }
    )
    ###
    init_filt_data <- eventReactive(
      req(data()$data$geno),
      {
        geno_data <- list(geno = data()$data$geno,
                          meta = data()$metadata$geno)

        shinybusy::show_modal_spinner('fading-circle', text = 'Initializating filtering stats...')
        # get min and max values of sliders
        prop_typed_loc <- adegenet::propTyped(geno_data$geno, by = 'loc')
        prop_typed_ind <- adegenet::propTyped(geno_data$geno, by = 'ind')
        maf <- adegenet::minorAllele(geno_data$geno)

        baseline <- list(miss_loc = prop_typed_loc,
                         miss_ind = prop_typed_ind,
                         maf = maf)
        shinybusy::remove_modal_spinner()
        return(baseline)
      }
    )

    observeEvent(
      c(req(data()$data$geno),
      req(init_filt_data())),
      {
        print("good data to filter!")

        baseline <- init_filt_data()
        geno_data <- list(geno = data()$data$geno,
                          meta = data()$metadata$geno)

        output$by_locus_filter_params <- renderUI( tags$span(
            checkboxInput(ns("mono_check"),
                          "Filter monomorphic:"),
            sliderInput(ns("miss_marker"),
                        " Minimum number of genotyped individuals by marker:",
                        round(min(baseline$miss_loc)*adegenet::nInd(geno_data$geno)),
                        round(max(baseline$miss_loc)*adegenet::nInd(geno_data$geno)),
                        step = 1,
                        value = round(min(baseline$miss_loc)*adegenet::nInd(geno_data$geno))),
            sliderInput(ns("maf_marker"),
                        " Minimum minor allele frequency (MAF):",
                        min = 0,
                        max = round(max(baseline$maf),2),
                        step = 0.01,
                        value = 0)
            # ADD multiple select to remove target locus
          ))

        output$by_ind_filter_params <- renderUI({
          tags$span(
            sliderInput(ns("miss_ind"),
                        " Minimum percentage of genotyped loci by individual:",
                        min(baseline$miss_ind),
                        max(baseline$miss_ind),
                        step = 0.01,
                        value = min(baseline$miss_ind))
          )
        })
    })

    filter_data <- eventReactive(input$filter_btn,
      {
        req(init_filt_data())
        baseline <- init_filt_data()

        filtered <- list(geno = data()$data$geno,
             meta = data()$metadata$geno)
        print(input$miss_marker)
        print(input$maf_marker)
        print(input$miss_ind)



        # Remove monomorphic markers if check
        if(input$mono_check){
          filtered <- cgiarGenomics::filter_monomorphic(filtered)

        }

        if(input$miss_marker > min(baseline$miss_loc)){
          threshold <- round(input$miss_marker/adegenet::nInd(filtered$geno),3)
          filtered <- cgiarGenomics::filter_missing_rate_by_marker(filtered,
                                                               threshold = threshold)
        }

        if(input$maf_marker > min(baseline$maf)){
          filtered <- cgiarGenomics::filter_MAF(filtered,
                                           threshold = input$maf_marker)
        }

        if(input$miss_ind > min(baseline$miss_ind)){
          filtered <- cgiarGenomics::filter_missing_rate_by_indv(filtered,
                                                                 threshold = input$maf_marker)
        }
        return(filtered)
      })

    output$dist_params <- renderUI({
      req(filter_data())

      filter_data <- filter_data()
      tabsetPanel(id = "filt_params",type = "tabs",
                  tabPanel("Miss by locus",
                           renderPlot({cgiarGenomics::plot_missingnes_by_marker(filter_data)})
                           ),
                  tabPanel("Miss by MAF",
                           renderPlot({cgiarGenomics::plot_MAF(filter_data)})
                           ),
                  tabPanel("Miss by individual",
                           renderPlot({cgiarGenomics::plot_overall_missingness(filter_data)})
                           )
                  )

    })


  })
}

## To be copied in the UI
# mod_qaGenoApp_ui("qaGenoApp_1")

## To be copied in the server
# mod_qaGenoApp_server("qaGenoApp_1")
