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
                            h4("By locus"),
                            uiOutput(ns("by_locus_filter_params")),
                            br(),
                            h4("By individual"),
                            uiOutput(ns("by_ind_filter_params")),
                          actionButton(inputId = ns('filter_btn'), "Filter")

                          ),
                 fluidRow(
                   column(width = 8,
                          titlePanel("Plots"),
                          tabsetPanel(id = "filt_distributions", type = "tabs",
                            tabPanel("Missing by locus",
                                     plotly::plotlyOutput(ns("dist_miss_loc"))
                                     ),
                            tabPanel("Missing by individual",
                                     plotly::plotlyOutput(ns("dist_miss_ind"))
                            ),
                            tabPanel("Heterozygosity by locus",
                                     plotly::plotlyOutput(ns("dist_ho_loc"))
                            ),
                            tabPanel("Heterozygosity by individual",
                                     plotly::plotlyOutput(ns("dist_ho_ind"))
                            ),
                            tabPanel("MAF",
                                     plotly::plotlyOutput(ns("dist_MAF"))
                            ),
                          ),
                          ),
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
        loc_missing <- cgiarGenomics::get_loc_missing(geno_data)
        ind_missing <- cgiarGenomics::get_ind_missing(geno_data)
        loc_he <- cgiarGenomics::get_loc_het(geno_data)
        ind_he <-cgiarGenomics::get_ind_het(geno_data)
        maf <- cgiarGenomics::get_maf(geno_data)


        baseline <- list(loc_missing = loc_missing,
                         ind_missing = ind_missing,
                         loc_he = loc_he,
                         ind_he = ind_he,
                         maf = maf)

        output$dist_miss_loc <- plotly::renderPlotly({
          fig <- plotly::plot_ly(alpha = 0.6)
          fig <- fig %>% plotly::add_histogram(x = baseline$loc_missing, histnorm = "probability", name="Missing data by locus" )
        })

        output$dist_miss_ind <- plotly::renderPlotly({
          fig <- plotly::plot_ly(alpha = 0.6)
          fig <- fig %>% plotly::add_histogram(x = baseline$ind_missing, histnorm = "probability", name="Missing data by indivudual" )
        })

        output$dist_ho_loc <- plotly::renderPlotly({
          fig <- plotly::plot_ly(alpha = 0.6)
          fig <- fig %>% plotly::add_histogram(x = baseline$loc_he, histnorm = "probability", name="Heterozygosity by locus" )
        })

        output$dist_ho_ind <- plotly::renderPlotly({
          fig <- plotly::plot_ly(alpha = 0.6)
          fig <- fig %>% plotly::add_histogram(x = baseline$ind_he, histnorm = "probability", name="Heterozygosity by indivudual" )
        })

        output$dist_MAF <- plotly::renderPlotly({
          fig <- plotly::plot_ly(alpha = 0.6)
          fig <- fig %>% plotly::add_histogram(x = baseline$maf, histnorm = "probability", name="Minor allele frequency" )
        })

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
            sliderInput(ns("loc_missing"),
                        " Minimum number of genotyped individuals by marker:",
                        round(min(baseline$loc_missing)*adegenet::nInd(geno_data$geno)),
                        round(max(baseline$loc_missing)*adegenet::nInd(geno_data$geno)),
                        step = 1,
                        value = round(min(baseline$miss_loc)*adegenet::nInd(geno_data$geno))),
            sliderInput(ns("loc_he"),
                        " Minimum heterozygosity by locus:",
                        round(min(baseline$loc_he),2),
                        round(max(baseline$loc_he),2),
                        step = 0.01,
                        value = round(min(baseline$miss_loc), 2)),
            sliderInput(ns("maf"),
                        " Minimum minor allele frequency (MAF):",
                        min = 0,
                        max = round(max(baseline$maf),2),
                        step = 0.01,
                        value = 0)
            # ADD multiple select to remove target locus
          ))

        output$by_ind_filter_params <- renderUI({
          tags$span(
            sliderInput(ns("ind_missing"),
                        " Minimum percentage of genotyped loci by individual:",
                        round(min(baseline$ind_missing),2),
                        round(max(baseline$ind_missing),2),
                        step = 0.01,
                        value = round(min(baseline$ind_missing), 2)),
            sliderInput(ns("ind_he"),
                        " Minimum percentage of heterozygosity by individual:",
                        round(min(baseline$ind_he), 2),
                        round(max(baseline$ind_he),2),
                        step = 0.01,
                        value = round(min(baseline$ind_he), 2))
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

  })
}

## To be copied in the UI
# mod_qaGenoApp_ui("qaGenoApp_1")

## To be copied in the server
# mod_qaGenoApp_server("qaGenoApp_1")
