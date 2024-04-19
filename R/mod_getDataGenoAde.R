geno_example  <- 'www/example/geno.hmp.txt'

# Currently polyploid supported formats
polyploid_support <- c("vcf", 'darttag', 'hapmap')
dartseq_formats <- c("dartseqpa", "dartsesnp")

mod_getDataGenoAde_ui <- function(id) {
  ns <- NS(id)
  tagList(fluidRow(
    column(
      width = 6,
      wellPanel(
        h4("File Input:"),
        selectInput(
          inputId = ns('adegeno_input'),
          label   = 'Genotypic SNPs Source*:',
          choices = list(
            'HapMap Upload' = 'hapmap',
            'VCF Upload' = 'vcf',
            'DartSeq PA' = 'dartseqpa',
            'DartSeq SNP' = 'dartsesnp',
            'DartTag SNP' = 'darttag'
          ),
          width   = '200px'
        ),

        tags$span(
          id = ns('adegeno_file_holder'),
          fileInput(
            inputId = ns('adegeno_file'),
            label   = NULL,
            width   = '400px',
            accept  = c('application/gzip', '.gz', '.txt', '.hmp', '.csv', '.vcf')
          ),
          # counts file input
          uiOutput(ns("darttag_dosage_file")),
        )
      )
    ),
    column(
      width = 6,
      wellPanel(
        h4("Reading Parameters:"),
        selectInput(
          inputId = ns("ploidlvl_input"),
          label = "Ploidity level:",
          choices = list("Diploid" = 2,
                         "Tetraploid" = 4),
        ),
        # Dartseq input params
        uiOutput(ns("dartseq_params")),
        actionButton(inputId = ns("load_geno"), "Load Data")
      )
    )
  ),
  hr(),
  # Outputs
  tableOutput(ns('summary_by_chrom'))
  )
}


mod_getDataGenoAde_server <-
  function(id, data = NULL, res_auth = NULL) {
    print('Is this executing?')
    moduleServer(id, function(input, output, session) {
      ns <- session$ns

      # reactive to format select
      observeEvent(input$adegeno_input,
                   if (length(input$adegeno_input) > 0) {
                     output$dartseq_params <- renderUI({return(NULL)})
                     output$darttag_dosage_file <- renderUI({return(NULL)})


                     if (input$adegeno_input %in% polyploid_support) {
                       print('Polyploid supported format')
                       updateSelectInput(session, "ploidlvl_input",
                                         choices = list("Diploid" = 2, "Tetraploid" = 4))
                     } else {
                       # ADD should limit the obtions of poly lvl to only diploid
                       print('Polyploid unsupported format')
                       updateSelectInput(session, "ploidlvl_input",
                                         choices = list("Diploid" = 2))
                     }
                     if (input$adegeno_input %in% dartseq_formats){
                       print('dartseq format')
                       output$dartseq_params = renderUI({
                         tags$span(
                           selectInput(
                             inputId = ns("dartseq_markerid"),
                             label = "Marker id:",
                             choices = list()
                           ),
                           selectInput(
                             inputId = ns("dartseq_chrom"),
                             label = "Chromosome",
                             choices = list()
                           ),
                           selectInput(
                             inputId = ns("dartseq_position"),
                             label = "Variant position",
                             choices = list()
                           )
                        )
                       })
                     }
                     if (input$adegeno_input == "darttag"){
                       print('dartag format')
                       output$darttag_dosage_file <- renderUI({
                         fileInput(
                           inputId = ns('darttag_dosage_file'),
                           label   = "Dosage file:",
                           width   = '400px',
                           accept  = c('.csv')
                         )
                       })
                     }
                   })

      # Reactive to select format and after upload file
      observeEvent(c(input$adegeno_input, input$adegeno_file),
                   {
                     req(input$adegeno_input, input$adegeno_file)
                     if(input$adegeno_input %in% dartseq_formats){
                       dart_cols <- dart_getCols(input$adegeno_file$datapath)
                       updateSelectInput(session, "dartseq_markerid",
                                         choices = dart_cols)
                       updateSelectInput(session, "dartseq_chrom",
                                         choices = dart_cols)
                       updateSelectInput(session, "dartseq_position",
                                         choices = dart_cols)
                     }
                   })


      # Reactive to file reading
      geno_data <-
        eventReactive(input$load_geno,
          {
          # function reactive to adegeno_file input
          req(input$adegeno_file)
          req(input$adegeno_input)

          genotype_file <- input$adegeno_file$datapath
          geno_data <- NULL

          switch(input$adegeno_input,
                 hapmap = {
                   shinybusy::show_modal_spinner('fading-circle', text = 'Loading...')
                   tryCatch({
                     geno_data = cgiarGenomics::read_hapmap(path = genotype_file,
                                                            ploidity = input$ploidlvl_input)
                   }, error = function(e){
                     print(e)
                     shinyWidgets::show_alert(title = 'Error !!', text = 'Not a valid file format :-(', type = 'error')
                   })
                   shinybusy::remove_modal_spinner()
                 },
                 vcf = {
                   shinybusy::show_modal_spinner('fading-circle', text = 'Loading...')
                   tryCatch({
                     geno_data = cgiarGenomics::read_vcf(path = genotype_file,
                                                         ploidity = input$ploidlvl_input)
                   }, error = function(e){
                     print(e)
                     shinyWidgets::show_alert(title = 'Error !!', text = 'Not a valid file format :-(', type = 'error')
                   })
                   shinybusy::remove_modal_spinner()
                 },
                 dartseqpa = {
                   shinybusy::show_modal_spinner('fading-circle', text = 'Loading...')
                   tryCatch({
                     geno_data = cgiarGenomics::read_DArTSeq_PA( dart_path = genotype_file,
                                                                 marker_id = input$dartseq_markerid,
                                                                 chr_name = input$dartseq_chrom,
                                                                 pos_name = input$dartseq_position)
                   }, error = function(e){
                     print(e)
                     shinyWidgets::show_alert(title = 'Error !!', text = 'Not a valid file format :-(', type = 'error')
                   })
                   shinybusy::remove_modal_spinner()
                 },
                 dartsesnp = {
                   shinybusy::show_modal_spinner('fading-circle', text = 'Loading...')
                   tryCatch({
                     geno_data = cgiarGenomics::read_DArTSeq_SNP( dart_path = genotype_file,
                                                                  snp_id = input$dartseq_markerid,
                                                                  chr_name = input$dartseq_chrom,
                                                                  pos_name = input$dartseq_position)
                   }, error = function(e){
                     print(e)
                     shinyWidgets::show_alert(title = 'Error !!', text = 'Not a valid file format :-(', type = 'error')
                   })
                   shinybusy::remove_modal_spinner()
                 },
                 darttag = {
                   shinybusy::show_modal_spinner('fading-circle', text = 'Loading...')
                   tryCatch({
                     geno_data = cgiarGenomics::read_DArT_Tag(counts.file = input$genotype_file,
                                                              dosage.file = input$darttag_dosage_file$datapath,
                                                              ploidy = input$ploidlvl_input)
                   }, error = function(e){
                     print(e)
                     shinyWidgets::show_alert(title = 'Error !!', text = 'Not a valid file format :-(', type = 'error')
                   })
                   shinybusy::remove_modal_spinner()
                 },
                 {
                   print("This should not execute D;")
                 })
          return(geno_data)
        })

      output$summary_by_chrom <- renderTable({
        req(geno_data())
        geno_data <- geno_data()
        data.frame(
          chrom = unique(geno_data$meta[,'CHROM']),
          min_pos = aggregate(POS ~ CHROM, data = geno_data$meta, FUN = min)[,2],
          max_pos = aggregate(POS ~ CHROM, data = geno_data$meta, FUN = max)[,2],
          snps_count = aggregate(POS ~ CHROM, data = geno_data$meta, FUN = length)[,2]
        )

      })

    })
  }


# Util functions pending to move ------------------------------------------

# get the column names of dart csv files
dart_getCols <- function(path){
  # remove all rows staring with * character
  top_rows <- read.csv(path,
                       comment.char = "*",
                       nrows = 20)
  return(colnames(top_rows))
}
