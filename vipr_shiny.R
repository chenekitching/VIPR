library(shiny)
library(shinydashboard)
library(readr)
library(reactable)
library(dplyr)
library(DT)
library(flexdashboard)
library(fresh)
library(shinyBS)
library(ggplot2)
library(tidyr)
source("functions.R")

#Conflicts
conflicted::conflicts_prefer(shinydashboard::box)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(shinydashboard::valueBoxOutput)
conflicted::conflicts_prefer(shinydashboard::renderValueBox)
conflicted::conflicts_prefer(shinydashboard::valueBox)
conflicted::conflicts_prefer(DT::dataTableOutput)

#Get cmdline args
args <- commandArgs(trailingOnly = TRUE)

# Check if the file path is provided, otherwise, set a default value
file_path <- ifelse(length(args) > 0, args[1], "results/prioritised_file.txt")
build <- ifelse(length(args) > 0, args[2], "hg19")
# Adjust column names based on the genome build
if (build == "hg19"){
    af_colname = "AF"
    af_male_colname = "AF_male"
    af_female_colname = "AF_female"
    af_raw_colname = "AF_raw"
    af_afr_colname = "AF_afr"
    af_sas_colname = "AF_sas"
    af_amr_colname = "AF_amr"
    af_eas_colname = "AF_eas"
    af_nfe_colname = "AF_nfe"
}
if (build == "hg38"){
    af_colname = "gnomad40_genome_AF"
    af_male_colname = "gnomad40_genome_AF_male"
    af_female_colname = "gnomad40_genome_AF_female"
    af_raw_colname = "gnomad40_genome_AF_raw"
    af_afr_colname = "gnomad40_genome_AF_afr"
    af_sas_colname = "gnomad40_genome_AF_sas"
    af_amr_colname = "gnomad40_genome_AF_amr"
    af_eas_colname = "gnomad40_genome_AF_eas"
    af_nfe_colname = "gnomad40_genome_AF_nfe"
}

#load ml file
ml_file <- read_tsv(file_path) %>%
  mutate_at(c(af_colname, 
                   af_male_colname,
                   af_female_colname,
                   af_raw_colname, 
                   af_afr_colname, 
                   af_sas_colname, 
                   af_amr_colname, 
                   af_eas_colname, 
                   af_nfe_colname,
              'non_topmed_AF_popmax',
              'non_neuro_AF_popmax',
              'non_cancer_AF_popmax',
              'non_neuro_AF_popmax',
              'FATHMM_score',
              'GERP++_RS',
              "CADD_phred",
              'GenoCanyon_score',
              'LRT_score',
              'M-CAP_score',
              'MetaLR_score', 
              'MetaSVM_score', 
              'MutationTaster_score',
              'MutationAssessor_score',
              'PROVEAN_score', 
              'SIFT_score', 
              'SiPhy_29way_logOdds',
              'fathmm-MKL_coding_score',
              'integrated_fitCons_score', 
              'integrated_confidence_value',
              'phastCons100way_vertebrate',
              'phastCons30way_mammalian', 
              'phyloP100way_vertebrate',
              'phyloP30way_mammalian'),
            as.numeric) %>%
            rename(AF = all_of(af_colname),
           AF_male = all_of(af_male_colname),
           AF_female = all_of(af_female_colname),
           AF_raw = all_of(af_raw_colname),
           AF_afr = all_of(af_afr_colname),
           AF_sas = all_of(af_sas_colname),
           AF_amr = all_of(af_amr_colname),
           AF_eas = all_of(af_eas_colname),
           AF_nfe = all_of(af_nfe_colname)) %>%
  rm_missing_rows()

ml_file[ml_file == "."] <- NA
ml_file[ml_file == ""] <- NA



ui <- dashboardPage(
  dashboardHeader(title = "VIPR"),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    use_theme(mytheme),
    fluidRow(
      tabBox(
        # tags$head(
        #   tags$style(HTML(" #tabBox { height:90vh !important; } "))
        # ),
        id="tabs",
        # title = "tabBox",
        width = 12,
        tabPanel(title = "Overview",
                 value = "overview",
                 h4("Table of prioritised variants. Click on the ID of a variant to view detailed information."),
                 fluidRow(
                   shinydashboard::box(
                     title = "Prioritised variants",
                     solidHeader = TRUE,
                     collapsible = TRUE,
                     status = "primary",
                     width = 12,
                     height = 12,
                     DT::dataTableOutput("var_pr")
                   )
                 )
        ),
        tabPanel(title = "Individual analysis",
                 value = "indiv_var",
                 fluidRow( 
                   box(title = "Overview",
                       solidHeader = TRUE,
                       collapsible = TRUE,
                       width = 8,
                       status = "primary",
                       DT::dataTableOutput("overview_table")),
                   box(title = "VIPR pathogenicity prediction",
                       solidHeader = TRUE,
                       collapsible = TRUE,
                       status = "primary",
                       width = 4,
                       height = 5,
                       gaugeOutput("ml"))
                   ),
                 fluidRow(
                   valueBoxOutput("clinvar"),
                   valueBoxOutput("intervar"),
                   box(title = "ACMG criteria met",
                       solidHeader = TRUE,
                       collapsible = TRUE,
                       status = "primary",
                       icon("info", id = "icon_acmg",
                            class = "about-icon fa-pull-right"),
                       bsTooltip("icon_acmg", 
                                 "As determined by InterVar", 
                                 placement = "bottom", 
                                 trigger = "hover",
                                 options = NULL),
                       width = 4,
                       height = 5,
                       textOutput("acmg_crit")
                   )
                 ),
                 fluidRow(
                   box(title = "Functional prediction scores",
                       solidHeader = TRUE,
                       collapsible = TRUE,
                       status = "primary",
                       icon("info", id = "icon_func",
                            class = "about-icon fa-pull-right"),
                       bsTooltip("icon_func", 
                                 "Functional scores were scaled to allow graph view", 
                                 placement = "bottom", 
                                 trigger = "hover",
                                 options = NULL),
                       
                       plotOutput("func_preds"
                       ),
                       width = 4
                      # height = 10 
                       ),
                   box(title = "Conservation scores",
                       solidHeader = TRUE,
                       collapsible = TRUE,
                       status = "primary",
                       icon("info", id = "icon_cons",
                            class = "about-icon fa-pull-right"),
                       bsTooltip("icon_cons", 
                                 "Conservation scores were scaled to allow graph view", 
                                 placement = "bottom", 
                                 trigger = "hover",
                                 options = NULL),
                       plotOutput("cons_scores"),
                       width = 4
                      # height = 10
                       ),
                   
                   box(title = "Population frequency",
                       solidHeader = TRUE,
                       collapsible = TRUE,
                       status = "primary",
                       icon("info", id = "icon_pop",
                            class = "about-icon fa-pull-right"),
                       bsTooltip("icon_pop", 
                                 "Allele frequencies from gnomAD database", 
                                 placement = "bottom", 
                                 trigger = "hover",
                                 options = NULL),
                       DT::dataTableOutput("pop_freq"),
                       width = 4
                       #height = 10
                       )
                 )
                 
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  
  reactive_ml <- reactive({
    #ov_cols <- overview_cols(ml_file)
    overview_df <- ml_file %>%
      select(varID,
             Chr,
             Start,
             End,
             Gene.refGene,
             AAChange.refGene,
             avsnp150, 
             CLNALLELEID,
             CLNSIG,
             InterVar_automated,
             `.pred_P_LP`) %>%
      mutate(AAChange.refGene = gsub(",.*", "", AAChange.refGene)) %>%
      mutate(across('.pred_P_LP',
                    round, 2)) %>%
      rename(
        #"Variant" = "AAChange.refGene",
        "RSID" = "avsnp150",
        "ClinVar ID" = "CLNALLELEID" ,
        "ClinVar classification" = "CLNSIG",
        "InterVar classification" = "InterVar_automated"
        # "Probability of pathogenicity" = ".pred_P_LP"
      ) %>%
      mutate(across(c(`ClinVar ID`,
                      RSID),
                    ~ replace(.x, is.na(.x), "")))
    
    #hyperlinks
    overview_df$`ClinVar ID` <- paste0("<a href='https://www.ncbi.nlm.nih.gov/clinvar/?term=",
                                       overview_df$`ClinVar ID`,
                                       "'target='_blank'>",
                                       overview_df$`ClinVar ID`,
                                       "</a>")
    
    overview_df$RSID <- paste0("<a href='https://www.ncbi.nlm.nih.gov/snp/?term=",
                               overview_df$RSID,
                               "'target='_blank'>",
                               overview_df$RSID,
                               "</a>")
    
    
    return(overview_df)
  })
  
  #table with probability predictions
  output$var_pr <- DT::renderDataTable({
    
    DT::datatable(
      reactive_ml(),
      # rownames = FALSE,
      options = list(order = list(list(11, "desc")),
                     scrollX = TRUE,
                     #Change displayed column name, underlying name stays same
                     columns = list(
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       list(title = "Probability of pathogenicity")
                     )),
      # colnames = c("Probability of pathogenicity" = 11),
      selection = list(mode = "single", target = "cell",
                       selectable = rbind(cbind(1:nrow(reactive_ml()), rep(1, nrow(reactive_ml()))))),
      # rownames = FALSE,
      escape = FALSE) %>% 
      formatStyle(
      'varID',
      cursor = 'pointer',
      color = 'teal'
    ) %>%
      color_gradient(".pred_P_LP")
  })
  
  #change to individual page when variant clicked
  observeEvent(input$var_pr_cells_selected, {
    # alternative: input$dt1_cells_selected
    req(input$var_pr_cells_selected)
    updateTabsetPanel(session, inputId = "tabs", selected = "indiv_var")

  })

  #value of clicked cell (variant)
  cell_clicked <- reactive({
    as.data.frame(reactive_ml())[input$var_pr_cells_selected] %>%
      as.character()
  })
  
  output$overview_table <- DT::renderDataTable({

    DT::datatable(
      indiv_table(as.character(cell_clicked()),
                  ml_file),
      options = list(dom = 't'),
      selection = "none",
      rownames= FALSE,
      escape = FALSE
    )

  })

  
  #Clinvar prediction
  clinvar_reactive <- reactive({
    get_class(as.character(cell_clicked()),
              ml_file,
              "CLNSIG")
  })
  
  #Intervar prediction
  intervar_reactive <- reactive({
    get_class(as.character(cell_clicked()),
              ml_file,
              "InterVar_automated")
  })
  
  #Colour of clinvar/intervar boxes
  col_reactive_c <- reactive({
    val_box_col(clinvar_reactive())
  })
  
  col_reactive_i <- reactive({
    val_box_col(intervar_reactive())
  })
  
  output$clinvar <- renderValueBox({
    valueBox(value = tags$p(clinvar_reactive(),
                            style = "font-size: 50%;"),
             "Clinvar",
             color = col_reactive_c(),
             icon = tags$i(class = "fas fa-list", style="font-size: 50px")
    )
  })
  
  output$intervar <- renderValueBox({
    valueBox(value = tags$p(intervar_reactive(),
                            style = "font-size: 50%;"),
             "Intervar",
             color = col_reactive_i(),
             icon = tags$i(class = "fas fa-magnifying-glass-chart", style="font-size: 50px")
    )
  })
  
  
  #Pathogenicity gauge
  prob_path_reactive <- reactive({
    get_class(as.character(cell_clicked()),
              reactive_ml(),
              ".pred_P_LP") %>%
      as.numeric() * 100
    
  })
  
  output$ml <- renderGauge({
    gauge(value = prob_path_reactive(),
          min = 0,
          max = 100,
          symbol  = '%',
          label = "Probability of pathogenicity",
          gaugeSectors(success = c(0,40),
                       warning = c(41, 60),
                       danger = c(61, 100)))
  })
  
  acmg_reactive <- reactive({
    get_criteria(as.character(cell_clicked()),
                 ml_file)
    
  })
  
  
  
  output$acmg_crit <- renderText({
    acmg_txt <- acmg_reactive()
    validate(
      need(nchar(acmg_txt) > 0, "No data to show")
    )
    acmg_txt
    })
  
  
  
  
  func_reactive <- reactive({
    func <- func_pred(as.character(cell_clicked()),
                      ml_file)
    
    scaled <- unlist(func) %>%
      scale() %>% 
      t() %>%
      as.data.frame()
    
      scaled %>%
      pivot_longer(cols = everything(),
                   names_to = "Functional prediction",
                   values_to = "Score")

    
  })
  
  output$func_preds <- renderPlot({
    func_reactive() %>%
      ggplot(aes(x = `Functional prediction`,
                 y = Score)) +
      geom_segment(aes(x = `Functional prediction`, xend = `Functional prediction`,
                       y = min(Score,
                               na.rm = TRUE), yend = Score),
                   color = "skyblue") +
      geom_point( color="blue", size=4, alpha=0.6) +
      theme_light() +
      coord_flip() +
      theme(text = element_text(size = 20)) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()
      )
  })
  
  cons_reactive <- reactive({
    cons <- cons_scores(as.character(cell_clicked()),
                        ml_file)
    
    scaled <- unlist(cons) %>%
      scale() %>% 
      t() %>%
      as.data.frame()
    
    scaled %>%
      pivot_longer(cols = everything(),
                   names_to = "Conservation score",
                   values_to = "Score")
    
  })
  
  output$cons_scores <- renderPlot({
    cons_reactive() %>%
      ggplot(aes(x = `Conservation score`,
                 y = Score)) +
      geom_segment(aes(x = `Conservation score`, xend = `Conservation score`,
                       y = min(Score,
                               na.rm = TRUE), yend = Score),
                   color = "skyblue") +
      geom_point( color="blue", size=4, alpha=0.6) +
      theme_light() +
      coord_flip() +
      theme(text = element_text(size = 20)) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()
      )
  })
  
  pop_freq_reactive <- reactive({
    freq <- pop_freq(as.character(cell_clicked()),
                     ml_file)
  })
  
  output$pop_freq <- DT::renderDataTable({
    
    pop_freq_tbl <- DT::datatable(pop_freq_reactive(),
                                  options = list(scrollX = TRUE,
                                                 dom = "t"),
                                  rownames = TRUE,
                                  colnames = rep("", ncol(pop_freq_reactive())))
    na_rows <- sum(!complete.cases(pop_freq_reactive()))
    validate(
      need(na_rows < 12, "No data to show")
    )
    pop_freq_tbl
    
    
  })
  
}

shinyApp(ui, server)