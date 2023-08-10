################################################
################# GBA1 PORTAL ##################
################################################

################## LIBRARIES ################
library(shiny)
library(plyr)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(plotly)
library(tippy)
library(r3dmol)
library(DT)
library(readr)
library(tidyverse)
library(vembedr)
library(RColorBrewer)
library(shiny.i18n)
library(shinyjs)
################################################

source("dependence/ui_dependence.R")

################# UI #################
ui <- 
  navbarPage(#fluid = FALSE, 
    windowTitle = "GBA1 Portal",   #c
    id = "TabDisplay",
    theme = "mytheme2.css",
    title = p(icon("dna"), landing_portal_title), #c
    
    header = tagList(useShinydashboard()),
    
    tabPanel(
      
      tags$style(HTML(paste0("
                      .navbar-default .navbar-brand {
                      background-color: ",schema_color_light,";
                      color: black
                    }"))),
      tags$style(HTML(paste0("
                      .navbar-default {
                      background-color: ",schema_color_light,";
                    }"))),
      
      
      title = "Welcome",
      value = "welcomeTab",
      tabName = "welcomeTab",
      
      div(style = "font-size:100%", fluidRow(
        column(1.5,
               style = "background-color: white; color: #676767",
               align = "center",
               br(), br(),
               img(src = landing_bannername1, width = "100%")   #insert your banner, saved in the www-folder
        )
      ),
      
      div(style = 'background-color: #F8FCFE',
          fluidRow(
            
            # Basic Information Tab
            div(width = "100%", 
                column(3,offset = 2, align = "center",
                       br(),
                       panel(width = 12, 
                             status = "success",
                             heading = "",
                             h2(tags$i(class = "fa fa-dna", style = "color:#676767")),
                             br(),
                             div(landing_tab1,  #c
                                 br(), br(),
                                 actionBttn(
                                   inputId = "infoBtn",
                                   label = "Basic Information",
                                   color = "success",
                                   block = TRUE,
                                   size = "md",
                                   style = "stretch")),
                             style =  "background-color: #f3faf4;",   #set color of the tab
                       ) 
                ),
                
                # Variant Analysis Tab
                column(3,  align = "center",
                       br(),
                       div(panel(
                         status = "info",
                         heading = "",
                         h2(tags$i(class = "fa fa-code-branch", style = "color:#676767")),
                         br(),
                         div(landing_tab2,br(), 
                             br(),
                             actionBttn(
                               inputId = "variantBtn",
                               label = "Variant Analysis",
                               color = "royal",
                               block = TRUE,
                               style = "stretch"
                             )),
                         style = "background-color: #f9f1fa;",    #set color of the tab
                       ))
                ),
                
                # Research Tab
                column(3,  align = "center",
                       br(),
                       div(panel(
                         status = "danger",
                         heading = "",
                         h2(tags$i(class = "fa fa-microscope", style = "color:#676767")),
                         br(),
                         div(landing_tab4,br(), 
                             br(),
                             actionBttn(
                               inputId = "researchBtn",
                               label = "Research",
                               color = "danger",
                               block = TRUE,
                               style = "stretch"
                             )),
                         style = "background-color: #fdf0f1;",
                         
                       ))
                ),
                
                column(12,
                       style = "background-color: white; color: #676767",
                       align = "center",
                       br(), br(),
                       img(src = landing_bannername, width = "25%")   #insert your banner, saved in the www-folder
                )
            )),
          fluidRow(
            column(4, offset = 2, align = "center",
                   
                   p("Visit our other Portals:",style = "font-size:14px;",align = "center",style = "font-size:15px;"),
                   p(shiny::a("SLC6A1-Portal", href="http://slc6a1-portal.broadinstitute.org/", target = '_blank'), "& ",
                     shiny::a("GRIN-Portal", href="http://grin-portal.broadinstitute.org/", target = '_blank'),
                     style = "font-size:15px;")
            ),
            column(4, align = "center",
                   
                   p("You want to join the project or provide feedback?", align = "center",style = "font-size:15px;"),
                   p("Please ", shiny::a("contact us!", href=contact_us), style = "font-size:15px;", align = "center"),
                   br()
            )))
      )), # end tab Panel 
    
    ##### BASIC INFORMATION #####
    tabPanel(title = "Basic Information", value = "infoTab",
             useShinyjs(),
             
             fluidRow(
               column(10, offset = 1,
                      panel(
                        heading = div("Basic Information", style = "color : white"),    #c
                        status = "success",
                        
                        #Introduce your gene/genes of interest and display basic information, add as much tabs as you want
                        ## Gene ##
                        tabsetPanel(
                          type  = "tabs",
                          id = "displaygene",
                          tabPanel(
                            title = p(em(gene1)),    #c
                            value = "basic_GBA1",
                            panel(status = "default", heading = gene1_info_header,   
                                  fluidRow(column(4, box(width = 12,
                                                         
                                                         title = p("History of ",em(gene1), " Research"),   #c
                                                         
                                                         timelineBlock(
                                                           reversed = FALSE,
                                                           width = 12,
                                                           
                                                           #Add as many blocks as desired for new publications
                                                           timelineLabel(1952, color = "teal"),
                                                           timelineItem(
                                                             title = div(strong("Discovery of the action potential.")),
                                                             icon = icon("dna"),
                                                             color = "olive",
                                                             #ADD link to the publication
                                                             time = shiny::a("Hodgkin et al.", href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1392413/", target = '_blank'),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(1978, color = "teal"),
                                                           timelineItem(
                                                             title = div("First description of Dravet Syndrome"),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = "Dravet",
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(1997, color = "teal"),
                                                           timelineItem(
                                                             title = div("Defining genetic epilepsy with febrile seizures plus (GEFS+)"),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = shiny::a("Scheffer et al.", href="https://pubmed.ncbi.nlm.nih.gov/9126059/", target = '_blank'),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(1998, color = "teal"),
                                                           timelineItem(
                                                             title = div("Seizure worsening due to sodium channel blockers in Dravet syndrome"),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = shiny::a("Guerrini et al.", href="https://pubmed.ncbi.nlm.nih.gov/9596203/", target = '_blank'),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(2000, color = "teal"),
                                                           timelineItem(
                                                             title = div(HTML(paste0("First association of ", em("SCN1A")," with GEFS+"))),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = shiny::a("Escayg et al.", href="https://pubmed.ncbi.nlm.nih.gov/10742094/", target = '_blank'),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(2000, color = "teal"),
                                                           timelineItem(
                                                             title = div("Stiripentol treatment for Dravet syndrome"),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = shiny::a("Chiron et al.", href="https://pubmed.ncbi.nlm.nih.gov/11089822/", target = '_blank'),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(2001, color = "teal"),
                                                           timelineItem(
                                                             title = div(HTML(paste0("De novo mutations in ",em("SCN1A")," cause Dravet syndrome"))),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = shiny::a("Claes et al.", href="https://pubmed.ncbi.nlm.nih.gov/11359211/", target = '_blank'),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(2006, color = "teal"),
                                                           timelineItem(
                                                             title = div(HTML(paste0("Selective loss of GABAergic interneurons results in hyperexcitability of pyramidal neurons, leading to seizures in ",em("SCN1A")," mice"))),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = shiny::a("Yu et al.", href="https://pubmed.ncbi.nlm.nih.gov/16921370/", target = '_blank'),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(2006, color = "teal"),
                                                           timelineItem(
                                                             title = div("The misconception of vaccine encephalopathy associated with Dravet syndrome"),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = shiny::a("Berkovic et al.", href="https://pubmed.ncbi.nlm.nih.gov/16713920/", target = '_blank'),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel("2007/2009/2012/2017", color = "teal"),
                                                           timelineItem(
                                                             title = div(HTML(paste0("Genotype phenotype associations in large ",em("SCN1A"), " cohorts"))),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = p(shiny::a("Harkin et al.", href="https://pubmed.ncbi.nlm.nih.gov/17347258/", target = '_blank'),
                                                                      br(),
                                                                      shiny::a("Depienne et al.", href="https://pubmed.ncbi.nlm.nih.gov/18930999/", target = '_blank'),
                                                                      br(),
                                                                      shiny::a("Brunklaus et al.", href="https://pubmed.ncbi.nlm.nih.gov/22719002/", target = '_blank'),
                                                                      br(),
                                                                      shiny::a("Ishii et al.", href="https://pubmed.ncbi.nlm.nih.gov/28012175/", target = '_blank')),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(2011, color = "teal"),
                                                           timelineItem(
                                                             title = div("First image of the crystal structure of a voltage-gated sodium channel"),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = shiny::a("Payandeh et al.", href="https://pubmed.ncbi.nlm.nih.gov/21743477/", target = '_blank'),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(2012, color = "teal"),
                                                           timelineItem(
                                                             title = div("Autistic-like behaviour in Scn1a+/- mice and rescue by enhanced GABA-mediated neurotransmission"),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = shiny::a("Han et al.", href="https://pubmed.ncbi.nlm.nih.gov/22914087/", target = '_blank'),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(2013, color = "teal"),
                                                           timelineItem(
                                                             title = div("Dravet iPSC-derived neuron work"),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = shiny::a("Liu et al.", href="https://pubmed.ncbi.nlm.nih.gov/23821540/", target = '_blank'),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(2017, color = "teal"),
                                                           timelineItem(
                                                             title = div("Cannabidiol treatment for Dravet syndrome"),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = shiny::a("Devinsky et al.", href="https://pubmed.ncbi.nlm.nih.gov/28813226/", target = '_blank'),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(2019, color = "teal"),
                                                           timelineItem(
                                                             title = div("Fenfluramine treatment for Dravet syndrome"),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = shiny::a("Lagae et al.", href="https://pubmed.ncbi.nlm.nih.gov/31862249/", target = '_blank'),
                                                             border = FALSE,
                                                           ),
                                                           timelineLabel(2020, color = "teal"),
                                                           timelineItem(
                                                             title = div("Antisense oligonucleotide treatment in Dravet mouse model"),
                                                             icon = icon("user"),
                                                             color = "aqua",
                                                             #ADD link to the publication
                                                             time = shiny::a("Han et al.", href="https://pubmed.ncbi.nlm.nih.gov/32848094/", target = '_blank'),
                                                             border = FALSE,
                                                           )
                                                           
                                                           
                                                         )
                                  )),
                                  column(8, box(title=basic_text_title_1, width = 12,   #c
                                                Gene1_basic_text)),
                                  column(8,
                                         box(title = p(basic_clinical_info_title),
                                             width = 12,
                                             
                                             tabsetPanel(
                                               #Example of a factor phenotype representation
                                               tabPanel(title = "Parkinson phenotype",
                                                        br(),
                                                        div(p("Number of patients with data in this category =" ,Patient_data.df %>% filter(Gene == gene1) %>% nrow(), ""),style = "font-size:15px;color:black;"),
                                                        br(),
                                                        fluidRow(align = "center",
                                                                 div(title = paste0("Phenotypes available for " ,Patient_data.df %>% filter(Gene == gene1) %>% nrow(), " patients"), width = 12),
                                                                 column(12, align = "center",plotOutput(outputId = "basic_legend_plot1_1", height = 20, width = 300)),
                                                                 column(12, box(title=div(basic_plot_title_fac1, style = "font-size: 15px"), width=12, plotlyOutput(outputId = "Phenotype_fac1_gene1"))),
                                                                 
                                                        )),
                                               
                                               #Example of a numeric phenotype representation
                                               tabPanel(title = "Gaucher phenotype",   #c
                                                        br(),
                                                        div(p("Number of patients with data in this category =" ,Patient_data.df %>% filter(Gene == gene1) %>% nrow(), ""),style = "font-size:15px;color:black;"),
                                                        br(),
                                                        fluidRow(align = "center",
                                                                 column(12, align = "center",plotOutput(outputId = "basic_legend_plot1_2", height = 20, width = 300)),
                                                                 column(12, box(title=div(basic_plot_title_num1, style = "font-size: 10px"), width=12, plotlyOutput(outputId = "Phenotype_fac2_gene1"))),
                                                                 
                                                                 # column(12, align="center", p(basic_abbreviations1, style=sub_style))
                                                        ))
                                      
                                             )
                                         ))
                                  ))
                          )
                          
                        ))))
             #          
    ), # end basic info tab panel

    # ##### VARIANT ANALYSIS #####
    
    tabPanel(title = "Variant Analysis", value = "variantTab",
             fluidRow(
               column(10, offset = 1,
                      panel(heading = "Analyse your variants", status = "info",
                            
                            fluidRow(column(12,style='padding:30px;',
                                            fluidRow(panel(status="default", heading = "Enter variant",
                                                           column(3,
                                                                  pickerInput(
                                                                    inputId = 'var_gene',
                                                                    label =  h5(strong(var_possible_genes_title)),
                                                                    width = "100%",
                                                                    choices = var_possible_genes,
                                                                    options = list(`style` = "default")
                                                                  )
                                                           ),
                                                     
                                                           column(3,
                                                                  numericInputIcon(
                                                                    inputId = "search_cDNA_pos",
                                                                    label = h5(strong("cDNA Position")),
                                                                    min = 1,
                                                                    max = 10000,
                                                                    value = 23,
                                                                    width = "100%"
                                                                  ),
                                                                  pickerInput(
                                                                    inputId = "search_Allele",
                                                                    label = "Ref",
                                                                    choices = c("G", "A", "C", "T"),
                                                                    selected = "C"
                                                                  ),
                                                                  pickerInput(
                                                                    inputId = "search_cDNA_alt",
                                                                    label = "Alt",
                                                                    choices = c("G", "A", "C", "T"),
                                                                    selected = "T"
                                                                  ),
                                                                  actionButton(inputId = "search_var_c", label = "Search")),
                                                           column(3,
                                                                  numericInputIcon(
                                                                    inputId = "search_AA_pos",
                                                                    label = h5(strong("Amino Acid Position")),
                                                                    min = 1,
                                                                    max = 10000,
                                                                    value = 8,
                                                                    width = "100%"
                                                                  ),
                                                                  pickerInput(
                                                                    inputId = "search_AA_ref",
                                                                    label = "Ref",
                                                                    choices = sort(unique(master.df$AA_ref)),
                                                                    selected = "Arg"
                                                                  ),
                                                                  pickerInput(
                                                                    inputId = "search_AA_alt",
                                                                    label = "Alt",
                                                                    choices = sort(unique(master.df$AA_alt)),
                                                                    selected = "Trp"
                                                                  ),
                                                                  actionButton(inputId = "search_var_p", label = "Search")),
                                                           column(2,
                                                                  pickerInput(
                                                                    inputId = "get_var_type",
                                                                    label = h5(strong("Variant Consequence")),
                                                                    choices = "Missense",
                                                                    selected = "Missense"))
                                            )))),
                            
                            fluidRow(
                              
                              column(12,style='padding:30px;',
                                     fluidRow(
                                       panel(
                                         status="default", heading = "Variant Information",
                                         fluidRow(
                                           div(width = "100%",valueBoxOutput("geneBox1")),
                                           div(width = "100%",valueBoxOutput("geneBox2")),
                                           div(width = "100%",valueBoxOutput("geneBox3"))
                                           
                                         ))))),

                            fluidRow(
                              column(12,style='padding:30px;',
                                     fluidRow(
                                       panel(
                                         status="default", heading = "Variant Information",
                                         tabsetPanel(
                                           tabPanel("Patient Information",br(),
                                                    box(title = p(var_patient_info_title,
                                                                  tippy(icon("question-circle"),
                                                                        tooltip = h6(HTML(paste0("DOI of Articles:")),
                                                                                     HTML(paste0("<ul><li>Article 1: 10.1093/brain/awv179</li>",
                                                                                                 "<li>Article 2: 10.1002/acn3.51164</li>",
                                                                                                 "<li>Article 3: 10.1002/mgg3.267</li>")),
                                                                                     
                                                                                     align = "left"),
                                                                        animation = "scale",
                                                                        theme = "light")
                                                    ),
                                                    width = 12,
                                                    DT::dataTableOutput(outputId = "patientTable"),
                                                    br(), p(var_patient_info_abb, style=sub_style, align = "center")))
                                     ))))),
                            
                            fluidRow(
                              column(12,style='padding:30px;',
                                     fluidRow(
                                       panel(
                                         status="default", heading = "Custom variant analysis",
                                         tabsetPanel(
                                           tabPanel("Comparative Information",
                                                    br(),
                                                    h4("Compare your variant with other variants in the ", em("GBA1", style="font-style: italic;"), " Portal"),
                                                    br(),
                                                    radioGroupButtons(
                                                      inputId = "compareButtons",
                                                      label = "Variants with the same:",
                                                      choices = c("Domain","Variant Type"),
                                                      justified = TRUE,
                                                      status = "default",
                                                      checkIcon = list(yes = icon("ok",lib = "glyphicon"))
                                                    ),
                                                    materialSwitch(inputId = "other_genes"),
                                                    br(),
                                                    div(width = "100%", plotlyOutput("comparePlot")),
                                                    br()
                                           )
                                         ))
                                     ))))))), # end variant analysis tab

    #         ##### RESEARCH #####
    tabPanel(title = "Research", value = "researchTab",
             
             fluidRow(
               column(10, offset = 1,
                      panel(heading = "Analyse your variants", status = "danger",
                            tabsetPanel(id = "Researchfilter",
                                        tabPanel(title = "Filter variants",
                                                 status = "default",
                                                 
                                                 fluidRow(
                                                   column(12,style='padding:15px;',
                                                          panel(heading= p("Filter variants form the ", em("GBA1", style="font-style: italic;"), "Portal"),
                                                          status="default",
                                                          fluidRow(
                                                            column(12,
                                                                   selectizeGroupUI(
                                                                     id = "research-filters",
                                                                     btn_label = "Reset filters",
                                                                     params = list(
                                                                       gene = list(
                                                                         inputId = "Gene",
                                                                         placeholder="all",
                                                                         title = p(strong("Gene")),
                                                                         choices = var_possible_genes,
                                                                         multiple = TRUE
                                                                       ),
                                                                       varcons = list(
                                                                         inputId = "Vartype",
                                                                         title = p(strong("Variant Type"),
                                                                                   tippy(icon("question-circle"),
                                                                                         tooltip = h6(HTML(paste0("PTV: Protein truncatiing variant")),
                                                                                                      align = "left"),
                                                                                         animation = "scale", 
                                                                                         theme = "light")),
                                                                         placeholder="all",
                                                                         choices = c("Missense", "PTV","frameshift insertion","synonymous")
                                                                       ),
                                                                       aachange = list(
                                                                         inputId = "AA_alt",
                                                                         title = p(strong("Amino Acid Change")),
                                                                         placeholder="all",
                                                                         choices = unique(Patient_data.df$AA_alt)
                                                                       ),
                                                                       domain = list(
                                                                         inputId = "Domain",
                                                                         title = p(strong("Domain")),
                                                                         placeholder="all",
                                                                         choices = unique(Patient_data.df$Domain)
                                                                       ),
                                                                 
                                                                       phenotype = list(
                                                                         inputId = "Phenotype",
                                                                         title = p(strong("Primary phenotype PD")),
                                                                         placeholder="all",
                                                                         choices = c("Pathogenic", "Pathogenic/Pathogenic", "Risk", "Unknown"),
                                                                         multiple = TRUE),
                                                                       
                                                                       phenotype_GD= list(
                                                                         inputId = "phenotype_GD",
                                                                         title = p(strong("Primary phenotype GD")),
                                                                         placeholder="all",
                                                                         choices = c("Benign", "Likely benign", "Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic", "Unknown"),
                                                                         multiple = TRUE),
                                                                       
                                                                       ID_x = list(
                                                                         inputId = "ID_x",
                                                                         title = p(strong("Severity of PD")),
                                                                         placeholder="all",
                                                                         choices = c(" Unknown", "Mild","Risk","Mild/severe",'Severe'),
                                                                         multiple = TRUE),
                                                                       
                                                                       ID = list(
                                                                         inputId = "ID",
                                                                         title = p(strong("Severity of GD")),
                                                                         placeholder="all",
                                                                         choices = c(" Unknown", "Mild","Severe (Null)","Risk Variant",'Severe'),
                                                                         multiple = TRUE)
                                                                      
                                                                     )
                                                                   ))
                                                          )))),
                                                 
                                                 fluidRow(
                                                   
                                                   column(12,style='padding:15px;',
                                                          panel(status = "default",heading = "Custom variant exploration",
                                                                div(textOutput(outputId = "filtered_n"),br(),"  "),
                                                                tabsetPanel(
                                                                  tabPanel(
                                                                    "Genotype Interface",
                                                                    fluidRow(
                                                                      column(
                                                                        12,
                                                                        br(),
                                                                        p(
                                                                          strong("Selected variants are displayed in 2D (lolliplot) and 3D (protein structure).")
                                                                        ),
                                                                        
                                                                        fluidRow(column(
                                                                          8,
                                                                          materialSwitch(
                                                                            inputId = "gnomad_m",
                                                                            label = "Reference population missense variants (gnomAD)",
                                                                            status = "primary",
                                                                            right = T,
                                                                            inline = T)),
                                                                          
                                                                        column(12,plotOutput(outputId = "Genotype_legend_plot", height = 30, width = 450),
                                                                               p("PTV: Protein truncating variant"))
                                                                        ),
                                                                        
                                                                        addSpinner(plotlyOutput(outputId = "Genotype_overview_plot"), color =
                                                                                     spinner_color),
                                                                        br(),
                                                                        p(research_geno_transcripts, align="center", style=sub_style),
                                                                        br()
                                                                      )),
                                                                    
                                                                    fluidRow(
                                                                      
                                                                      column(
                                                                        4,
                                                                      
                                                                        p(h4("Protein Data Bank (PDB)")),
                                                                        div(align = "center", style = sub_style),
                                                                        addSpinner(color = spinner_color,
                                                                                   r3dmolOutput(
                                                                                     outputId = "threeDmolGene1",
                                                                                     width = "100%",
                                                                                     height = "400px"
                                                                                   )),
                                                                        div("UniProt: ", align="center", style=sub_style,
                                                                            shiny::a("6TJK", href="https://www.rcsb.org/structure/6TJ",  target="_blank")),
                                                                      ),
                                                                      column(
                                                                        4,
                                                                        
                                                                        p(h4("Alphafold")),
                                                                        div(align = "center", style = sub_style),
                                                                        addSpinner(color = spinner_color,
                                                                                   r3dmolOutput(
                                                                                     outputId = "threeDmolGene2",
                                                                                     width = "100%",
                                                                                     height = "400px"
                                                                                   )),
                                                                        div("UniProt: ", align="center", style=sub_style,
                                                                            shiny::a("42 structures in PDB for P04062", href="https://alphafold.ebi.ac.uk/entry/P04062", target="_blank")),
                                                                      )
                                                                    )
                                                                  ),
                                                                  
                                                                  tabPanel("Phenotype Interface",br(),
                                                                           fluidRow(
                                                                             column(12,align="justify", plotlyOutput("research_phenotype1"))
                                                                             ),
                                                                           fluidRow(
                                                                             column(6, plotlyOutput("research_phenotype2"), align="center",
                                                                             ),
                                                                             column(6, plotlyOutput("research_phenotype3")),
                                                                           
                                                                             column(6, plotlyOutput("research_phenotype4")),
                                                                             column(6, plotlyOutput("research_phenotype5")),
                                                                             column(6, plotlyOutput("research_phenotype6")),
                                                                             column(6, plotlyOutput("research_phenotype7"))
                                                                           ))
                                                                  ),
                                                                br()
                                                          )
                                                   )))
                                        
                            )
                      )))
    ), # end research tab

    ##### ABOUT #####
    
    tabPanel(title = "About", value = "aboutTab",
             
             fluidRow(
               column(10, offset = 1,
                      panel(heading = "About", status = "primary",
                            tabsetPanel(
                              tabPanel("General information",
                                       panel(heading ="Portal version", status = "default",
                                             p(strong("This is the alpha version of the", em("GBA1", style="font-style: italic;"),"Portal.")),
                                             p("The", em("GBA1", style="font-style: italic;"),"Portal is a coalition of investigators seeking to aggregate and harmonize data generated to study",
                                               em("GBA1", style="font-style: italic;"),"related disorders, and to make summary data interactively accessible for the wider scientific community,
                                    while providing educational resources for everyone. The goals of this project are: "),br(),
                                      
                                             
                                             fluidRow(column(6,
                                                             p(HTML("<ul><li> Providing information on <em>GBA1</em> related disorders")),
                                                             p(HTML("</li><li> Supporting research on <em>GBA1</em> related disorders")),
                                                             p(HTML("</li><li> Providing support in variant interpretation and classification"))
                                             )),
                                             br(),
                                             footer = div("The", em("GBA1", style="font-style: italic;"),"Portal is an ongoing project and interested collaborators are invited to reach out to join the project.")),
                                       panel(heading = "Teams and People", status = "default",
                                             p("The current version of the", em("GBA1", style="font-style: italic;"),"Portal has been developed by an international team of researchers and clinicians:"),
                                             fluidRow(
                                               
                                               column(12, panel(heading = "Team Leaders",
                                                                
                                                                p(strong("Patrick May"), "(Luxembourg Centre for Systems Biomedicine, Luxembourg):",em("GBA1"),"Clinical & genetic data"),
                                                                p(strong("Rejko KRÜGER"), "(Luxembourg Centre for Systems Biomedicine, Luxembourg):",em("GBA1"),"Clinical & genetic data")    
                                               )),
                                               
                                               column(2, panel(style="height: 230px;",heading = "Clinical & Genetic Data",
                                                               div(style="height: 100%;",
                                                                   p(strong("Rejko KRÜGER")),
                                                                   p(strong("Sinthuja Pachchek Peiris")))
                                               )),
                                               
                                               column(2, panel(style="height: 230px;",heading = "Molecular Data",
                                                               div(style="height: 100%",p(strong("Rejko KRÜGER")))
                                               )),
                                               
                                               
                                               column(2, panel(style="height: 230px;",heading = "Web Development",
                                                               div(style="height: 100%",
                                                                   p(strong("Danial Tabbakh")))
                                               )),
                                               
                                               
                                               column(2, panel(style="height: 230px;",heading = "Bioinformatics",
                                                               div(style="height: 100%",
                                                                   p(strong("Sinthuja Pachchek Peiris")),
                                                                   p(strong("Danial Tabbakh")))
                                               ))
                                           )),
                                       
                                       panel(heading = "Impressum", status = "default", p("We object to any commercial use and disclosure of data."),
                                             p(strong("Copyright and use:"), "The authors grants you the right of use to make a private copy for personal purposes.
                                        However, you are not entitled to change and/or pass on the materials or publish them yourself.
                                        Upon request, you will receive free information about the personal data stored about you.
                                        To do so, please contact the administrator."),
                                             p(strong("No liability:"), "The contents of this web project have been carefully checked and created to
                                        the best of our knowledge. But for the information presented here is no claim to completeness,
                                        timeliness, quality and accuracy. No responsibility can be accepted for any damage caused by reliance
                                        on or use of the contents of this website."))
                              ),
                              tabPanel("Terms and Data Information",
                                       panel(heading = "Terms of Use", status = "default",
                                             p(about_terms_of_use,
                                               shiny::a(href=contact_us, "Contact us"),
                                               "that we can improve.")),
                                       panel(heading = "Data Generation", status = "default",
                                             about_data_generation))
                            ))))
             
    )
  ) # end ui