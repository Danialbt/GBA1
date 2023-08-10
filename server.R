################################################
################# GBA1 PORTAL ##################
################################################

################## LIBRARIES ################

library(shiny)
library(plyr)
library(plotly)
library(readr)
library(r3dmol)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(tidyverse)
library(seqinr)
library(bio3d)
library(RColorBrewer)
library(shiny.i18n)
library(shinyjs)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bio3d)
library(readxl)
library(patchwork)

################################

source("dependence/functions.R")
source("dependence/server_data.R")

############### SERVER ###############
options(shiny.debug = TRUE)
server <- function(input, output,session) {

  ############### PANEL DECISION ###############
  observeEvent(input$infoBtn, {
    updateTabsetPanel(session, "TabDisplay", selected = "infoTab")
  })
  
  observeEvent(input$variantBtn, {
    updateTabsetPanel(session, "TabDisplay", selected = "variantTab")
  })
  
  observeEvent(input$researchBtn, {
    updateTabsetPanel(session, "TabDisplay", selected = "researchTab")
  })

  ####Basic Information graphs
  output$Phenotype_fac1_gene1 <- renderPlotly({

    data.df <- res_mod() %>%filter(!is.na(AA_pos)) %>%
      rename(Phenotype1 = phenotype) %>%
      group_by(Gene, Phenotype1) %>% summarise(n = n(), .groups = "drop")  %>% ungroup()

    validate(
      # Old Code
      need(nrow(data.df) > 0, 'No data exists, for selection!')
      # With parentheses fix
    )
    Gene_colors1 <- c("#2ca02c")


    data.df %>%
      plot_ly(
        x=~Phenotype1, y=~n, split=~Gene, type="bar", color=~Gene, alpha = 0.8, colors= Gene_colors1) %>%
      layout(title=research_phenotype2_title ,font=plotly_font, yaxis = list(showline = T, title = "Number of patients"),
             margin = list(t =50),
             xaxis = list(title = "",  showline = T,tickangle = 45)
      ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  })

  output$Phenotype_fac2_gene1 <- renderPlotly({

    data.df <- res_mod() %>%filter(!is.na(AA_pos)) %>%
      rename(Phenotype2 = phenotype_GD) %>%
      group_by(Gene, Phenotype2) %>% summarise(n = n(), .groups = "drop")  %>% ungroup()

    validate(
      # Old Code
      need(nrow(data.df) > 0, 'No data exists, for selection!')
      # With parentheses fix
    )
    Gene_colors2 <- c("#9467bd")


    data.df %>%
      plot_ly(
        x=~Phenotype2, y=~n, split=~Gene, type="bar", color=~Gene, alpha = 0.8, colors= Gene_colors2) %>%
      layout(title="Phenotype GD" ,font=plotly_font, yaxis = list(showline = T, title = "Number of patients"),
             margin = list(t =50),
             xaxis = list(title = "",  showline = T,tickangle = 45)
      ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  })


  ####Variant Analysing selection
  #Updates after gene change 
  
  observeEvent(input$var_gene, {
    
    filtered_cDNA_pos <- all_exchanges.df %>% filter(cDNA_pos == input$search_cDNA_pos, Gene == input$var_gene)
    
    
    updatePickerInput(session = session, inputId = "search_Allele",
                      choices = unique(filtered_cDNA_pos$cDNA_ref))
    
    updatePickerInput(session = session, inputId = "search_cDNA_alt",
                      choices = if(all(is.na(filtered_cDNA_pos$cDNA_alt))){NA}else{unique(filtered_cDNA_pos$cDNA_alt[!is.na(filtered_cDNA_pos$cDNA_alt)])},
                      selected = filtered_cDNA_pos$cDNA_alt[1])
    
    updatePickerInput(session = session, inputId = "search_AA_alt",
                      choices = all_exchanges.df %>% filter(Gene == input$var_gene, AA_pos == filtered_cDNA_pos$AA_pos[1]) %>% .$AA_alt %>% unique(),
                      selected = filtered_cDNA_pos %>% filter(cDNA_alt == cDNA_alt[1]) %>% .$AA_alt)
    
    updatePickerInput(session = session, inputId = "search_AA_ref",
                      choices = filtered_cDNA_pos %>% filter(cDNA_alt == cDNA_alt[1]) %>% .$AA_ref,
                      selected = filtered_cDNA_pos %>% filter(cDNA_alt == cDNA_alt[1]) %>% .$AA_ref)
    
  })
  
  #Updates after DNA change 
  
  observeEvent(input$search_cDNA_pos, {
    
    filtered_cDNA_pos <- all_exchanges.df %>% filter(cDNA_pos == input$search_cDNA_pos, Gene == input$var_gene)
    
    # check that filtered_cDNA_pos is not empty before filtering further
    if (nrow(filtered_cDNA_pos) > 0) {
      filtered_cDNA_pos_aa <- all_exchanges.df %>% filter(AA_pos == filtered_cDNA_pos$AA_pos, Gene == input$var_gene)
      
      updatePickerInput(session = session, inputId = "search_Allele",
                        choices = unique(filtered_cDNA_pos$cDNA_ref))
      
      updatePickerInput(session = session, inputId = "search_cDNA_alt",
                        choices = if(all(is.na(filtered_cDNA_pos$cDNA_alt))){NA}else{unique(filtered_cDNA_pos$cDNA_alt[!is.na(filtered_cDNA_pos$cDNA_alt)])})
      
      updateNumericInputIcon(session = session, inputId = "search_AA_pos",
                             value = if(all(is.na(filtered_cDNA_pos$AA_pos))){NA}else{unique(filtered_cDNA_pos$AA_pos[!is.na(filtered_cDNA_pos$AA_pos)])})
      
      updatePickerInput(session = session, inputId = "search_AA_ref",
                        choices = unique(filtered_cDNA_pos_aa$AA_ref))
      
      updatePickerInput(session = session, inputId = "search_AA_alt",
                        choices = unique(filtered_cDNA_pos_aa$AA_alt))
    }
  })

  
  observeEvent(input$search_cDNA_alt, {
    
    filtered_cDNA_pos <- all_exchanges.df %>% filter(cDNA_pos == input$search_cDNA_pos, Gene == input$var_gene, cDNA_alt == input$search_cDNA_alt)
    
    updatePickerInput(session = session, inputId = "search_AA_alt",
                      selected = unique(filtered_cDNA_pos$AA_alt))
    
  })
  
  #Updates after AA change 
  observeEvent(input$search_AA_alt, {
    filtered_AA_alt <- all_exchanges.df %>% filter(AA_pos == input$search_AA_pos, Gene == input$var_gene, AA_alt == input$search_AA_alt)
    
    
    updatePickerInput(session = session,
                      inputId = c("get_var_type"),
                      choices = c(unique(filtered_AA_alt$Vartype)),
                      selected = c(unique(filtered_AA_alt$Vartype)))
    
  })
  
  observeEvent(input$search_AA_pos, {
    
    filtered_AA_pos <- all_exchanges.df %>% filter(AA_pos == input$search_AA_pos, Gene == input$var_gene)
    
    updatePickerInput(session = session, inputId = "search_AA_ref",
                      choices = sort(unique(filtered_AA_pos$AA_ref)))
    
    updatePickerInput(session = session, inputId = "search_AA_alt",
                      choices = unique(filtered_AA_pos$AA_alt))
    
  })
  
  ## Once the search button is pressed 
  
  varFilterInput <- reactiveValues(data=NULL)

  observeEvent(input$search_var_c, {
    varFilterInput$data <- all_exchanges.df %>% filter(Gene==input$var_gene) %>% filter(cDNA_pos==input$search_cDNA_pos) %>%
      filter(cDNA_ref==input$search_Allele) %>% filter(cDNA_alt==input$search_cDNA_alt)
  })
  
  observeEvent(input$search_var_p, {
    varFilterInput$data <- all_exchanges.df %>% filter(Gene==input$var_gene) %>% filter(AA_pos==input$search_AA_pos) %>%
      filter(AA_ref==input$search_AA_ref) %>% filter(AA_alt==input$search_AA_alt)

  })
  
  # update cDNA input according to user search input , not 100% precise as many option may be possible 
  
  
  observeEvent(input$search_var_p, {
    updateNumericInputIcon(
      session = session,
      inputId = c("search_cDNA_pos"),
      value = c(unique(varFilterInput$data$cDNA_pos)[1])
    )
    
    updatePickerInput(
      session = session,
      inputId = c("search_Allele"),
      choices = c(unique(varFilterInput$data$cDNA_ref)[1])
    )
    
    updatePickerInput(
      session = session,
      inputId = c("search_cDNA_alt"),
      selected = c(unique(varFilterInput$data$cDNA_alt))
    )
    
  })

  #### Variant Information ####
  
  
  output$geneBox1 <- renderValueBox({
    valueBox(
      value = tags$p(paste0("Gene: ",unique(varFilterInput$data$Gene)),style = "font-size:60%"),
      div("Transcript: ",unique(varFilterInput$data$Transcript),br(), br(), "   "), icon = icon("dna"),
      color = "yellow"
    )
  })
  
  output$geneBox2 <- renderValueBox({
    
    if(!is.null(varFilterInput$data)){
      
      
      p_old <- unique(varFilterInput$data$AA_ref)
      p_old_aa <- seqinr::a(p_old)
      p_new <- unique(varFilterInput$data$AA_alt)
      p_new_aa <- ifelse(p_new != "Stop",seqinr::a(p_new),"*")
      p_domain <- Domain_gene.df %>% filter(Gene == varFilterInput$data$Gene[1], AA_pos == varFilterInput$data$AA_pos[1]) %>% .$Domain %>% unique()
      p_pos <- unique(varFilterInput$data$AA_pos)
      
    }else{
      
      p_old <- ""
      p_old_aa <- ""
      p_new <- ""
      p_new_aa <- ""
      p_domain <- ""
      p_pos <- ""
      
    }
    
    valueBox(
      value = tags$p(paste0("Domain: ", p_domain),style = "font-size:60%"),
      div("Amino Acid Position: ",p_pos, br(),
          paste0("Amino Acid Change: ", p_old, " (",p_old_aa, ") "), "-" ,
          paste0(p_new, " (",ifelse(p_new != "Stop",p_new_aa,"*"),")")), icon = icon("dna"),
      color = "light-blue"
    )
  })
  
  output$geneBox3 <- renderValueBox({
    
    if(is.null(varFilterInput$data)){
      
      valueBox(
        value = tags$p(paste0("Control variants: "),style = "font-size:60%"),
        div(paste0("gnomAD Allele Count at this aminoacid position: "),br(),
            paste0("gnomAD Allele Frequency for this aminoacid exchange: ")), icon = icon("dna"),
        color = "green")
      
    }else{
      
      gnomad_count_exchange <- extract_gnomad_features(Control_data.df,varFilterInput$data,"Allele count","exchange")
      
      gnomad_freq_exchange <- extract_gnomad_features(Control_data.df,varFilterInput$data,"Allele freq","exchange")
      
      gnomad_count_position <- extract_gnomad_features(Control_data.df,varFilterInput$data,"Allele count","position")
      
      valueBox(
        value = tags$p(paste0("Control variants: ",  gnomad_count_exchange),style = "font-size:80%"),
        div(paste0("gnomAD Allele Count at this aminoacid position: ", gnomad_count_position),br(),
            paste0("gnomAD Allele Frequency for this aminoacid exchange: ", gnomad_freq_exchange)), icon = icon("dna"),
        color = "green"
      )
    }
  })

  # Table with patient information ####
  output$patientTable <- DT::renderDataTable({
    
    validate(
      need(!plyr::empty(varFilterInput$data),
           "There is no data that matches your filters."))
    
    datatable(Patient_data_missense_only.df %>% filter(Gene==varFilterInput$data$Gene[1]) %>%
                filter(AA_pos==varFilterInput$data$AA_pos[1]) %>%
                filter(AA_alt == varFilterInput$data$AA_alt[1]) %>%
                select(Transcript,Gene, Domain, cDNA, Original_AA_change, phenotype, phenotype_GD,Mean_EnzymeActivity,Refrence) %>%
                mutate(Refrence = sprintf('<a href="%s">%s</a>', Refrence, Refrence)),
              colnames = c("Transcript","Gene", "Domain"," cDNA", "Protein",phenotype_name1, phenotype_name2,"Mean_EnzymeActivity (Article)", "Refrence"),
              options = list(dom = 't', scrollY = TRUE), escape=FALSE)
  })
  
  
  
  output$paraTable <- DT::renderDataTable({
    
    validate(
      need(!plyr::empty(varFilterInput$data),
           "There is no data that matches your filters."))
    
    aln_pos <- all_exchanges.df$Aln_pos[which(all_exchanges.df$AA_pos == varFilterInput$data$AA_pos[1] & all_exchanges.df$Gene ==varFilterInput$data$Gene)]
    
    
    datatable(Patient_data_missense_only.df %>%
                filter(Aln_pos== aln_pos[1],
                       Gene!=varFilterInput$data$Gene[1]) %>%
                select(Transcript,Gene, Domain, cDNA, Protein, phenotype, phenotype_GD,Mean_EnzymeActivity,Refrence) %>%
                mutate(Refrence = sprintf('<a href="%s">%s</a>', Refrence, Refrence)),
              colnames = c("Transcript","Gene", "Domain"," cDNA", "Protein",phenotype_name1, phenotype_name2, "Mean_EnzymeActivity (Article)","Refrence"),
              options = list(dom = 't', scrollY = TRUE), escape=FALSE)
    
  })
  
  
  # table with similar variants according to user input
  
  
  output$compareTable <- DT::renderDataTable({
    
    validate(
      need(!plyr::empty(varFilterInput$data),
           "There is no data that matches your filters."))
    
    if (input$other_genes == F) {
      z <- Patient_data_missense_only.df %>%
        filter(case_when(
          input$compareButtons =="Variant Type" ~ Vartype ==varFilterInput$data$Vartype,
          input$compareButtons =="Domain" ~ Domain==varFilterInput$data$Domain))
      select(Transcript,Gene, Domain, cDNA, Protein,phenotype,phenotype_GD, Mean_EnzymeActivity,Refrence) %>%
        filter(Gene == varFilterInput$data$Gene)
    }
    
    if (input$other_genes == T) {
      
      z <- Patient_data_missense_only.df %>%
        filter(case_when(
          input$compareButtons =="Variant Type" ~ Vartype ==varFilterInput$data$Vartype,
          input$compareButtons =="Domain" ~ Domain==varFilterInput$data$Domain))
      select(Transcript,Gene, Domain, cDNA, Protein, phenotype,phenotype_GD, Mean_EnzymeActivity,Refrence)
    }
    
    
    if (input$hide_unknown == TRUE) {
      z <-
        z %>% filter(
          functional_effect != "Unknown"
        )
    }
    
    
    datatable(z, extensions = "Buttons",
              colnames = c("Transcript","Gene", "Domain"," cDNA", "Protein", phenotype_name1, phenotype_name2, "Mean_EnzymeActivity (Article)","Refrence"),
              options = list(dom = 'Brtip',
                             buttons = c('csv', 'excel'), pageLength=100, scrollY = "350px"), escape = FALSE)
  })
  
  output$comparePlot <- renderPlotly({
    
    validate(need(
      nrow(varFilterInput$data) >0,
      "There is no data that matches your filters."
    ))

    
    if (input$other_genes == F) {
      z <- Patient_data_missense_only.df %>%
        filter(case_when(
          input$compareButtons =="Variant Type" ~ Vartype ==varFilterInput$data$Vartype,
          input$compareButtons =="Domain" ~ Domain==varFilterInput$data$Domain))%>%
        select(Transcript,Gene, Domain, cDNA, Protein, phenotype,phenotype_GD, Mean_EnzymeActivity,Refrence) %>% 
        filter(Gene == varFilterInput$data$Gene)
    }
    
    if (input$other_genes == T) {
      
      z <- Patient_data_missense_only.df %>%
        filter(case_when(
          input$compareButtons =="Variant Type" ~ Vartype ==varFilterInput$data$Vartype,
          input$compareButtons =="Domain" ~ Domain==varFilterInput$data$Domain))
      
    }

    
    plotty1 <- plot_ly(data = z  %>%
                         rename(selected_pheno = "phenotype") %>%
                         add_count(selected_pheno, name = "n_pheno") %>% distinct(selected_pheno,n_pheno) %>%
                         assign("save",.,envir = .GlobalEnv),
                       x = ~ selected_pheno,
                       y = ~ n_pheno,
                       color = ~ selected_pheno,
                       colors = basic_phenotype_colors[1:nrow(save)],
                       type = "bar",
                       hoverinfo = "text", showlegend = FALSE,
                       text= ~ paste0(n_pheno, " individuals")) %>%
      layout(title="",
             font=plotly_font,
             xaxis = list(title="",showline = T, tickangle = 45),
             yaxis = list(title="N of individuals",showline = T),
             margin = list(b = 160)) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    plotty2 <- plot_ly(data = z  %>%
                         rename(selected_pheno = "phenotype_GD") %>%
                         add_count(selected_pheno, name = "n_pheno") %>% distinct(selected_pheno,n_pheno) %>%
                         assign("save",.,envir = .GlobalEnv),
                       x = ~ selected_pheno,
                       y = ~ n_pheno,
                       color = ~ selected_pheno,
                       colors = basic_phenotype_colors[1:nrow(save)],
                       type = "bar",
                       hoverinfo = "text", showlegend = FALSE,
                       text= ~ paste0(n_pheno, " individuals")) %>%
      layout(title="",
             font=plotly_font,
             xaxis = list(title="",showline = T, tickangle = 45),
             yaxis = list(title="N of individuals",showline = T),
             margin = list(b = 160)) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    

    plotty <- subplot(plotty1,plotty2, nrows = 1,titleY = T) %>%
      layout(annotations = list(
        list(x = 0.1 , y = 1.1, text = "Phenotype PD", showarrow = F, xref='paper', yref='paper'),
        list(x = 0.8 , y = 1.1, text = "Phenotype GD", showarrow = F, xref='paper', yref='paper'))
      )
    
    
    return(plotty)
    
    
  })

  #####Research #####

  res_mod <- callModule(
    module = selectizeGroupServer,
    id = "research-filters",
    data = Patient_data.df %>% mutate(ID_x = case_when(ID_x == "No"~ " No",
                                                       TRUE~ID_x)),
    vars = c("Gene","Vartype",  "AA_alt", "Domain", "Phenotype", "phenotype_GD", "ID_x","ID")
  )
  
  output$filtered_n <- renderText({
    x <- nrow(res_mod())
    x <- paste("Number of individuals:", x)
    return(x)
  })
  
  # Table with displayed variants
  output$subsetTable <- DT::renderDataTable({
    req(res_mod())
    
    z <- res_mod() #%>% filter(!is.na(variant_p)) %>% filter(origin != "biparental") %>%

    patient_table <- datatable(z %>% select(Transcript,Gene, Domain, cDNA, Protein,phenotype) , extensions = "Buttons",
                               colnames = c("Transcript","Gene", "Domain"," cDNA", "Protein", phenotype_name1), escape=FALSE,
                               options = list(dom = 'Brtip',buttons = c('csv', 'excel'), pageLength=300, scrollY = "350px"))
    
    
    annotation_table <-  datatable(z %>% select(Transcript,Gene, Domain, cDNA, Protein,phenotype), extensions = "Buttons",
                                   colnames = c("Transcript","Gene", "Domain"," cDNA", "Protein", phenotype_name1), escape=FALSE,
                                   options = list(dom = 'Brtip',buttons = c('csv', 'excel'), pageLength=300, scrollY = "350px"))
    
    if (input$patientFunSwitch == FALSE) {
      return(patient_table)
    } else {
      return(annotation_table)
    }
  })
  
  # 2D lolliplot with GRIN variants
  ### Genotype Interface
  
  output$Genotype_overview_plot <- renderPlotly({
    
    
    # 2D lolliplot variants
    
    research_genotype_domain.df <- all_exchanges.df %>% 
      filter(Gene %in% names(Gene_colors)) %>%
      distinct(Domain,Gene,AA_pos,Domain_color) %>%
      group_by(Gene,Domain,Domain_color) %>% 
      summarise(start = min(AA_pos),
                end = max(AA_pos)) 
    
    filtered_research_genotype_domain_signal.df <- research_genotype_domain.df %>%
      filter(Domain_color == "Signal peptide")
    
    filtered_research_genotype_domain_immune.df <- research_genotype_domain.df %>%
      filter(Domain_color == "immunoglobulin like structure")
    
    filtered_research_genotype_domain_beta.df <- research_genotype_domain.df %>%
      filter(Domain_color == "antiparallel beta sheet structure")
    
    filtered_research_genotype_domain_tim.df <- research_genotype_domain.df %>%
      filter(Domain_color == "TIM structure")
    
    filtered_research_genotype_domain_un.df <- research_genotype_domain.df %>%
      filter(Domain_color == "UNKNOWN")

    g <- ggplot(data=all_exchanges.df %>%
                  filter(Gene %in% names(Gene_colors)) %>%
                  distinct(AA_pos,Gene,Domain,Domain_color) %>%
                  left_join(res_mod() %>%
                              filter(!is.na(AA_pos)) %>%
                              select(Gene,Protein,AA_pos,Vartype) %>%
                              group_by(Protein,AA_pos,Gene,Vartype) %>%
                              summarise(var_count = n(), .groups = "keep") %>%
                              ungroup() %>%
                              mutate(Protein_count = paste0(" ",Protein,", Variant count ",var_count)) %>%
                              group_by(Gene,AA_pos,Vartype) %>%
                              summarise(Protein_final = paste(Protein_count, collapse = ";"),multiple = "all")))+
      
      geom_segment(aes(x=AA_pos, xend=AA_pos,y=4,yend=ifelse(Vartype=="Missense", 7,8)), colour="black")+
      geom_point(aes(x=AA_pos, y=ifelse(Vartype=="Missense", 7,8), color=Vartype, text=Protein_final))+
      # geom_rect(data=research_genotype_domain.df, aes(xmin=start, xmax=end, ymin=3, ymax=4, fill=Domain_color, text=Domain))+
      geom_rect(data=filtered_research_genotype_domain_signal.df, aes(xmin=start, xmax=end, ymin=4, ymax=4.6, fill=Domain_color, text=Domain)) +
      geom_rect(data=filtered_research_genotype_domain_immune.df, aes(xmin=start, xmax=end, ymin=3, ymax=3.6, fill=Domain_color, text=Domain)) +
      geom_rect(data=filtered_research_genotype_domain_beta.df, aes(xmin=start, xmax=end, ymin=2, ymax=2.6, fill=Domain_color, text=Domain)) +
      geom_rect(data=filtered_research_genotype_domain_tim.df, aes(xmin=start, xmax=end, ymin=1, ymax=1.6, fill=Domain_color, text=Domain)) +
      # geom_rect(data=filtered_research_genotype_domain_un.df, aes(xmin=start, xmax=end, ymin=3, ymax=4, fill=Domain_color, text=Domain)) +
      theme_classic()+
      ylim(c(1,10))+
      labs(x= "Mutiple sequence")+
      scale_color_manual(values = lolliplot_fill_scheme)+
      scale_fill_manual(values = lolliplot_fill_scheme)+
      facet_grid(Gene ~ .)+
      theme(
        text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none"
      )
    
    
    if (input$gnomad_m == TRUE) {
      g <- g + geom_point(data=Control_data.df %>% filter(Gene %in% names(Gene_colors)),
                          size=1, aes(x=AA_pos, y=2, alpha=0.1*Allele_count, text=paste0("Position: ",AA_pos,", Allele count: ", Allele_count)))
    }
    
    g <- ggplotly(g, tooltip = "text") %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  %>%
      layout(title="",font=plotly_font,
             xaxis = list(title = "Mutiple sequence")
      )
    
  })
  
  output$Genotype_legend_plot <- renderPlot({
    legend <- data.frame(x=c(1,11,21), y=c(1, 1, 1), text=c("Missense", "PTV", "synonymous"))
    plot <- ggplot(legend, aes(x=x, y=y, color=text))+
      geom_point(size = 6)+
      scale_color_manual(values = c("Missense"="#D55E00","PTV"="#0072B2","synonymous" ="grey"))+
      ylim(c(0,2))+
      xlim(c(0,40))+
      theme_void()+
      geom_text(aes(label=text), hjust=-0.4, color="black")+
      theme(legend.position = "none")

    return(plot)
  })
  

  # # 3D plot
  output$threeDmolGene1 <- renderR3dmol({
    
    print("Structuregene1")
    
    map_var_3d(res_mod(),"GBA1",input$gnomad_m,pdb_sel_gene1,structure_coordinates_gene1,2)
    
    
  })

  output$threeDmolGene2 <- renderR3dmol({
    
    print("Structuregene2")
    
    map_var_3d(res_mod(),"GBA1",input$gnomad_m,pdb_sel_gene2,structure_coordinates_gene2,3)
    
    
  })

  output$research_phenotype1 <- renderPlotly({
    
    gene_sel.df <- res_mod()
    
    if(gene_sel.df%>% unique() %>% length() != 4){
      
      control_sel.df <- Control_data.df %>%
        filter(Gene %in% unique(gene_sel.df$Gene),
               Domain  %in% unique(gene_sel.df$Domain)
        ) %>%
        select(Domain) %>%
        mutate(Gene = "Control")
      # 
      # 
    }else{
      
      control_sel.df <- Control_data.df %>%
        filter(Gene %in% unique(gene_sel.df$Gene),
               Domain  %in% unique(gene_sel.df$Domain)
        ) %>%
        select(Domain) %>%
        mutate(Gene = "Control")
    }
    
    plot <- res_mod() %>%
      filter(!is.na(AA_pos)) %>%
      select(Gene,Domain) %>%
      bind_rows(control_sel.df) %>%
      arrange(Domain) %>%
      mutate(Domain = factor(Domain, levels = unique(Domain))) %>%
      group_by(Gene, Domain) %>%
      summarise(n = n(), .groups = "drop")  %>%
      ungroup() %>%
      distinct() %>%
      plot_ly(
        x=~Domain, y=~(n), split=~Gene, type="bar", color=~Gene, alpha = 0.8, colors= Gene_colors) %>%
      layout(title=research_phenotype1_title,font=plotly_font,
             margin = list(t =50),
             yaxis = list(type = "log",
                          title = "Number of Patient/Controls",
                          ticktext = list("1","5", "10","20","50", "100", "200","500","1000","10000"),
                          tickvals = list(1, 5,10,20,50,100,200,500,1000,10000),
                          showline = T,
                          tickmode = "array"
             ),
             xaxis = list(title = "", showline = T, tickangle = 45)) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    return(plot)
    
  })
  # 
  output$research_phenotype2 <- renderPlotly({
    
    data.df <- res_mod() %>%filter(!is.na(AA_pos)) %>%
      rename(Phenotype1 = phenotype) %>%
      group_by(Gene, Phenotype1) %>% summarise(n = n(), .groups = "drop")  %>% ungroup()
    
    validate(
      # Old Code
      need(nrow(data.df) > 0, 'No data exists, for selection!')
      # With parentheses fix
    )
    Gene_colors1 <- c("#2ca02c")
    
    
    data.df %>%
      plot_ly(
        x=~Phenotype1, y=~n, split=~Gene, type="bar", color=~Gene, alpha = 0.8, colors= Gene_colors1) %>%
      layout(title=research_phenotype2_title ,font=plotly_font, yaxis = list(showline = T, title = "Number of patients"),
             margin = list(t =50),
             xaxis = list(title = "",  showline = T,tickangle = 45)
      ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    
  })
  
  output$research_phenotype3 <- renderPlotly({
    
    data.df <- res_mod() %>%filter(!is.na(AA_pos)) %>%
      rename(Phenotype2 = phenotype_GD) %>%
      group_by(Gene, Phenotype2) %>% summarise(n = n(), .groups = "drop")  %>% ungroup()
    
    validate(
      # Old Code
      need(nrow(data.df) > 0, 'No data exists, for selection!')
      # With parentheses fix
    )
    Gene_colors2 <- c("#9467bd")
    
    
    data.df %>%
      plot_ly(
        x=~Phenotype2, y=~n, split=~Gene, type="bar", color=~Gene, alpha = 0.8, colors= Gene_colors2) %>%
      layout(title="Phenotype GD" ,font=plotly_font, yaxis = list(showline = T, title = "Number of patients"),
             margin = list(t =50),
             xaxis = list(title = "",  showline = T,tickangle = 45)
      ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  })
  
  output$research_phenotype4 <- renderPlotly({
    
    data.df <- res_mod() %>% filter(!is.na(AA_pos),!is.na(ID_x)) %>%
      rename(Phenotype3 = ID_x) %>%
      filter(Phenotype3 != "Unknown") %>%
      mutate(Phenotype3 = ifelse(Phenotype3 == "Yes"," Yes",Phenotype3)) %>%
      group_by(Gene, Phenotype3) %>% summarise(n = n(), .groups = "drop")  %>% ungroup()
    
    validate(
      # Old Code
      need(nrow(data.df) > 0, 'No data exists, for selection!')
      # With parentheses fix
    )
    Gene_colors2 <- c("#ff7f0e")
    
    data.df %>%
      plot_ly(
        x=~Phenotype3, y=~n, split=~Gene, type="bar", color=~Gene, alpha = 0.8, colors= Gene_colors2) %>%
      layout(title="Severity of PD" ,font=plotly_font, yaxis = list(showline = T, title = "Number of patients"),
             margin = list(t =50),
             xaxis = list(title = "",  showline = T,tickangle = 45)
      ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
  })
  output$research_phenotype5 <- renderPlotly({
    
    data.df <- res_mod() %>% filter(!is.na(AA_pos),!is.na(ID)) %>%
      rename(Phenotype4 = ID) %>%
      filter(Phenotype4 != "Unknown") %>%
      mutate(Phenotype4 = ifelse(Phenotype4 == "Yes"," Yes",Phenotype4)) %>%
      group_by(Gene, Phenotype4) %>% summarise(n = n(), .groups = "drop")  %>% ungroup()
    
    validate(
      # Old Code
      need(nrow(data.df) > 0, 'No data exists, for selection!')
      # With parentheses fix
    )
    Gene_colors2 <- c("#2ca02c","#ff7f0e")
    
    data.df %>%
      plot_ly(
        x=~Phenotype4, y=~n, split=~Gene, type="bar", color=~Gene, alpha = 0.8, colors= Gene_colors2) %>%
      layout(title="Severity of GD" ,font=plotly_font, yaxis = list(showline = T, title = "Number of patients"),
             margin = list(t =50),
             xaxis = list(title = "",  showline = T,tickangle = 45)
      ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    # phenotype_research(data.df,"Developmental delay GD ")
    
  })
  
  output$research_phenotype6 <- renderPlotly({
    validate(
      need(!is.null(paraz_mtr_input.df), 'No data exists for selection!')
    )
    
    data.df <- paraz_mtr_input.df %>%
      arrange(`Mean `)
    
    # Filter out the trace you want to remove based on a condition
    filtered_data <- data.df %>%
      filter(`Mean ` %in% paraz_mtr_input.df$`Mean `)
    
    # Convert 'Mean' column to numeric
    filtered_data$`Mean ` <- as.numeric(as.character(filtered_data$`Mean `))
    
    # Calculate max and median values
    max_value <- max(filtered_data$`Mean `, na.rm = TRUE)
    median_value <- median(filtered_data$`Mean `, na.rm = TRUE)
    
    # Create box plot
    p <- plot_ly(filtered_data, y = ~`Mean `, x = "Enzyme activity", type = "box",
                 alpha = 0.8, colors = Gene_colors, name = "Box Plot") %>%
      # Add scatter plot for Variants and Mean values
      add_trace(y = ~`Mean `, x ="Enzyme activity", type = "scatter",
                mode = "markers", marker = list(size = 6, color = "#000000"),
                # text = ~`Variants`, hoverinfo = "text", name = "Variants") %>%
                text = ~paste("Variants: ", `Variants`, "<br>Enzyme activity: ", `Mean `),hoverinfo = "text", name = "Variants") %>%
      layout(
        title = "Glucocerebrosidase activity of GBA1",
        font = plotly_font,
        margin = list(t = 50),
        xaxis = list(
          title = "GBA1",
          tickangle = 45,
          showline = TRUE,
          tickvals = filtered_data$Variants,  # Update tick values
          ticktext = filtered_data$Variants  # Update tick labels
        ),
        yaxis = list(
          title = "β-glucocerebrosidase activity (μmol/L/h)",
          showline = TRUE
        )
      ) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  })

  output$research_phenotype7 <- renderPlotly({
    validate(
      need(!is.null(paraz_mtr_input.df), 'No data exists for selection!')
    )
    merged_data <- left_join(paraz_mtr_input.df, paraz_mtr_input_control.df, by = "vartype", multiple = "all")
    
    # Convert 'Mean' and 'Control' columns to numeric
    merged_data$`Mean ` <- as.numeric(as.character(merged_data$`Mean `))
    merged_data$Control <- as.numeric(as.character(merged_data$Control))
    
    # Remove non-finite values
    filtered_data <- na.omit(merged_data)
    
    # Calculate the mean of the 'Control' column
    control_mean <- mean(filtered_data$Control, na.rm = TRUE)
    boxplot_colors <- c("Unclassified" = "blue", "Risk" = "red", "Mild" = "green", "Severe" = "purple")
    
    # Define a vector of colors
    p <- ggplot(filtered_data, aes(x = factor(vartype, levels = c("Unclassified", "Risk", "Mild", "Severe")), y = `Mean `, fill = factor(vartype))) +
      geom_boxplot(alpha = 0.8, color = "black") +
      geom_hline(yintercept = control_mean, linetype = "dashed", color = "black") +  # Add horizontal line for the Control value
      scale_fill_manual(values = boxplot_colors) +
      labs(title = "β-glucocerebrosidase activity and mutation type",
           x = "GBA",
           y = "β-glucocerebrosidase activity (μmol/L/h)",
           fill = "Mutation Type") +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black")
      )
    
    plotly::plotly_build(p) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  })

  # ##### Compare variants ####
  # 
  res_mod1 <- callModule(
    module = selectizeGroupServer,
    id = "research-group1",
    data = Patient_data.df %>% mutate(ID_x = case_when(ID_x == "No"~ " No",
                                                       TRUE~ID_x)),
    vars = c("Gene","Vartype", "Domain", "Phenotype","ID_x")
  )
  
  output$filtered_n_compare1 <- renderText({
    x <- nrow(res_mod1())
    x <- paste("Number of individuals:",x)
    return(x)
  })
  # 
  # 
  res_mod2 <- callModule(
    module = selectizeGroupServer,
    id = "research-group2",
    data = Patient_data.df %>% mutate(ID_x = case_when(ID_x == "No"~ " No",
                                                       TRUE~ID_x)),
    vars = c("Gene","Vartype", "Domain", "Phenotype","ID_x")
  )
  # 
  output$filtered_n_compare2 <- renderText({
    x <- nrow(res_mod2())
    x <- paste("Number of individuals:", x)
    return(x)
  })
  res_mod3 <- callModule(
    module = selectizeGroupServer,
    id = "research-group1",
    data = Patient_data.df %>% mutate(ID = case_when(ID == "No"~ " No",
                                                     TRUE~ID)),
    vars = c("Gene","Vartype", "Domain","phenotype_GD", "ID")
  )
  
  output$filtered_n_compare3 <- renderText({
    x <- nrow(res_mod3())
    x <- paste("Number of individuals:",x)
    return(x)
  })
  res_mod4 <- callModule(
    module = selectizeGroupServer,
    id = "research-group2",
    data = Patient_data.df %>% mutate(ID = case_when(ID == "No"~ " No",
                                                     TRUE~ID)),
    vars = c("Gene","Vartype", "Domain", "phenotype_GD","ID")
  )
  
  output$filtered_n_compare4 <- renderText({
    x <- nrow(res_mod4())
    x <- paste("Number of individuals:",x)
    return(x)
  })
  # 
  ##overview plot
  output$Genotype_legend_plot_compare <- renderPlot({
    
    colorGroup1 <- color_list_research_compare[input$colorGroup1] %>% as.vector()
    colorGroup2 <- color_list_research_compare[input$colorGroup2] %>% as.vector()
    
    legend <- data.frame(x=c(1,11), y=c(1, 1), text=c("Group 1", "Group 2"))
    plot <- ggplot(legend, aes(x=x, y=y, color=text))+
      geom_point(size = 6)+
      scale_color_manual(values = c("Group 1"= colorGroup1,"Group 2"= colorGroup2))+
      ylim(c(0,2))+
      xlim(c(0,40))+
      theme_void()+
      geom_text(aes(label=text), hjust=-0.4, color="black")+
      theme(legend.position = "none")
    
    
    
    return(plot)
    
  })
  # 
  # 
  output$Genotype_overview_plot_compare <- renderPlotly({
    
    colorGroup1 <- color_list_research_compare[input$colorGroup1] %>% as.vector()
    colorGroup2 <- color_list_research_compare[input$colorGroup2] %>% as.vector()
    
    
    data_pre_sel1.df <- res_mod1() %>%
      select(Gene,Phenotype,Domain,AA_alt,Vartype,AA_pos,ID_x)
    
    data_pre_sel2.df <- res_mod2() %>%
      select(Gene,Phenotype,Domain,AA_alt,Vartype,AA_pos,ID_x)
    
    data.df <- rbind(Patient_data.df %>%
                       #rename(Phenotype1 = Phenotype) %>%
                       filter(!is.na(Phenotype)) %>%
                       filter(Gene %in% data_pre_sel1.df$Gene,
                              Domain %in% data_pre_sel1.df$Domain,
                              Phenotype %in% data_pre_sel1.df$Phenotype,
                              AA_alt %in% data_pre_sel1.df$AA_alt,
                              Vartype %in% data_pre_sel1.df$Vartype,
                              ID_x %in% data_pre_sel1.df$ID_x) %>%
                       # functional_effect %in% data_pre_sel1.df$functional_effect) %>%
                       mutate(group = "Group 1"),
                     Patient_data.df %>%
                       #rename(Phenotype1 = Phenotype) %>%
                       filter(!is.na(Phenotype)) %>%
                       filter(Gene %in% data_pre_sel2.df$Gene,
                              Domain %in% data_pre_sel2.df$Domain,
                              Phenotype %in% data_pre_sel2.df$Phenotype,
                              AA_alt %in% data_pre_sel2.df$AA_alt,
                              Vartype %in% data_pre_sel2.df$Vartype,
                              ID_x %in% data_pre_sel2.df$ID_x) %>%
                       # functional_effect %in% data_pre_sel2.df$functional_effect) %>%
                       mutate(group = "Group 2"))
    
    lolliplot_fill_scheme <- c(lolliplot_fill_scheme,"Group 1" = colorGroup1,"Group 2" = colorGroup2)
    
    research_genotype_domain.df <- all_exchanges.df %>%
      filter(Gene %in% names(Gene_colors)) %>%
      distinct(Domain,Gene,AA_pos,Domain_color) %>%
      group_by(Gene,Domain,Domain_color) %>%
      summarise(start = min(AA_pos),
                end = max(AA_pos))
    
    g <- ggplot(data=all_exchanges.df %>%
                  filter(Gene %in% names(Gene_colors)) %>%
                  distinct(AA_pos,Gene,Domain,Domain_color) %>%
                  left_join(data.df  %>%
                              filter(!is.na(AA_pos)) %>%
                              group_by(Protein,AA_pos,Gene,group) %>%
                              summarise(var_count = n(), .groups = 'keep') %>% 
                              ungroup() %>%
                              mutate(Protein_count = paste0(" ",Protein,", Variant count ",var_count)) %>%
                              group_by(Gene,AA_pos,group) %>%
                              reframe(Protein_final = paste(Protein_count, collapse = ";"), multiple = "all"),
                            by = c('AA_pos', 'Gene'),multiple = "all"))+
      geom_segment(aes(x=AA_pos, xend=AA_pos, y=4, yend=ifelse(group=="Group 1", 7,8)), colour="black")+
      geom_point(aes(x=AA_pos, y=ifelse(group=="Group 1", 7,8), color=group, text = Protein_final))+
      geom_rect(data=research_genotype_domain.df, aes(xmin=start, xmax=end, ymin=3, ymax=4, fill=Domain_color, text=Domain))+
      theme_classic()+
      ylim(c(1,10))+
      labs( x= "Mutiple sequence")+
      scale_color_manual(values = lolliplot_fill_scheme)+
      scale_fill_manual(values = lolliplot_fill_scheme)+
      facet_grid(Gene ~ .)+
      theme(
        text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none"
      )
    
    
    if (input$gnomad_m == TRUE) {
      g <- g + geom_point(data=Control_data.df %>% filter(Gene %in% names(Gene_colors)),
                          size=1, aes(x=AA_pos, y=2, alpha=0.1*Allele_count, text=paste0("Position: ",AA_pos,", Allele count: ", Allele_count)))
    }

    g <- ggplotly(g, tooltip = "text") %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  %>%
      layout(title="",font=plotly_font,
             xaxis = list(title = "Mutiple sequence")
      )
    
  })
  # 
  #display structure
  output$structure_legend_plot_compare <- renderPlot({
    
    colorGroup1 <- color_list_research_compare[input$colorGroup1] %>% as.vector()
    colorGroup2 <- color_list_research_compare[input$colorGroup2] %>% as.vector()
    
    legend <- data.frame(x=c(1,18,27), y=c(1, 1, 1), text=c("Variants in both groups", "Group 1", "Group 2"))
    plot <- ggplot(legend, aes(x=x, y=y, color=text))+
      geom_point(size = 6)+
      scale_color_manual(values = c("Variants in both groups"="gray","Group 1" = colorGroup1,"Group 2" = colorGroup2))+
      ylim(c(0,2))+
      xlim(c(-5,38))+
      theme_void()+
      geom_text(aes(label=text, x = c(-2,18,27)), hjust=-0.4, color="black")+
      theme(legend.position = "none")
    
    return(plot)
  })
  # 
  
  output$threeDmol_compare <- renderR3dmol({
    
    data_pre_sel1.df <- res_mod1() %>%
      select(Gene,Phenotype,Domain,AA_alt,Vartype,AA_pos,ID_x)
    
    data_pre_sel2.df <- res_mod2() %>% #filter(Gene == "SCN2A") %>%
      select(Gene,Phenotype,Domain,AA_alt,Vartype,AA_pos,ID_x)
    
    data.df <- rbind(Patient_data.df %>%
                       #rename(Phenotype1 = Phenotype) %>%
                       filter(!is.na(Phenotype)) %>%
                       filter(Gene %in% data_pre_sel1.df$Gene,
                              Domain %in% data_pre_sel1.df$Domain,
                              Phenotype %in% data_pre_sel1.df$Phenotype,
                              AA_alt %in% data_pre_sel1.df$AA_alt,
                              Vartype %in% data_pre_sel1.df$Vartype,
                              ID_x %in% data_pre_sel1.df$ID_x) %>%
                       # functional_effect %in% data_pre_sel1.df$functional_effect
                       mutate(group = "Group 1"),
                     Patient_data.df %>%
                       #rename(Phenotype1 = Phenotype) %>%
                       filter(!is.na(Phenotype)) %>%
                       filter(Gene %in% data_pre_sel2.df$Gene,
                              Domain %in% data_pre_sel2.df$Domain,
                              Phenotype %in% data_pre_sel2.df$Phenotype,
                              AA_alt %in% data_pre_sel2.df$AA_alt,
                              Vartype %in% data_pre_sel2.df$Vartype,
                              ID_x %in% data_pre_sel2.df$ID_x) %>%
                       # functional_effect %in% data_pre_sel2.df$functional_effect) %>%
                       mutate(group = "Group 2"))
    
    
    variant.df <- data.df %>%
      filter(Vartype == "Missense") %>%
      mutate(label = "pathogenic") %>%
      group_by(AA_pos,AA_ref,Gene,group) %>%
      summarise(n_occ = n()) %>%
      select(AA_pos,AA_ref,n_occ,Gene,group)
    
    gnomad.df <- Control_data.df %>%
      filter(Gene %in% variant.df$Gene) %>%
      group_by(AA_pos,AA_ref,Gene,group) %>%
      summarise(n_occ = n()) %>%
      select(AA_pos,AA_ref,n_occ,Gene,group)
    
    structure.df <- read_delim("data/pdb/6tjk.txt",delim = "\t") %>%
      mutate(Aminoacid = aaa(Aminoacid)) %>%
      select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
    
    variant.df$AA_ref <- unlist(variant.df$AA_ref)
    
    # For variant.df
    variant.df <- variant.df %>%
      left_join(structure.df %>% unnest(Aminoacid), by = c("AA_pos" = "Uniprot_position", "AA_ref" = "Aminoacid", "Gene" = "gene"), suffix = c(".variant", ".structure"), multiple = "all") %>%
      mutate(struc_cov = ifelse(is.na(Position_in_structure), "no", "yes")) %>% 
      filter(struc_cov == "yes") %>%
      distinct(Position_in_structure, group.variant, Gene.variant, .keep_all = TRUE) %>%
      group_by(Position_in_structure, Gene.variant) %>%
      summarise(var_mut = ifelse(n() > 1, "multiple", group.variant))
    
    # For gnomad.df
    gnomad.df <- gnomad.df %>%
      full_join(structure.df %>% unnest(Aminoacid), by = c("AA_pos" = "Uniprot_position", "AA_ref" = "Aminoacid", "Gene" = "gene"), suffix = c(".gnomad", ".structure"), multiple = "all") %>%
      mutate(struc_cov = ifelse(is.na(Position_in_structure), "no", "yes")) %>% 
      filter(struc_cov == "yes") %>%
      distinct(Position_in_structure, Gene.gnomad, .keep_all = TRUE) %>%
      group_by(Position_in_structure, Gene.gnomad) %>%
      summarise(var_mut = ifelse(n() > 1, "multiple", Gene.gnomad))
    
    
    
    
    colorGroup1 <- color_list_research_compare[input$colorGroup1] %>% as.vector()
    colorGroup2 <- color_list_research_compare[input$colorGroup2] %>% as.vector()
    
    
    sub_color <- c("#e2f970","#6fbbf7","#ee6c71","#ffbc5a","#bf73cc")
    sub_scale <- c(1.2,0.8)
    struc_color <- "white"
    
    rot = 270
    rot_axis = "x"
    spin_axis = "vy"
    
    #Specify yourself- color of the cartoon per subunit
    subunit_color <- c("wheat","white") #first color for GRIN1 second or GRIN2A
    
    #Model for the protein complex
    
    modelo <- r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 50,
        upperZoomLimit = 1000
      )
    )
    
    modelo <- modelo %>% m_add_model(data = "data/pdb/GBA.pdb1", format = "pdb")
    
    # Zoom to encompass the whole scene
    modelo <-modelo %>% m_zoom_to() %>%
      # Set color o cartoon representation
      m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
      # Set subunit colors
      m_set_style(
        sel = m_sel(chain = c("A")),
        style = m_style_cartoon(color = subunit_color[2])
      ) %>%
      # visualize variants all
      m_set_style(
        sel = m_sel(resi = variant.df %>% filter(var_mut == "mutiple") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = "gray",
                               scale = sub_scale[1])
      ) %>%
      #visualize variants 
      m_set_style(
        sel = m_sel(resi = variant.df %>% filter(var_mut == "Group 1") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = colorGroup1,
                               scale = sub_scale[1])
      ) %>%
      #visualize variants 
      m_set_style(
        sel = m_sel(resi = variant.df %>% filter(var_mut == "Group 2") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = colorGroup2,
                               scale = sub_scale[1])
      )
    
    
    
    if  (input$gnomad_m == TRUE) {
      
      modelo <- modelo %>% m_set_style(
        sel = m_sel(resi = gnomad.df$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = "#333333",
                               scale = sub_scale[2]))
      
    }
    return(modelo)
    
  })
  
  output$research_compare_pheno1 <- renderPlotly({
    
    data_pre_sel1.df <- res_mod1() %>%  rename(Phenotype1 = Phenotype) %>%
      select(Gene,Phenotype1,Domain,AA_alt,Vartype,ID_x)
    
    data_pre_sel2.df <- res_mod2() %>%  rename(Phenotype1 = Phenotype) %>%
      select(Gene,Phenotype1,Domain,AA_alt,Vartype,ID_x)
    
    data.df <- rbind(Patient_data.df %>%
                       rename(Phenotype1 = Phenotype) %>%
                       filter(!is.na(Phenotype1)) %>%
                       filter(Gene %in% data_pre_sel1.df$Gene,
                              Domain %in% data_pre_sel1.df$Domain,
                              Phenotype1 %in% data_pre_sel1.df$Phenotype1,
                              AA_alt %in% data_pre_sel1.df$AA_alt,
                              Vartype %in% data_pre_sel1.df$Vartype,
                              ID_x %in% data_pre_sel1.df$ID_x) %>%
                       # functional_effect %in% data_pre_sel1.df$functional_effect
                       mutate(group = "Group 1"),
                     Patient_data.df %>%
                       rename(Phenotype1 = Phenotype) %>%
                       filter(!is.na(Phenotype1)) %>%
                       filter(Gene %in% data_pre_sel2.df$Gene,
                              Domain %in% data_pre_sel2.df$Domain,
                              Phenotype1 %in% data_pre_sel2.df$Phenotype1,
                              AA_alt %in% data_pre_sel2.df$AA_alt,
                              Vartype %in% data_pre_sel2.df$Vartype,
                              ID_x %in% data_pre_sel2.df$ID_x) %>%
                       # functional_effect %in% data_pre_sel2.df$functional_effect
                       mutate(group = "Group 2")) %>%
      group_by(Gene, Phenotype1,group) %>%
      summarise(n = n())  %>%
      ungroup()
    
    validate(
      # Old Code
      need(nrow(data.df) > 0, 'No data exists, for selection!')
      # With parentheses fix
    )
    
    Gene_colors <- c(color_list_research_compare[input$colorGroup1] %>% as.vector(),color_list_research_compare[input$colorGroup2] %>% as.vector())
    
    data.df %>%
      plot_ly(
        x=~Phenotype1, y=~n, type="bar", color=~group, alpha = 0.8, colors= Gene_colors) %>%
      layout(title=research_phenotype2_title ,font=plotly_font, yaxis = list(showline = T, title = "Number of patients"),
             margin = list(t =50),
             xaxis = list(title = "",  showline = T,tickangle = 45)
      ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    
  })
  
  output$research_compare_pheno2 <- renderPlotly({
    
    data_pre_sel1.df <- res_mod3() %>%  rename(Phenotype2 = phenotype_GD) %>%
      select(Gene,Phenotype2,Domain,AA_alt,Vartype,ID)
    
    data_pre_sel2.df <- res_mod4() %>%  rename(Phenotype2 = phenotype_GD) %>%
      select(Gene,Phenotype2,Domain,AA_alt,Vartype,ID)
    
    data.df <- rbind(Patient_data.df %>%
                       rename(Phenotype2 = phenotype_GD) %>%
                       filter(!is.na(Phenotype2)) %>%
                       filter(Gene %in% data_pre_sel1.df$Gene,
                              Domain %in% data_pre_sel1.df$Domain,
                              Phenotype2 %in% data_pre_sel1.df$Phenotype2,
                              AA_alt %in% data_pre_sel1.df$AA_alt,
                              Vartype %in% data_pre_sel1.df$Vartype,
                              ID %in% data_pre_sel1.df$ID) %>%
                       # functional_effect %in% data_pre_sel1.df$functional_effect
                       mutate(group = "Group 1"),
                     Patient_data.df %>%
                       rename(Phenotype2 = phenotype_GD) %>%
                       filter(!is.na(Phenotype2)) %>%
                       filter(Gene %in% data_pre_sel2.df$Gene,
                              Domain %in% data_pre_sel2.df$Domain,
                              Phenotype2 %in% data_pre_sel2.df$Phenotype2,
                              AA_alt %in% data_pre_sel2.df$AA_alt,
                              Vartype %in% data_pre_sel2.df$Vartype,
                              ID %in% data_pre_sel2.df$ID) %>%
                       # functional_effect %in% data_pre_sel2.df$functional_effect
                       mutate(group = "Group 2")) %>%
      group_by(Gene, Phenotype2,group) %>%
      summarise(n = n())  %>%
      ungroup()
    
    validate(
      # Old Code
      need(nrow(data.df) > 0, 'No data exists, for selection!')
      # With parentheses fix
    )
    
    Gene_colors <- c(color_list_research_compare[input$colorGroup1] %>% as.vector(),color_list_research_compare[input$colorGroup2] %>% as.vector())
    
    data.df %>%
      plot_ly(
        x=~Phenotype2, y=~n, type="bar", color=~group, alpha = 0.8, colors= Gene_colors) %>%
      layout(title="Phenotype GD" ,font=plotly_font, yaxis = list(showline = T, title = "Number of patients"),
             margin = list(t =50),
             xaxis = list(title = "",  showline = T,tickangle = 45)
      ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    
  })
  
  
  output$research_compare_pheno3 <- renderPlotly({
    
    data.df <- rbind(res_mod1(),res_mod2()) %>% filter(!is.na(ID_x))
    
    validate(
      # Old Code
      need(nrow(data.df) > 0, 'No data exists, for selection!')
      # With parentheses fix
    )
    
    phenotype_compare_research(res_mod1(),res_mod2(),"ID_x","Severity of PD",color_list_research_compare[input$colorGroup2] %>% as.vector(),color_list_research_compare[input$colorGroup1] %>% as.vector())
    
    
  })
  
  output$research_compare_pheno4 <- renderPlotly({
    
    data.df <- rbind(res_mod3(),res_mod4()) %>% filter(!is.na(ID))
    
    validate(
      # Old Code
      need(nrow(data.df) > 0, 'No data exists, for selection!')
      # With parentheses fix
    )
    
    phenotype_compare_research1(res_mod3(),res_mod4(),"ID","Severity of GD",color_list_research_compare[input$colorGroup2] %>% as.vector(),color_list_research_compare[input$colorGroup1] %>% as.vector())
    
    
  })

}

