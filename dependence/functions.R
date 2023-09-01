
############## FUNCTIONS, STYLE ############
# takes a string as input and extracts any numeric values from it.
numextract <- function(string) {
  as.numeric(str_extract(string, "\\-*\\d+\\.*\\d*"))
}

# takes an amino acid letter code as input and converts it to a full name or a special code if the input is "*",
#which represents a stop codon in the genetic code.
convert_aa <- function(aa){
  
  return(ifelse(aa == "*","Stop",aaa(aa)))
  
}

# takes a numeric input and performs a log10 transformation on it. The function first checks whether the input is greater than 0.
#If it is, the function takes the log10 of the input and assigns the result to the variable "x".
log_fun <- function(x) { 
  x = case_when(
    as.numeric(x) > 0 ~ log10(x), #log10 transformation to make the visuals clean 
    as.numeric(x) == 0 ~ 0,
    as.numeric(x) < 0 ~ log10(x * -1) * -1, 
    TRUE ~ -9999)
  return(x)
  
}

## Modify the CSS style of a given selector
modifyStyle <- function(selector, ...) {
  
  values <- as.list(substitute(list(...)))[-1L]
  parameters <- names(values)
  
  args <- Map(function(p, v) paste0("'", p,"': '", v,"'"), parameters, values)
  jsc <- paste0("$('",selector,"').css({", paste(args, collapse = ", "),"});")
  
  shinyjs::runjs(code = jsc)
  
}

plotly_font <- list(
  family = "sans-serif",
  size = 15)

goodbye <- c("zoomIn2d", "zoomOut2d", "hoverClosestGeo", "hoverClosestGl2d",
             "hoverClosestCartesian","lasso2d","select2d","resetScale2d",
             "hoverCompareCartesian", "hoverClosestPie","toggleSpikelines")

# coloring for colorblindness 
lolliplot_fill_scheme <-
  c("Signal peptide" = "red",
    "immunoglobulin like structure" = "yellow",
    "antiparallel beta sheet structure" = "#6666FF",
    "TIM structure" = "green",
    "no" = "#FFFFFF",
    "UNKNOWN" = "darkred",
    "Missense" = "#D55E00",
    "Synonymous" = "#D55E00",
    "PTV" = "#0072B2",
    "Other" = "#0072B2"
  )


######Functions######

####Basic Information 
Phenotype_fac_1.fun <- function(select_gene,phenotype,color_sel){
  
  plot_ly(Patient_data.df %>% 
            filter(Gene == select_gene) %>% 
            rename(phenotype_sel = phenotype) %>% 
            select(phenotype_sel) %>%
            arrange(phenotype_sel) %>% 
            mutate(phenotype_sel = factor(phenotype_sel, levels = unique(phenotype_sel))) %>% 
            group_by(phenotype_sel) %>% 
            summarise(n = n()) %>% 
            assign("save",.,envir = .GlobalEnv), 
          values=~n, labels=~phenotype_sel, type="pie",
          textposition="inside", 
          textinfo="label+percent",
          insidetextfont = list(color = '#333333'),
          hoverinfo = "text",
          text= ~ paste(n, "individuals"),
          marker=list(colors= basic_phenotype_colors[1:nrow(save)],
                      line = list(color = '#FFFFFF', width = 1)), showlegend = FALSE) %>% 
    layout(title=paste(""),font=plotly_font) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  
}

##GBA1, only missense 
##########################
Phenotype_fac_2.fun <- function(select_gene,phenotype,color_sel){
  
  plot_ly(Patient_data.df %>% 
            filter(Gene == select_gene) %>% 
            rename(phenotype_fac = phenotype) %>% 
            filter(Vartype %in% c("Missense","PTV"),
                   !is.na(phenotype),
                   phenotype_fac != "Unknown")%>% 
            select(phenotype_fac,Vartype) %>%
            arrange(phenotype_fac) %>% 
            mutate(phenotype_sel = factor(phenotype_fac, levels = unique(phenotype_fac))) %>% 
            group_by(phenotype_fac,Vartype) %>% 
            summarise(n = n()) %>% 
            mutate(n_gene = sum(n))%>% 
            assign("save",.,envir = .GlobalEnv),
          x = ~ phenotype_fac, 
          y = ~ round(n, digits = 2), 
          color = ~Vartype,
          colors = c("#BEBADA"),
          type = "bar", 
          hoverinfo = "text", showlegend = FALSE,
          text= ~ paste0(round(n, digits = 2), " (" ,n," individuals)")) %>% 
    layout(title = "",
           font=plotly_font,
           xaxis = list(title="",showline = T, tickangle = 45),
           yaxis = list(title="N of individuals",showline = T),
           margin = list(b = 80, t = 50)) %>%
    config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  
}


Phenotype_fac_3.fun <- function(select_gene,phenotype){
  
  plot_ly(Patient_data.df %>% 
            filter(Gene == select_gene) %>% 
            rename(phenotype_fac = phenotype) %>% 
            filter(Vartype %in% c("Missense", "PTV"),
                   !is.na(phenotype),
                   phenotype_fac != "Unknown")%>% 
            select(phenotype_fac,Vartype) %>%
            arrange(phenotype_fac) %>% 
            mutate(phenotype_sel = factor(phenotype_fac, levels = unique(phenotype_fac))) %>% 
            group_by(phenotype_fac,Vartype) %>% 
            summarise(n = n()) %>% 
            mutate(n_gene = sum(n))%>% 
            assign("save",.,envir = .GlobalEnv),
          x = ~ phenotype_fac, 
          y = ~ round(n, digits = 2), 
          color = ~Vartype,
          colors = c("#BEBADA","#FDB462"),
          type = "bar", 
          hoverinfo = "text", showlegend = FALSE,
          text= ~ paste0(round(n, digits = 2), " (" ,n," individuals)")) %>% 
    layout(title = "",
           font=plotly_font,
           xaxis = list(title="",showline = T, tickangle = 45),
           yaxis = list(title="N of individuals",showline = T),
           margin = list(b = 80, t = 50)) %>%
    config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  
}

phenotype_research <- function(data.df,title_sel){
  
  plot <- data.df %>%  
    plot_ly(
      x=~Phenotype2, y=~n, split=~Gene, type="bar", color=~Gene, alpha = 0.8, colors= Gene_colors) %>%
    layout(title= title_sel ,font=plotly_font, yaxis = list(showline = T, title = "Number of patients"),
           margin = list(t =50),
           xaxis = list(title = "",  showline = T,tickangle = 45)
    ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  
  return(plot)
  
}


phenotype_compare_research <- function(group1.df,group2.df,var,title_sel, color1,color2){
  
  data_pre_sel1.df <- group1.df %>% 
    select(Gene,phenotype,Domain,AA_alt,Vartype,ID_x)
  
  data_pre_sel2.df <- group2.df %>%  
    select(Gene,phenotype,Domain,AA_alt,Vartype,ID_x)
  
  
  data.df <- rbind(Patient_data.df %>% 
                     filter(!is.na(Phenotype)) %>% 
                     filter(Gene %in% data_pre_sel1.df$Gene,
                            Domain %in% data_pre_sel1.df$Domain,
                            phenotype %in% data_pre_sel1.df$phenotype,
                            AA_alt %in% data_pre_sel1.df$AA_alt,
                            Vartype %in% data_pre_sel1.df$Vartype,
                            ID_x %in% data_pre_sel1.df$ID_x) %>% 
                     mutate(group = "Group 1"),
                   Patient_data.df %>% 
                     filter(!is.na(phenotype)) %>% 
                     filter(Gene %in% data_pre_sel2.df$Gene,
                            Domain %in% data_pre_sel2.df$Domain,
                            phenotype %in% data_pre_sel2.df$phenotype,
                            AA_alt %in% data_pre_sel2.df$AA_alt,
                            Vartype %in% data_pre_sel2.df$Vartype,
                            ID_x %in% data_pre_sel2.df$ID_x) %>% 
                     mutate(group = "Group 2")) %>% 
    rename(Phenotype2 = var) %>% 
    filter(!is.na(AA_pos),!is.na(Phenotype2)) %>% 
    filter(Phenotype2 != "Unknown") %>% 
    mutate(Phenotype2 = ifelse(Phenotype2 == "Yes"," Yes",Phenotype2)) %>% 
    group_by(Phenotype2,group) %>% summarise(n = n())  %>% ungroup()
  
  if(any(data.df$group == "Group 1") & any(data.df$group == "Group 2")){
    
    Group_colors <- c(color1,color2)
    
  }else if(any(data.df$group == "Group 1")){
    
    Group_colors <- c(color1)
    
  }else{
    Group_colors <- c(color2)
  }
  
  plot <- data.df %>%
    plot_ly(
      x=~Phenotype2, y=~n, type="bar", color=~group, alpha = 0.8, colors= Group_colors) %>%
    layout(title=title_sel ,font=plotly_font, yaxis = list(showline = T, title = "Number of patients"),
           margin = list(t =50),
           xaxis = list(title = "",  showline = T,tickangle = 45)
    ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  
  return(plot)
}

phenotype_compare_research1 <- function(group1.df,group2.df,var,title_sel, color1,color2){
  
  data_pre_sel1.df <- group1.df %>% 
    select(Gene,phenotype_GD,Domain,AA_alt,Vartype,ID)
  
  data_pre_sel2.df <- group2.df %>%  
    select(Gene,phenotype_GD,Domain,AA_alt,Vartype,ID)
  
  
  data.df <- rbind(Patient_data.df %>% 
                     filter(!is.na(phenotype_GD)) %>% 
                     filter(Gene %in% data_pre_sel1.df$Gene,
                            Domain %in% data_pre_sel1.df$Domain,
                            phenotype_GD %in% data_pre_sel1.df$phenotype_GD,
                            AA_alt %in% data_pre_sel1.df$AA_alt,
                            Vartype %in% data_pre_sel1.df$Vartype,
                            ID%in% data_pre_sel1.df$ID) %>% 
                     mutate(group = "Group 1"),
                   Patient_data.df %>% 
                     filter(!is.na(phenotype_GD)) %>% 
                     filter(Gene %in% data_pre_sel2.df$Gene,
                            Domain %in% data_pre_sel2.df$Domain,
                            phenotype_GD %in% data_pre_sel2.df$phenotype_GD,
                            AA_alt %in% data_pre_sel2.df$AA_alt,
                            Vartype %in% data_pre_sel2.df$Vartype,
                            ID %in% data_pre_sel2.df$ID) %>% 
                     mutate(group = "Group 2")) %>% 
    rename(Phenotype2 = var) %>% 
    filter(!is.na(AA_pos),!is.na(Phenotype2)) %>% 
    filter(Phenotype2 != "Unknown") %>% 
    mutate(Phenotype2 = ifelse(Phenotype2 == "Yes"," Yes",Phenotype2)) %>% 
    group_by(Phenotype2,group) %>% summarise(n = n())  %>% ungroup()
  
  if(any(data.df$group == "Group 1") & any(data.df$group == "Group 2")){
    
    Group_colors <- c(color1,color2)
    
  }else if(any(data.df$group == "Group 1")){
    
    Group_colors <- c(color1)
    
  }else{
    Group_colors <- c(color2)
  }
  
  plot <- data.df %>%
    plot_ly(
      x=~Phenotype2, y=~n, type="bar", color=~group, alpha = 0.8, colors= Group_colors) %>%
    layout(title=title_sel ,font=plotly_font, yaxis = list(showline = T, title = "Number of patients"),
           margin = list(t =50),
           xaxis = list(title = "",  showline = T,tickangle = 45)
    ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  
  return(plot)
}

basic_onset_legend <- function(){
  legend <- data.frame(x=c(1,3), y=c(3,3), text=c("Missense", "PTV"))
  plot <- ggplot(legend, aes(x=x, y=y, color=text))+
    geom_point(size = 8)+
    scale_color_manual(values = c("Missense"="#BEBADA","PTV"="#FDB462"))+
    ylim(c(2.9,3.1))+
    xlim(c(0.8,5.2))+
    theme_void()+
    geom_text(aes(label=text), hjust=-0.2, color="black", size =5)+
    theme(legend.position = "none")
  
  
  return(plot)
}

####Whatever gnomAD
extract_gnomad_features <- function(Control_data.df,selected.df,variable,label){
  
  if(label == "exchange"){
    
    control_int.df <- Control_data.df %>% filter(Gene == selected.df$Gene[1], AA_pos == selected.df$AA_pos[1], AA_alt == selected.df$AA_alt[1])
    
    control_int2.df <- Control_data.df %>% filter(Gene == selected.df$Gene[1], AA_pos == selected.df$AA_pos[1])
    
    if(variable == "Allele count"){
      
      out <- ifelse(nrow(control_int.df) == 0,0,control_int.df$Allele_count)
      
    }else if(variable == "Allele freq"){
      
      out <- ifelse(nrow(control_int.df) == 0,0,control_int.df$Allele_freq)
      
    }
  }else if(label == "position"){
    
    control_int.df <- Control_data.df %>% filter(Gene == selected.df$Gene[1], AA_pos == selected.df$AA_pos[1])
    
    if(variable == "Allele count"){
      
      out <- ifelse(nrow(control_int.df) == 0,0,sum(control_int.df$Allele_count))
      
    }else if(variable == "Allele freq"){
      
      out <- ifelse(nrow(control_int.df) == 0,0,sum(control_int.df$Allele_freq))
      
    }
    
  }
  
  return(out)
  
  
}

####3D mappigns 
map_var_3d <- function(data,Gene_sel,gnomad_bool,pdb_sel,structure_coordinates,sub_color_i){
  
  variant.df <- data %>%
    filter(Gene == Gene_sel) %>%
    filter(Vartype == "Missense") %>%
    mutate(label = "pathogenic") %>%
    group_by(AA_pos,AA_ref,Gene) %>%
    summarise(n_occ = n()) %>%
    select(AA_pos,AA_ref,n_occ,Gene)
  
  gnomad.df <- Control_data.df %>%
    filter(Gene == Gene_sel) %>%
    group_by(AA_pos,AA_ref,Gene) %>%
    summarise(n_occ = n()) %>%
    select(AA_pos,AA_ref,n_occ,Gene)
  
  structure.df <- read_delim(structure_coordinates,delim = "\t") %>%
    mutate(Aminoacid = aaa(Aminoacid)) %>%
    select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
  
  variant.df <- variant.df %>%
    mutate(AA_ref = unlist(AA_ref)) %>%
    left_join(structure.df %>%
                unnest(Aminoacid),
              by = c("AA_pos" = "Uniprot_position", "AA_ref" = "Aminoacid", "Gene" = "gene")) %>%
    mutate(struc_cov = ifelse(is.na(Position_in_structure), "no", "yes")) %>%
    filter(struc_cov == "yes")
  
  gnomad.df <- gnomad.df %>%
    left_join(structure.df %>%
                unnest(Aminoacid),
              by = c("AA_pos" = "Uniprot_position", "AA_ref" = "Aminoacid", "Gene" = "gene")) %>%
    mutate(struc_cov = ifelse(is.na(Position_in_structure), "no", "yes")) %>%
    filter(struc_cov == "yes")
  
  sub_color <- c("#e2f970","#6fbbf7","#ee6c71","#ffbc5a","#bf73cc")
  sub_scale <- c(1.2,0.8)
  struc_color <- "wheat"
  
  rot = 270
  rot_axis = "x"
  spin_axis = "vy"
  
  #Specify yourself- color of the cartoon per subunit
  subunit_color <- c("wheat","white") #first color for GRIN1 second or GRIN2A
  
  print(variant.df)
  print("now gnomad")
  print(gnomad.df)
  
  #Model for the protein complex
  modelo <- r3dmol(
    viewer_spec = m_viewer_spec(
      cartoonQuality = 20,
      lowerZoomLimit = 50,
      upperZoomLimit = 1000
    )
  )
  
  modelo <- modelo %>% m_add_model(data = pdb_sel, format = "pdb")
  
  # Zoom to encompass the whole scene
  modelo <- modelo %>% m_zoom_to() %>%
    # Set color o cartoon representation
    m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
    # Set subunit colors
    m_set_style(
      sel = m_sel(chain = c("A")),
      style = m_style_cartoon(color = subunit_color[2])
    ) %>%
    # visualize variants grin1
    m_set_style(
      sel = m_sel(resi = variant.df$Position_in_structure,
                  atom = "CA",
                  chain = c("A")),
      style = m_style_sphere(colorScheme = NULL,
                             color = sub_color[sub_color_i],
                             scale = sub_scale[1])
    )
  
  if  (gnomad_bool == TRUE) {
    
    modelo <- modelo %>% m_set_style(
      sel = m_sel(resi = gnomad.df$Position_in_structure,
                  atom = "CA",
                  chain = c("A")),
      style = m_style_sphere(colorScheme = NULL,
                             color = "#333333",
                             scale = sub_scale[2]))
    
  }
  
  return(modelo) 
}      


#to 3 aa letter

aaa_new <- function(input){
  out <- sapply(input, function(x){
    if(!is.na(x)){
      aaa(x)
    }else{
      NA
    }
  })
  return(out)
}


aaa <- function(aa) {
  aa <- toupper(aa)
  aa_map <- sapply(aa, function(x) switch(x,
                                          A = "Ala",
                                          R = "Arg",
                                          N = "Asn",
                                          D = "Asp",
                                          C = "Cys",
                                          E = "Glu",
                                          Q = "Gln",
                                          G = "Gly",
                                          H = "His",
                                          I = "Ile",
                                          L = "Leu",
                                          K = "Lys",
                                          M = "Met",
                                          F = "Phe",
                                          P = "Pro",
                                          S = "Ser",
                                          T = "Thr",
                                          W = "Trp",
                                          Y = "Tyr",
                                          V = "Val",
                                          X = "Unknown"))
  return(aa_map)
}

##### Functions #####
#Basic Information
basic_gene1 = "GBA1"
basic_phenotype_fac = "Phenotype"
basic_phenotype_fac2 = "Phenotype_GD"
phenotype_name1 <- "Phenotype PD"
phenotype_name2 <- "Phenotype GD"
basic_phenotype_num = "#B2DF8A"
basic_phenotype_num1 = "blue"

basic_phenotype_colors <- rep(RColorBrewer::brewer.pal(12,"Set3"),2) #create more colors

colors_gene1 <-  c("#BEBADA","#FDB462","#BEBADA","#FDB462","#BEBADA","#FDB462","#BEBADA","#FDB462","#BEBADA","#FDB462")
##### Variant Analysis variable #####
variant_title1 <- "Phenotypes"

#####Research variable #####
color_list_research_compare <- c("green" = "#A6CEE3" ,"blue" = "#B2DF8A","red" = "#FB9A99","orange" = "#FDBF6F","purple" = "#CAB2D6")

Gene_colors <-  c("GBA1"="#6fbbf7", 
                  "Control" = "#333333",
                  "Other" = "#333333")

research_phenotype1_title <- "Number of patient variants per unit"
research_phenotype2_title <- "Phenotype PD"
research_phenotype3_title <- "Phenotype GD"
####################################################
# Onset_days_1.fun <- function(select_gene, colors_sel) {
# 
#   plot_ly(Patient_data.df %>%
#             filter(Gene == select_gene,
#                    !is.na(Earlyonset)) %>%
#             mutate(p_variant = Original_AA_change,
#                    p_variant = ifelse(is.na(p_variant), Vartype, p_variant)) %>%
#             filter(Vartype %in% c("Missense", "PTV")) %>%
#             group_by(Earlyonset) %>%
#             arrange(Earlyonset) %>%
#             ungroup(),
#           y = ~Earlyonset, type = "scatter", mode = "markers",  # Use scatter plot
#           marker = list(size = 10, opacity = 0.6, color = colors_sel),  # Adjust marker settings
#           hoverinfo = "text",
#           text = ~paste0(round(Earlyonset, digits = 2), " Days, ", p_variant)) %>%
#     layout(font = plotly_font,
#            title = "",
#            xaxis = list(title = "", tickangle = 45, showline = TRUE),
#            yaxis = list(type = "log",
#                         title = "Early-Onset",
#                         tickvals = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#                         tickmode = "array",
#                         showline = TRUE
#            )) %>%
#     config(modeBarButtonsToRemove = c("zoom2d", "pan2d", "select2d", "lasso2d", "autoScale2d", "resetScale2d", "toggleSpikelines", "zoomIn2d", "zoomOut2d", "hoverClosestCartesian", "hoverCompareCartesian"),
#            displaylogo = FALSE)
# 
# }

# Onset_days_1.fun <- function(select_gene, colors_sel) {
#   
#   plot_ly(Patient_data.df %>%
#             filter(Gene == select_gene,
#                    !is.na(Earlyonset)) %>%
#             mutate(p_variant = Original_AA_change,
#                    p_variant = ifelse(is.na(p_variant), Vartype, p_variant)) %>%
#             filter(Vartype %in% c("Missense", "PTV")) %>%
#             group_by(Earlyonset) %>%
#             arrange(Earlyonset) %>%
#             ungroup(),
#           x = ~Original_AA_change,
#           y = ~Earlyonset, type = "scatter", mode = "markers",  # Use connected scatter plot
#           marker = list(size = 10, opacity = 0.6, color = colors_sel),  # Adjust marker settings
#           hoverinfo = "text",
#           text = ~paste0(round(Earlyonset, digits = 2), " Days, ", p_variant)) %>%
#     layout(font = plotly_font,
#            title = "",
#            xaxis = list(title = "Early Onset of PD", tickangle = 45, showline = TRUE),
#            yaxis = list(type = "log",
#                         title = "Early-Onset",
#                         tickvals = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#                         tickmode = "array",
#                         showline = TRUE
#            )) %>%
#     config(modeBarButtonsToRemove = c("zoom2d", "pan2d", "select2d", "lasso2d", "autoScale2d", "resetScale2d", "toggleSpikelines", "zoomIn2d", "zoomOut2d", "hoverClosestCartesian", "hoverCompareCartesian"),
#            displaylogo = FALSE)
#   
# }
Onset_days_1.fun <- function(select_gene, colors_sel) {
  
  plot_ly(early_late_onset.df %>%
            filter(Gene == select_gene,
                   !is.na(Earlyonset)) %>%
            mutate(p_variant = Original_AA_change,
                   p_variant = ifelse(is.na(p_variant), Vartype, p_variant)) %>%
            filter(Vartype %in% c("Missense", "PTV")) %>%
            arrange(Earlyonset) %>%
            ungroup(),
          x = ~Original_AA_change,
          y = ~ifelse(Earlyonset == 0, 0, log(Earlyonset + 1)),  # Adjust y-axis representation
          type = "scatter", mode = "markers",  # Use connected scatter plot
          marker = list(size = 10, opacity = 0.6, color = colors_sel),  # Adjust marker settings
          hoverinfo = "text",
          text = ~paste0(round(Earlyonset, digits = 2), " Days, ", p_variant)) %>%
    layout(font = plotly_font,
           title = "",
           xaxis = list(title = "Early Onset of PD", tickangle = 45, showline = TRUE),
           yaxis = list(type = "linear",  # Use linear scale for y-axis
                        title = "Early-Onset",
                        showline = TRUE
           )) %>%
    config(modeBarButtonsToRemove = c("zoom2d", "pan2d", "select2d", "lasso2d", "autoScale2d", "resetScale2d", "toggleSpikelines", "zoomIn2d", "zoomOut2d", "hoverClosestCartesian", "hoverCompareCartesian"),
           displaylogo = FALSE)
  
}



Onset_days_2.fun <- function(select_gene, colors_sel) {
  
  plot_ly(early_late_onset.df %>%
            filter(Gene == select_gene,
                   !is.na(Lateonset)) %>%
            mutate(p_variant = Original_AA_change,
                   p_variant = ifelse(is.na(p_variant), Vartype, p_variant)) %>%
            filter(Vartype %in% c("Missense", "PTV")) %>%
            arrange(Lateonset) %>%
            ungroup(),
          x = ~Original_AA_change,
          y = ~ifelse(Lateonset == 0, 0, log(Lateonset + 1)),  # Adjust y-axis representation
          type = "scatter", mode = "markers",  # Use connected scatter plot
          marker = list(size = 10, opacity = 0.6, color = colors_sel),  # Adjust marker settings
          hoverinfo = "text",
          text = ~paste0(round(Lateonset, digits = 2), " Days, ", p_variant)) %>%
    layout(font = plotly_font,
           title = "",
           xaxis = list(title = "Late Onset of PD", tickangle = 45, showline = TRUE),
           yaxis = list(type = "linear",  # Use linear scale for y-axis
                        title = "Late-Onset",
                        showline = TRUE
           )) %>%
    config(modeBarButtonsToRemove = c("zoom2d", "pan2d", "select2d", "lasso2d", "autoScale2d", "resetScale2d", "toggleSpikelines", "zoomIn2d", "zoomOut2d", "hoverClosestCartesian", "hoverCompareCartesian"),
           displaylogo = FALSE)
  
}
