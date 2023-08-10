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


################# CSS/STYLE/INFO #################
landing_panel <- "color: #333333;
      height: 200px;
      width: 260px"

spinner_color <- "#2196f3"

sub_style <- "color:gray; font-style: italic; font-size: 14px;"

schema_color <- "#B2DAB6"#"#7DCEA0"

schema_color_strong <- "#E9F7EF"#"#1E8449"

schema_color_light <- "gainsboro"#"#B2DAB6"

##### General Variable #####

contact_us <- "mailto:patrick.may@uni.lu"


##### Landing Page Variable #####
#General
landing_portal_title = "GBA1 Portal"

landing_bannername = "v2.gif"
landing_bannername1 = "DNA.png"


#Tabs
#top color
landing_top_color1 = "success"
landing_top_color2 = "warning"
landing_top_color3 = "royal"



#body color
landing_body_color1 = "#F8FCFE"
landing_body_color2 = "#fff8ef"
landing_body_color3 = "#f9f1fa"


landing_tab1 = "GBA1 genes, their function and associated disorders"
landing_tab2 = "Comprehensive information on variant interpretation"


#####Basic Information Variable #####
#General
basic_clinical_info_title <- p("Summary of clincal information curated in the", em("GBA1", style="font-style: italic;"),"Portal")


#Gene 1
gene1 = "GBA1"
gene1_info_header = div("Information about ", em("GBA1", style="font-style: italic;"), style="font-size: 17px;")

basic_text_title_1 <- p(em("GBA1", style="font-style: Italic;"),"related disorders")


#History not automized yet

Gene1_basic_text <- p(em("GBA1", style="font-style: Italic;")," ")

##### Variant Analysis variable #####

master.df <- read_delim("data/master_table4.txt", delim = "\t")
aa_pos_miss.df <- read_delim("data/aa_pos_miss.txt", delim = "\t")

var_possible_genes_title <- "Select gene"

var_possible_genes <- c("GBA1")

var_possible_phenotype <- c("Severe", "Mild", "Risk","Unknown", "Mild/severe")

var_patient_info_title <- h4("Individuals with the same variant in the", em("GBA1", style="font-style: italic;"), "Portal")
var_patient_info_abb <- "GBA: glucocerebrosidase gene; PD: Parkinson’s Disease ; GD: Gaucher disease;
DOI of Article 1: 10.1093/brain/awv179; DOI of Article 2: 10.1002/acn3.51164; DOI of Article 3: 10.1002/mgg3.267"


var_paralog_info_abb <- "Explain Abbreviations"
variants_compareplots_abb <- "DOI of Article 1: 10.1093/brain/awv179; DOI of Article 2: 10.1002/acn3.51164; DOI of Article 3: 10.1002/mgg3.267"

#####Research Variable #####

research_geno_transcripts <- "The following transcripts were used: GBA1: NM_001005742"

##### Registry Variable#####

Patient_data.df <- read_delim("data/master_table4.txt", delim = "\t")

#####################################
Patient_data_missense_only.df <- Patient_data.df %>%
  filter(Vartype == "Missense")


##### About Variable #####
about_terms_of_use <- p("All data here are publicly for the benefit of the wider biomedical community.
                               You can freely explore the data, and we encourage the use and publication of results generated
                               from these data. However, we encourage you to contact us before embarking on analyses to
                               check if your proposed analysis overlaps with work currently underway by our team. Further,
                               we request that you actively acknowledge and give attribution to the", em("GBA1", style="font-style: italic;") ,"Portal project, and
                               link back to the relevant page, wherever possible. All users of our data agree to not attempt
                               to reidentify participants. Our data set has been subjected to extensive quality control,
                               but may be imperfect so errors may remain. If you spot any results that seem impossible,
                               or have suggestions for", em("GBA1", style="font-style: italic;") ,"Portal improvements: ")
about_data_generation <- "Data has been curated in a community effort."


