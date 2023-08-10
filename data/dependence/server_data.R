############## DATA ############
#Load all possible exchanges 
all_exchanges.df <- read_delim("data/exchange.txt",delim = "\t") %>% 
  mutate(AA_ref = aaa(AA_ref),
         AA_alt = convert_aa(AA_alt))

#Load master table 
master.df <- read_delim("data/master_table4.txt", delim = "\t")

#load domain gene df
Domain_gene.df <- master.df %>% 
  distinct(Gene,AA_pos,Domain,Domain_color) %>% 
  ungroup()

#load Enzyme activity
paraz_mtr_input.df <- read_delim("data/enzyme_boxplot.txt", delim = "\t")
paraz_mtr_input_control.df <- read_delim("data/enzyme_control.txt", delim = "\t")

############################---------
#Load patient and control data

Patient_data.df <- read_delim("data/master_table4.txt", delim = "\t") %>% 
  select(-Transcript) %>% 
  mutate(AA_pos = as.numeric(AA_pos)) %>% 
  left_join(master.df %>% distinct(Transcript,Gene,AA_pos), by = c("AA_pos" = "AA_pos","Gene" = "Gene"), multiple = "all") %>% 
  left_join(all_exchanges.df %>% distinct(Gene,AA_pos = type.convert(AA_pos, as.is = TRUE)), by = c("AA_pos" = "AA_pos","Gene" = "Gene"), multiple = "all") %>%  
  mutate(AA_ref = ifelse(!is.na(AA_ref),AA_ref,"XXX") %>% aaa(), ##warnings due to none matching aminoacids are fine 
         AA_alt = ifelse(Vartype == "Missense",AA_alt,NA),
         cDNA = ifelse(!is.na(cDNA_pos), paste0("c.",cDNA_pos,cDNA_ref,">",cDNA_alt), "Not available"),
         Protein = AA_ori,
         Inheritance = ifelse(is.na(Inheritance),"Not available",Inheritance)) %>% 
  
  mutate(Phenotype = ifelse(str_detect(phenotype, "Early onset"), "EOEE", phenotype),
         phenotype = case_when(
           
           phenotype %in% c("Unknown", "Pathogenic/Pathogenic", "Pathogenic", "Risk") ~ Phenotype,
           TRUE ~ Phenotype
         ))


Patient_data_missense_only.df <- Patient_data.df %>% 
  filter(Vartype == "Missense") %>% 
  mutate(AA_alt_old = AA_alt,
         AA_alt = aaa_new(AA_alt))

Control_data.df <- read_delim("data/gnomad_variants.txt", delim = "\t") %>% 
  mutate(AA_ref = aaa(AA_ref),
         AA_alt = convert_aa(AA_alt)) 


### 3D_mapping Genotype interface ++
#PDB
pdb_sel_gene1 = "data/pdb/GBA.pdb1"
structure_coordinates_gene1 <- "data/pdb/6tjk.txt"

#Alphafold
pdb_sel_gene2 = "data/pdb/AF-A0A812PPN5-F1-model_v4.pdb1"
structure_coordinates_gene2 <- "data/pdb/alphafopld3.txt"

