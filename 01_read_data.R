library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)

### list all files in the data directory

all_data_files <- list.files("data/", full.names = TRUE)


## read in all files and combine into table

all_data_terrestrial <- list()

for(f in all_data_files[str_detect(all_data_files, "Terrestrial")]){
  
  all_data_terrestrial[[f]] <- read.delim(f, sep = "|") %>% 
    mutate(aq_ter = "Terrestrial") %>% 
    as.data.table()
  
}

all_df_terrestrial <- rbindlist(all_data_terrestrial)

## clean column names

colnames_clean_terrestrial <- colnames(all_df_terrestrial) %>% 
  str_remove("X") %>% 
  str_replace_all("\\.\\.", "_") %>% 
  str_replace_all("\\.", "_") %>% 
  str_remove("^_") %>% 
  str_remove("_$")

colnames(all_df_terrestrial) <- colnames_clean_terrestrial

### filter papers for case study 

filtered_papers <- all_df_terrestrial[Chemical_Name == "(2E)-1-[(6-Chloro-3-pyridinyl)methyl]-N-nitro-2-imidazolidinimine" & Species_Scientific_Name == "Apis mellifera",]

write.csv(filtered_papers, "data_clean/apis_imidacloprid.csv", row.names = FALSE)

### select only relevant columns and add common names for neonicotinoids

terrestrial_all <- all_df_terrestrial[, .(CAS_Number, Chemical_Name, Species_Scientific_Name, Organism_Lifestage, Organism_Age_Mean, Organism_Age_Mean_Op,
                                  Exposure_Type, Media_Type, Test_Location, Number_of_Doses, 
                                  Observed_Response_Mean, Observed_Response_Min, Observed_Response_Max, Observed_Response_Units, Effect, Endpoint, Response_Site, aq_ter, Title, Publication_Year,
                                  Observed_Duration_Days, Observed_Duration_Units_Days),] %>% 
  mutate(Observed_Duration_Days = as.numeric(Observed_Duration_Days))

chemical_table <- data.frame(Chemical_Name = unique(terrestrial_all$Chemical_Name), Chemical_short = c("Dinotefuran", "Imidacloprid", "Imidaclothiz", "Nitenpyram", "Nithiazine", "Thiacloprid", 
                                                                                     "Thiamethoxam", "Clothianidin", "Acetamiprid"))
terrestrial_all <- terrestrial_all %>% 
  left_join(chemical_table)

## save compiled neonic data 

saveRDS(terrestrial_all, "data_clean/terrestrial_all.rds")

####
