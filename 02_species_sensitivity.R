library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(taxize)
library(purrr)
library(tidyr)
library(meta)
library(broom)

##### Sensitivity variation for all bees ####

##### get taxonomic classification for all species #####

spnames <- unique(terrestrial_all$Species_Scientific_Name)
out <- classification(spnames, db='ncbi')

df_filter <- lapply(out, FUN = function(x) is.data.frame(x))

out_only_df <- out[unlist(df_filter)]

species_family_names <- out_only_df %>% 
  map_df(~filter(.x, rank %in% c( 'family', 'genus')), .id = 'species') %>% 
  select(Species_Scientific_Name = species, name, rank)  %>% 
  pivot_wider(names_from = 'rank', values_from = 'name')

write.csv(species_family_names, "data_clean/species_names.csv", row.names = FALSE)

###### LD50 for all bees #######

## load compiled neonic data  ##

terrestrial_all <- readRDS("data_clean/terrestrial_all.rds")

## load species taxonomic data  ##

species_family_names <- read.csv("data_clean/species_names.csv")

## Filter families in Apidae, Megachilidae, and Halictidae, as well as LD50 endpoints 

bee_studies <- terrestrial_all %>% 
  left_join(species_family_names) %>% 
  filter(family %in% c('Apidae', 'Megachilidae', 'Halictidae')) %>% 
  filter(Endpoint == "LD50") %>% 
  mutate(Observed_Response_Units = ifelse(Title == 'Impacts of Chronic Sublethal Exposure to Clothianidin on Winter Honeybees',
                                          "AI ng/org", Observed_Response_Units)) %>% 
  mutate(Observed_Response_Mean = as.numeric(str_remove(Observed_Response_Mean, "\\/"))) %>% 
  mutate(Observed_Response_Min = as.numeric(str_remove(Observed_Response_Min, "\\/"))) %>% 
  mutate(Observed_Response_Max = as.numeric(str_remove(Observed_Response_Max, "\\/"))) 


###### standardize units ####

bee_studies_std <- bee_studies %>% 
  ## replace bee for org
  mutate(Observed_Response_Units = str_replace(Observed_Response_Units, "bee", "org")) %>% 
  ## change from ug to ng
  mutate(Observed_Response_Units_std = str_replace(Observed_Response_Units, "ug/org", "ng/org")) %>% 
  mutate(Observed_Response_Mean_std = ifelse(str_detect(Observed_Response_Units, "ug/org"), Observed_Response_Mean *1000, Observed_Response_Mean)) %>% 
  mutate(Observed_Response_Min_std = ifelse(str_detect(Observed_Response_Units, "ug/org"), Observed_Response_Min *1000, Observed_Response_Min)) %>% 
  mutate(Observed_Response_Max_std = ifelse(str_detect(Observed_Response_Units, "ug/org"), Observed_Response_Max *1000, Observed_Response_Max)) %>% 
  ## change from pg to ng
  mutate(Observed_Response_Units_std = str_replace(Observed_Response_Units_std, "pg/org", "ng/org")) %>% 
  mutate(Observed_Response_Mean_std = ifelse(str_detect(Observed_Response_Units, "pg/org"), Observed_Response_Mean *0.001, Observed_Response_Mean_std)) %>% 
  mutate(Observed_Response_Min_std = ifelse(str_detect(Observed_Response_Units, "pg/org"), Observed_Response_Min *0.001, Observed_Response_Min_std)) %>% 
  mutate(Observed_Response_Max_std = ifelse(str_detect(Observed_Response_Units, "pg/org"), Observed_Response_Max *0.001, Observed_Response_Max_std)) %>% 
  ## change from mg to ng
  mutate(Observed_Response_Units_std = str_replace(Observed_Response_Units_std, "mg/org", "ng/org")) %>% 
  mutate(Observed_Response_Mean_std = ifelse(str_detect(Observed_Response_Units, "mg/org"), Observed_Response_Mean *1000000, Observed_Response_Mean_std)) %>% 
  mutate(Observed_Response_Min_std = ifelse(str_detect(Observed_Response_Units, "mg/org"), Observed_Response_Min *1000000, Observed_Response_Min_std)) %>% 
  mutate(Observed_Response_Max_std = ifelse(str_detect(Observed_Response_Units, "mg/org"), Observed_Response_Max *1000000, Observed_Response_Max_std)) %>% 
  filter(Observed_Response_Units_std %in% c("ng/org", "AI ng/org")) %>%   
  filter(Chemical_short %in% c("Clothianidin", "Acetamiprid", "Imidacloprid", "Thiamethoxam"))

write.csv(bee_studies_std, "data_clean/all_bee_responses.csv")

################ READ DATA FROM HERE ########

######### adding sample sizes ########

### we added sample sizes and some standard errors manually

# data with all sample sizes

all_bee_responses <- read.csv("data_clean/all_bee_responses_RSS 2.csv") %>% 
  select(-ld50_se, -ld50_sd, -ld50_variance, -notes_ld50, -Observed_Response_Min, -Observed_Response_Max)

# data with better se

all_bee_responses_se <- read.csv("data_clean/bee_neonics_lethality_20230719_adding_extra_ci - bee_neonics_lethality_20230719_adding_extra_ci.csv") %>% 
  select(-sample_size)

joined_bee_responses <- all_bee_responses %>% 
  full_join(all_bee_responses_se) %>% 
  unique()
  
## re-standardizing the CI for the new added CI

bee_responses_with_ci <- joined_bee_responses %>% 
  ## change from ug to ng
  mutate(Observed_Response_Min_std = ifelse(str_detect(Observed_Response_Units, "ug/org"), Observed_Response_Min *1000, Observed_Response_Min)) %>% 
  mutate(Observed_Response_Max_std = ifelse(str_detect(Observed_Response_Units, "ug/org"), Observed_Response_Max *1000, Observed_Response_Max)) %>% 
  ## change from pg to ng
  mutate(Observed_Response_Min_std = ifelse(str_detect(Observed_Response_Units, "pg/org"), Observed_Response_Min *0.001, Observed_Response_Min_std)) %>% 
  mutate(Observed_Response_Max_std = ifelse(str_detect(Observed_Response_Units, "pg/org"), Observed_Response_Max *0.001, Observed_Response_Max_std)) %>% 
  ## change from mg to ng
  mutate(Observed_Response_Min_std = ifelse(str_detect(Observed_Response_Units, "mg/org"), Observed_Response_Min *1000000, Observed_Response_Min_std)) %>% 
  mutate(Observed_Response_Max_std = ifelse(str_detect(Observed_Response_Units, "mg/org"), Observed_Response_Max *1000000, Observed_Response_Max_std)) 

## add simplified exposure 

bee_studies_std_ex <- bee_responses_with_ci %>% 
  mutate(exposure_clean = case_when(Exposure_Type %in% 
                                      c("Topical, general", "Dermal", "Environmental, unspecified",
                                        "Spray, unspecified") ~ 
                                      "Topical",
                                    Exposure_Type %in% 
                                      c("Food", "Diet, unspecified", "Dropwise") ~ 
                                      "Consumption")) %>% 
  mutate(Observed_Duration_Days = round(Observed_Duration_Days)) %>% 
  mutate(result_id = 1:n())

### get proper se

# total 252 rows

## with confidence limits

# 162 with confidence intervals original

all_bee_responses_with_confidence <- bee_studies_std_ex %>%
  filter(!is.na(Observed_Response_Min_std))
  
## With standard errors
# 1 observation

all_bee_responses_with_se <- bee_studies_std_ex %>%
  filter(is.na(Observed_Response_Min_std)) %>% 
  filter(!is.na(`ld50_se_...CI`) & `ld50_se_...CI`!="") %>% ## change units if needed 
  ## change from ug to ng
  mutate(ld50_se_std = ifelse(str_detect(Observed_Response_Units, "ug/org"), `ld50_se_...CI` *1000, `ld50_se_...CI`)) %>% 
  ## change from pg to ng
  mutate(ld50_se_std = ifelse(str_detect(Observed_Response_Units, "pg/org"), `ld50_se_...CI` *0.001, ld50_se_std)) %>% 
  ## change from mg to ng
  mutate(ld50_se_std = ifelse(str_detect(Observed_Response_Units, "mg/org"), `ld50_se_...CI` *1000000, ld50_se_std)) 
  
## with standard deviation 
# 18 observations

all_bee_responses_with_sd <- bee_studies_std_ex %>%
  filter(is.na(Observed_Response_Min_std)) %>% 
  filter(!is.na(`ld50_sd_...CI`) & `ld50_sd_...CI`!="") %>% 
  ## change from ug to ng
  mutate(ld50_sd_std = ifelse(str_detect(Observed_Response_Units, "ug/org"), `ld50_sd_...CI` *1000, `ld50_sd_...CI`)) %>% 
  ## change from pg to ng
  mutate(ld50_sd_std = ifelse(str_detect(Observed_Response_Units, "pg/org"), `ld50_sd_...CI` *0.001, ld50_sd_std)) %>% 
  ## change from mg to ng
  mutate(ld50_sd_std = ifelse(str_detect(Observed_Response_Units, "mg/org"), `ld50_sd_...CI` *1000000, ld50_sd_std)) %>% 
  mutate(ld50_se_std_log = log(ld50_sd_std)/sample_size) %>%  ## calcualte standard error by dividing standard deviation by sample size 
  mutate(ld50_se_std = ld50_sd_std/sample_size) ## calcualte standard error by dividing standard deviation by sample size


#### with chi squared and pvalue 
## 29 responses 

all_bee_responses_with_pvalue <-bee_studies_std_ex %>%
  filter(is.na(Observed_Response_Min_std)) %>% 
  filter(!is.na(`Chisq`)) %>% 
  mutate(p_value = pchisq(Chisq, DF, lower.tail = FALSE)) %>% 
  mutate(p_value = ifelse(p_value ==1, p_value - 0.0001, p_value)) %>% 
  mutate(p_value = ifelse(p_value < 0.0001, 0.0001, p_value))


## rejoin data 

all_data_ci <- bind_rows(all_bee_responses_with_confidence, all_bee_responses_with_se, all_bee_responses_with_sd, all_bee_responses_with_pvalue)

#### run the meta analysis ####

########## re-run meta analysis separately for mode of exposure to look at duration and genus
# 8 models, 2 mode of exposure and 4 neonics
### run second model just with apis and then look at interaction between mode of exposure and neonic + fixed effect of duration

combinations_chemical_exposure <- all_data_ci %>% 
  select(Chemical_short, exposure_clean) %>% unique()

m.gen.reg.all.ci <- list()

for(i in 1:nrow(combinations_chemical_exposure)){
  
  data_to_use <- filter(all_data_ci, Chemical_short == combinations_chemical_exposure$Chemical_short[i]& 
                          exposure_clean == combinations_chemical_exposure$exposure_clean[i]) 
  
  m.gen <- metagen(TE = log(Observed_Response_Mean_std),
                   seTE = ld50_se_std_log,
                   studlab = Title,
                   lower = log(Observed_Response_Min_std),
                   upper = log(Observed_Response_Max_std),
                   pval = p_value, 
                   df = DF,
                   data = data_to_use,
                   sm = "MLN", ## log transformed means
                   title = "LD50", 
                   n.e = sample_size)
  
  if(length(unique(data_to_use$genus)) > 1){
    model_output <- metareg(m.gen, ~genus + Observed_Duration_Days)
  }else{
    model_output <- metareg(m.gen, ~ Observed_Duration_Days)
  }
  
  m.gen.reg.all.ci[[i]] <- tidy(model_output) %>% 
    mutate(significant = ifelse(`p.value` < 0.05, TRUE, FALSE), n = nrow(data_to_use), Chemical_short = combinations_chemical_exposure$Chemical_short[i],
           exposure_clean = combinations_chemical_exposure$exposure_clean[i])
  
  
}

m.gen.reg.all.ci


###### looking at apis only and the interaction between chemical, exposure, and duration

data_to_use <- filter(all_data_ci, genus == 'Apis')

m.gen <- metagen(TE = log(Observed_Response_Mean_std),
                 seTE = ld50_se_std,
                 studlab = Title,
                 lower = log(Observed_Response_Min_std),
                 upper = log(Observed_Response_Max_std),
                 data = data_to_use,
                 sm = "MLN", ## log transformed means
                 title = "LD50", 
                 n.e = sample_size)

model_output <- metareg(m.gen, ~ Observed_Duration_Days*Chemical_short + exposure_clean*Chemical_short)

### looking at apis life stage

data_to_use <- filter(all_data_ci, genus == 'Apis') %>% 
  mutate(larva = ifelse(Organism_Lifestage == "Larva", 'Larva', 'Not Larva')) %>% 
  filter(Chemical_short != 'Clothianidin') %>% 
  filter(Observed_Duration_Days < 5) %>% 
  filter(Chemical_short != 'Thiamethoxam') %>% 
  filter(exposure_clean == 'Consumption')

life_stage_models <- list()

for(c in unique(data_to_use$Chemical_short)){
  
  data_to_use_life <- data_to_use %>% 
    filter(Chemical_short == c)
  
  m.gen <- metagen(TE = log(Observed_Response_Mean_std),
                   seTE = ld50_se_std,
                   studlab = Title,
                   lower = log(Observed_Response_Min_std),
                   upper = log(Observed_Response_Max_std),
                   data = data_to_use_life,
                   sm = "MLN", ## log transformed means
                   title = "LD50", 
                   n.e = sample_size)
  
  life_stage_models[[c]] <- metareg(m.gen, ~ larva)
}




############## FIGURES ############

######## differences in sensitivity between bee genera ######

Chemicals <- bee_studies_std_ex$Chemical_short %>% unique()

## read in the data frame above but with sample sizes

list_plots_topical <- list()

list_plots_consumption <- list()

for(c in Chemicals){
  
  ## filter for pesticides 
  
  bee_chemical <- bee_studies_std_ex %>% 
    filter(Chemical_short == c)
  
  ## visualize with boxplot
  
  list_plots_topical[[c]]  <- bee_chemical %>% 
    filter(exposure_clean == "Topical") %>% 
    ggplot(aes(x = genus, y = Observed_Response_Mean_std, colour = as.factor(Observed_Duration_Days), size = sample_size)) + 
    geom_jitter() +
    scale_y_log10() +
    theme_cowplot() +
    xlab("") +
    ylab("") +
    ggtitle(c) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  list_plots_consumption[[c]]  <- bee_chemical %>% 
    filter(exposure_clean == "Consumption") %>% 
    ggplot(aes(x = genus, y = Observed_Response_Mean_std, colour = as.factor(Observed_Duration_Days), size = sample_size)) + 
    geom_jitter() +
    scale_y_log10() +
    theme_cowplot() +
    xlab("") +
    ylab("") +
    ggtitle(c) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}

all_plots_topical <- plot_grid(list_plots_topical[[1]], list_plots_topical[[2]], 
                               list_plots_topical[[3]], list_plots_topical[[4]])


all_plots_consumption <- plot_grid(list_plots_consumption[[1]], list_plots_consumption[[2]], 
                                   list_plots_consumption[[3]], list_plots_consumption[[4]])


ggsave(all_plots_topical, filename = "figures/topical.pdf", width = 15)

ggsave(all_plots_consumption, filename = "figures/consumption.pdf", width = 15)



######## differences in sensitivity between larva and not larvae in imidacloprid apis ######

data_to_use_larva <- filter(all_data_ci, genus == 'Apis') %>% 
  mutate(larva = ifelse(Organism_Lifestage == "Larva", 'Larva', 'Not Larva')) %>% 
  filter(Chemical_short %in% c('Imidacloprid', "Acetamiprid")) %>% 
  #filter(Observed_Duration_Days < 5) %>% 
 filter(exposure_clean == "Consumption") 

data_to_use_larva %>% 
  ggplot(aes(x = larva, y = Observed_Response_Mean_std, colour = Observed_Duration_Days, size = sample_size)) + 
  geom_jitter() +
  facet_wrap(~Chemical_short) +
  scale_y_log10() +
  theme_cowplot() +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## notes

## studies with replicate units
bee_studies %>% 
  filter(Title == 'Synergistic Mortality Between a Neonicotinoid Insecticide and an Ergosterol-Biosynthesis-Inhibiting Fungicide in Three Bee Species') %>% 
  select(Observed_Response_Mean, Observed_Response_Units, Chemical_short, genus, Observed_Duration_Days) %>% 
  arrange(genus, Chemical_short, Observed_Response_Units, Observed_Duration_Days) %>% nrow()

bee_studies %>% 
  filter(Title == 'Spray Toxicity and Risk Potential of 42 Commonly Used Formulations of Row Crop Pesticides to Adult Honey Bees (Hymenoptera: Apidae)') %>% 
  select(Observed_Response_Mean, Observed_Response_Units, Chemical_short) %>% 
  arrange(Chemical_short, Observed_Response_Units) %>% nrow()


######## compare plot with estimated from risk assessment point #######

epa_estimates <- read.csv("data_clean/epa_estimates.csv") %>% 
  ## change from ug to ng
  mutate(value = value *1000) %>% 
  mutate(new_units = 'ng/bee')

## read in the data frame above but with sample sizes


all_data_ci_f <- all_data_ci %>% 
  mutate(with_proper_ci = TRUE)

Chemicals <- all_data_ci_f$Chemical_short %>% unique()


list_plots_topical <- list()

list_plots_consumption <- list()

for(c in Chemicals){
  
  ## filter for pesticides 
  
  bee_chemical <- bee_studies_std_ex %>% 
    filter(Chemical_short == c)
  
  ## visualize with boxplot
  
  list_plots_topical[[c]]  <- bee_chemical %>% 
    filter(exposure_clean == "Topical") %>% 
    ggplot(aes(x = genus, y = Observed_Response_Mean_std, colour = as.factor(Observed_Duration_Days), size = sample_size)) + 
    geom_jitter() +
    scale_color_viridis_d() +
    scale_y_log10() +
    theme_cowplot() +
    ylab("") +
    ggtitle(c) +
    geom_point(data = filter(epa_estimates, Chemical_short == c, exposure_clean == 'contact', 
                             lifestage == 'adult'),
               aes(x="Apis", y = value), colour = 'red', size = 5, shape = 23, fill = 'red') +
    labs(color = "Duration (days)", size = "Sample Size")
  
  list_plots_consumption[[c]]  <- bee_chemical %>% 
    filter(exposure_clean == "Consumption") %>% 
    ggplot(aes(x = genus, y = Observed_Response_Mean_std, colour = as.factor(Observed_Duration_Days), size = sample_size)) + 
    geom_jitter() +
    scale_y_log10() +
    theme_cowplot() +
    scale_color_viridis_d() +
    ylab("") +
    ggtitle(c) +
    geom_point(data = filter(epa_estimates, Chemical_short == c, exposure_clean == 'oral', 
                             lifestage == 'adult'),
               aes(x="Apis", y = value), colour = 'red', size = 5, shape = 23, fill = 'red') +
    labs(color = "Duration (days)", size = "Sample Size")
  
}

all_plots_topical <- plot_grid(list_plots_topical[[1]], list_plots_topical[[2]], 
                               list_plots_topical[[3]], list_plots_topical[[4]])

ggsave(all_plots_topical, filename  = "figures/epa_comparison_all_topical.pdf" , width = 15)

all_plots_consumption <- plot_grid(list_plots_consumption[[1]], list_plots_consumption[[2]], 
                                   list_plots_consumption[[3]], list_plots_consumption[[4]])


ggsave(all_plots_consumption, filename  = "figures/epa_comparison_all_consumption.pdf" , width = 15)



######## compare plot with estimated from risk assessment line plot #######

epa_estimates <- read.csv("data_clean/epa_estimates.csv") %>% 
  ## change from ug to ng
  mutate(value = value *1000) %>% 
  mutate(new_units = 'ng/bee')
  
## read in the data frame above but with sample sizes


all_data_ci_f <- all_data_ci %>% 
  mutate(with_proper_ci = TRUE)

Chemicals <- all_data_ci_f$Chemical_short %>% unique()


list_plots_topical <- list()

list_plots_consumption <- list()

for(c in Chemicals){
  
  ## filter for pesticides 
  
  bee_chemical <- bee_studies_std_ex %>% 
    filter(Observed_Duration_Days < 3) %>% 
    filter(Organism_Lifestage %in% c('Adult', "Not reported")) %>% 
    filter(Species_Scientific_Name %in% c("Apis mellifera", "Apis mellifera ligustica")) %>% 
    filter(Exposure_Type %in% c('Dermal', 'Diet, unspecified', 'Topical, general', 'Food', 'Spray, unspecified')) %>% 
    filter(Chemical_short == c)
  
  ## visualize with boxplot
  
  list_plots_topical[[c]]  <- bee_chemical %>% 
    filter(exposure_clean == "Topical") %>% 
    ggplot(aes(x = sample_size, y = Observed_Response_Mean_std)) + 
    geom_jitter() +
    scale_y_log10() +
    theme_cowplot() +
    ylab("") +
    xlab("Sample size") +
    ggtitle(c) +
    geom_hline(data = filter(epa_estimates, Chemical_short == c, exposure_clean == 'contact', 
                             lifestage == 'adult'),
               aes(yintercept = value), colour = 'red')
  
  list_plots_consumption[[c]]  <- bee_chemical %>% 
    filter(exposure_clean == "Consumption") %>% 
    ggplot(aes(x = sample_size, y = Observed_Response_Mean_std)) + 
    geom_jitter() +
    scale_y_log10() +
    theme_cowplot() +
    ylab("") +
    xlab("Sample size") +
    ggtitle(c) +
    geom_hline(data = filter(epa_estimates, Chemical_short == c, exposure_clean == 'oral', 
                             lifestage == 'adult'),
               aes(yintercept = value), colour = 'red')
  
}

all_plots_topical <- plot_grid(list_plots_topical[[1]], list_plots_topical[[2]], 
                               list_plots_topical[[3]], list_plots_topical[[4]])

ggsave(all_plots_topical, filename  = "figures/epa_comparison_subset_all_topical.pdf" , width = 15)


all_plots_consumption <- plot_grid(list_plots_consumption[[1]], list_plots_consumption[[2]], 
                                   list_plots_consumption[[3]])


ggsave(all_plots_consumption, filename  = "figures/epa_comparison_subset_all_consumption.pdf" , width = 15)


epa_estimates <- epa_estimates %>% 
  mutate(exposure_clean = ifelse(exposure_clean == 'oral', 'Consumption', "Topical")) %>% 
  filter(lifestage == 'adult')

bee_studies_std_ex %>% 
  filter(Observed_Duration_Days < 3) %>% 
  filter(Organism_Lifestage %in% c('Adult', "Not reported")) %>% 
  filter(Species_Scientific_Name %in% c("Apis mellifera", "Apis mellifera ligustica")) %>% 
  filter(Exposure_Type %in% c('Dermal', 'Diet, unspecified', 'Topical, general', 'Food', 'Spray, unspecified')) %>% 
  group_by(Chemical_short, exposure_clean) %>% 
  summarise(mean = mean(Observed_Response_Mean_std, na.rm = TRUE), median = median(Observed_Response_Mean_std, na.rm = TRUE), 
            weighted_mean = sum(Observed_Response_Mean_std*sample_size, na.rm = TRUE)/sum(sample_size, na.rm = TRUE)) %>% 
  left_join(epa_estimates) %>% 
  select(Chemical_short, exposure_clean, mean, median, weighted_mean, epa_estimate = value)


#### Summary results for differences in order of magnitude Apis #### 

bee_studies_std_ex %>% 
  filter(genus == 'Apis') %>% 
  group_by(Chemical_short, exposure_clean) %>% 
  summarise(min = min(Observed_Response_Mean_std, na.rm = TRUE), max = max(Observed_Response_Mean_std, na.rm = TRUE)) %>% 
  mutate(exponent_min = as.numeric(str_extract(formatC(min, format = "e", digits = 2), "\\+\\d+|\\-\\d+"))) %>% 
  mutate(exponent_max = as.numeric(str_extract(formatC(max, format = "e", digits = 2), "\\+\\d+|\\-\\d+"))) %>% 
  mutate(differences_order_magnitude = exponent_min - exponent_max) %>% View()
  

#### Summary results for underestimation for other bees ####

bee_studies_std_ex %>% 
  group_by(Chemical_short, exposure_clean, genus) %>% 
  summarise(median = median(Observed_Response_Mean_std, na.rm = TRUE)) %>% 
  filter(genus != "Apis") %>% 
  left_join(epa_estimates %>% 
              filter(length == 'acute', lifestage == "adult") %>% 
              mutate(exposure_clean = ifelse(exposure_clean == "oral", "Consumption", "Topical")) %>% 
              select(Chemical_short, exposure_clean, epa_value = value)) %>% 
  mutate(exponent_median = as.numeric(str_extract(formatC(median, format = "e", digits = 2), "\\+\\d+|\\-\\d+"))) %>% 
  mutate(exponent_epa = as.numeric(str_extract(formatC(epa_value, format = "e", digits = 2), "\\+\\d+|\\-\\d+"))) %>% 
  mutate(differences_order_magnitude = exponent_median - exponent_epa) %>% 
  filter(differences_order_magnitude < 0) 

formatC(numb, format = "e", digits = 2)
