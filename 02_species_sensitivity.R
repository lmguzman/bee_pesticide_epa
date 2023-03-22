library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)

## load compiled neonic data  ##

terrestrial_all <- readRDS("data_clean/terrestrial_all.rds")

##### Sensitivity variation within apis #####

## filter studies that are done with Apis and measure LD50

apis_studies <- terrestrial_all %>% 
  filter(str_detect(Species_Scientific_Name, "Apis") & Endpoint == "LD50") %>% 
  mutate(Observed_Response_Mean = as.numeric(str_remove(Observed_Response_Mean, "\\/"))) %>% 
  mutate(Observed_Response_Min = as.numeric(str_remove(Observed_Response_Min, "\\/"))) %>% 
  mutate(Observed_Response_Max = as.numeric(str_remove(Observed_Response_Max, "\\/"))) 

## standardize units for ng/org

apis_studies_std_ld50<- apis_studies %>% 
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
  filter(Observed_Response_Units_std %in% c("ng/org", "AI ng/org")) 

#### plot the species sentivitiy curve within apis ###

apis_studies_plot <- apis_studies_std_ld50 %>% 
  filter(Chemical_short %in% c("Acetamiprid", "Clothianidin", 
                               "Imidacloprid", "Thiamethoxam")) %>% 
  mutate(Species_Scientific_Name = str_replace(Species_Scientific_Name, "Apis", "A.")) %>% 
  group_by(Chemical_short) %>% 
  arrange(Observed_Response_Mean_std) %>% 
  mutate(rank = 1:n()) %>%  
  ggplot() + 
  geom_point(aes(x = Observed_Response_Mean_std, y = rank, colour = Species_Scientific_Name), size = 3) +
  geom_errorbarh(aes(xmin = Observed_Response_Min_std, xmax = Observed_Response_Max_std, y = rank, 
                     colour = Species_Scientific_Name), height = 0) +
  scale_x_log10() +
  facet_wrap(~Chemical_short, ncol = 1) +
  theme_cowplot() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15), 
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 15), 
        legend.position = 'bottom') +
  scale_color_discrete(name = "") +
  ylab("Study") +
  xlab("Compound concentration log ng/org") +
  guides(color=guide_legend(nrow=4,byrow=TRUE))

ggsave(apis_studies_plot, filename = "figures/apis_ssd.jpeg", height = 12)

### summarise apis results ##

apis_studies_std_ld50 %>% 
  filter(Chemical_short %in% c("Acetamiprid", "Clothianidin", 
                               "Imidacloprid", "Thiamethoxam")) %>% 
  group_by(Chemical_short) %>% 
  summarise(n(), mean = round(mean(Observed_Response_Mean_std),2), var_means = round(sd(Observed_Response_Mean_std),2)) %>% View()


##### Sensitivity variation for all bees #####

library(taxize)
library(purrr)
library(tidyr)

## get taxonomic classification for all species ##

spnames <- unique(terrestrial_all$Species_Scientific_Name)
out <- classification(spnames, db='ncbi')

df_filter <- lapply(out, FUN = function(x) is.data.frame(x))

out_only_df <- out[unlist(df_filter)]

species_family_names <- out_only_df %>% 
  map_df(~filter(.x, rank %in% c( 'family', 'genus')), .id = 'species') %>% 
  select(Species_Scientific_Name = species, name, rank)  %>% 
  pivot_wider(names_from = 'rank', values_from = 'name')

## filter for bees and for end point in LD50 

bee_studies <- terrestrial_all %>% 
  left_join(species_family_names) %>% 
  filter(family %in% c('Apidae', 'Megachilidae', 'Halictidae')) %>% 
  filter(Endpoint == "LD50") %>% 
  mutate(Observed_Response_Mean = as.numeric(str_remove(Observed_Response_Mean, "\\/"))) %>% 
  mutate(Observed_Response_Min = as.numeric(str_remove(Observed_Response_Min, "\\/"))) %>% 
  mutate(Observed_Response_Max = as.numeric(str_remove(Observed_Response_Max, "\\/"))) %>%   
  filter(Chemical_short %in% c("Clothianidin", "Acetamiprid", "Imidacloprid", "Thiamethoxam"))

### standardize units

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
  filter(Observed_Response_Units_std %in% c("ng/org", "AI ng/org")) 


#### plot the bee sensitivity aggregating by species ###

bee_studies_plot <- bee_studies_std%>% 
  group_by(Chemical_short, Species_Scientific_Name) %>%
  summarise(mean_response = mean(Observed_Response_Mean_std), sd_response = sd(Observed_Response_Mean_std)) %>% 
  mutate(sd_response = ifelse(is.na(sd_response), 0, sd_response)) %>% 
  ungroup() %>% 
  group_by(Chemical_short) %>% 
  arrange(mean_response) %>% 
  left_join(species_family_names) %>% 
  mutate(rank = 1:n()) %>%  
  ggplot() + 
  geom_point(aes(x = mean_response, y = rank, colour = genus, shape = family), size = 3) +
  geom_errorbarh(aes(xmin = mean_response-sd_response, xmax = mean_response+sd_response, y = rank, 
                     colour = genus), height = 0) +
  scale_x_log10() +
  facet_wrap(~Chemical_short, ncol = 1) +
  theme_cowplot() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15), 
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 15), 
        legend.position = 'bottom')+
  scale_color_discrete(name = "") + scale_shape(name = "") +
  ylab("Species") +
  xlab("Compound concentration log ng/org")+
  guides(color=guide_legend(nrow=4,byrow=TRUE), 
         shape = guide_legend(nrow = 2, byrow = TRUE))

ggsave(bee_studies_plot, filename = "figures/bee_ssd.jpeg", height = 15)


## combine both plots for Figure 3 

sentitivity_all <- plot_grid(apis_studies_plot, bee_studies_plot, labels = c("A.", "B."))

ggsave(sentitivity_all, filename = "figures/sentitivity_all.jpeg", height = 15, width = 12)


##### Differences in sensitivity between Apis #######


## run anova for each pesticide

Chemicals <- c("Acetamiprid", "Clothianidin", 
               "Imidacloprid", "Thiamethoxam")

for(c in Chemicals) {
  
  # filter for each pesticide
  bee_chemical <- apis_studies_std_ld50 %>% 
    filter(Chemical_short == c)
  
  # run anova
  species_model <- aov(log(Observed_Response_Mean_std) ~ Species_Scientific_Name, data = bee_chemical)
  print(c)
  print(summary(species_model))
  
  # run tukey test
  tukey_diff <- TukeyHSD(species_model, conf.level=.95) 
  
  tukey_diff$Species_Scientific_Name %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column('comparison') %>% 
    filter(`p adj` < 0.05) %>% 
    print()
  
}

## visualize observations among apis

list_boxplots_apis <- list()

for(c in Chemicals){
  
  ## filter for pesticide
  
  bee_chemical <- apis_studies_std_ld50 %>% 
    filter(Chemical_short == c)
  
  ## visualize with boxplot
  
  list_boxplots_apis[[c]]  <- bee_chemical %>% 
    mutate(species = str_replace(Species_Scientific_Name, "Apis", "A.")) %>% 
    ggplot(aes(x = species, y = Observed_Response_Mean_std)) + 
    geom_boxplot() +
    scale_y_log10() +
    theme_cowplot() +
    xlab("") +
    ylab("") +
    ggtitle(c) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}

all_boxplots_apis <- plot_grid(list_boxplots_apis[[1]], list_boxplots_apis[[2]], 
                          list_boxplots_apis[[3]], list_boxplots_apis[[4]])

ggsave(all_boxplots_apis, filename = "figures/apis_boxplots.jpeg", height = 10)


######## differences in sensitivity between bee genera ######

## run anova for each pesticide

Chemicals <- bee_studies_std$Chemical_short %>% unique()

for(c in Chemicals) {
  
  ## filter for each pesticide
  
  bee_chemical <- bee_studies_std %>% 
    filter(Chemical_short == c)
  
  ## run anova
  
  genus_model <- aov(log(Observed_Response_Mean_std) ~ genus, data = bee_chemical)
  print(c)
  print(summary(genus_model))
  
  ## run tukey test
  
  tukey_diff <- TukeyHSD(genus_model, conf.level=.95) 
  
  tukey_diff$genus %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column('comparison') %>% 
    filter(`p adj` < 0.05) %>% 
    print()
  
}


## visualize differences between bee genera

list_boxplots <- list()

for(c in Chemicals){
  
  ## filter for pesticides 
  
  bee_chemical <- bee_studies_std %>% 
    filter(Chemical_short == c)
  
  ## visualize with boxplot
  
  list_boxplots[[c]]  <- bee_chemical %>% 
    ggplot(aes(x = genus, y = Observed_Response_Mean_std)) + 
    geom_boxplot() +
    scale_y_log10() +
    theme_cowplot() +
    xlab("") +
    ylab("") +
    ggtitle(c) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}

all_boxplots <- plot_grid(list_boxplots[[1]], list_boxplots[[2]], 
                          list_boxplots[[3]], list_boxplots[[4]])

ggsave(all_boxplots, filename = "figures/genus_boxplots.jpeg")



