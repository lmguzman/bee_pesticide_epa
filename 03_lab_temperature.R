library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)

## read the data

lab <- read.csv("data_clean/Mutiple_effects_toxicity - New_lab.csv")

## check how many did not record temperature not recorded

lab %>% 
  filter(Temperature_C =="NR" | Temperature_C == "" | Temperature_C == "*") %>% nrow()

## standardize temperatures characters

lab_temperature_clean <- lab %>%
  filter(Temperature_C != "" & Temperature_C !="NR") %>% 
  ### get mean and sd for values 25 ± 3
  mutate(t_p_m = as.numeric(str_extract(Temperature_C, '^\\d+\\.\\d+\\s|^\\d+\\s'))) %>% 
  mutate(t_p_sd1 = as.numeric(str_extract(str_extract(Temperature_C, "±\\s\\d+"), "\\d+"))) %>% 
  mutate(t_p_range1 = as.numeric(str_remove(str_extract(Temperature_C, '\\d+\\-'), "-")), 
         t_p_range2 = as.numeric(str_remove(str_extract(Temperature_C, '\\-\\d+'), "-"))) %>% 
  rowwise() %>% 
  mutate(t_p_m2 = mean(c(t_p_range1, t_p_range2)), tpsd2 = sd(c(t_p_range1, t_p_range2))) %>% 
  mutate(t_p_m3 = as.numeric(Temperature_C)) %>% 
  mutate(mean_temperature = case_when(!is.na(t_p_m) ~ t_p_m, 
                                      !is.na(t_p_m2) ~ t_p_m2,
                                      TRUE ~t_p_m3))  %>%
  mutate(sd_temperature = case_when(!is.na(t_p_sd1) ~ t_p_sd1, 
                                    !is.na(tpsd2) ~ tpsd2,
                                    TRUE ~NA)) %>% 
  select(-t_p_m, -t_p_sd1, -t_p_range1, -t_p_range2, -t_p_m2, -t_p_m3,-tpsd2 ) 


## get summary for temperatures 

lab_temperature_clean %>% 
  filter(!is.na(mean_temperature)) %>% 
  select(Endpoint, Observed_Response_Mean, Observed_Response_Units, Temperature_C, mean_temperature, sd_temperature) %>% 
  group_by(Endpoint) %>% 
  summarise(n(), min(mean_temperature), max (mean_temperature))

###### filter for studies looking at LD50 and standardize units ######

ld50_lab <- lab_temperature_clean %>% 
  filter(!is.na(mean_temperature)) %>% 
  select(Effect_Measurement, Endpoint, Observed_Response_Mean, Observed_Response_Units, Temperature_C, mean_temperature, sd_temperature) %>% 
  filter(Endpoint == "LD50" & Effect_Measurement == 'Mortality') %>% 
  mutate(Observed_Response_Mean = as.numeric(Observed_Response_Mean)) %>% 
  ## replace bee for org
  mutate(Observed_Response_Units = str_replace(Observed_Response_Units, "bee", "org")) %>% 
  ## change from ug to ng
  mutate(Observed_Response_Units_std = str_replace(Observed_Response_Units, "ug/org", "ng/org")) %>% 
  mutate(Observed_Response_Mean_std = ifelse(str_detect(Observed_Response_Units, "ug/org"), Observed_Response_Mean *1000, Observed_Response_Mean)) %>% 
  ## change from pg to ng
  mutate(Observed_Response_Units_std = str_replace(Observed_Response_Units_std, "pg/org", "ng/org")) %>% 
  mutate(Observed_Response_Mean_std = ifelse(str_detect(Observed_Response_Units, "pg/org"), Observed_Response_Mean *0.001, Observed_Response_Mean_std)) %>% 
  ## change from mg to ng
  mutate(Observed_Response_Units_std = str_replace(Observed_Response_Units_std, "mg/org", "ng/org")) %>% 
  mutate(Observed_Response_Mean_std = ifelse(str_detect(Observed_Response_Units, "mg/org"), Observed_Response_Mean *1000000, Observed_Response_Mean_std)) %>% 
  filter(str_detect(Observed_Response_Units_std, "ng/org"))
  
### look at studies that evaluated LC50 and standardize units 

lc50_lab <- lab_temperature_clean %>% 
  filter(!is.na(mean_temperature)) %>% 
  select(Effect_Measurement, Endpoint, Observed_Response_Mean, Observed_Response_Units, Temperature_C, mean_temperature, sd_temperature) %>% 
  filter(Endpoint == "LC50") %>% 
  filter(str_detect(Observed_Response_Units, "ng/org")) %>% 
  mutate(Observed_Response_Mean = as.numeric(Observed_Response_Mean)) %>% 
  mutate(Observed_Response_Mean_std = Observed_Response_Mean)

## join LD50 and LC50 results 

mortality <- bind_rows(ld50_lab, lc50_lab)

## compare quadratic and linear model 

model_quadratic <- lm(data = mortality, Observed_Response_Mean_std~ poly(mean_temperature,2))
anova(model_quadratic)

model_linear <- lm(data = mortality, Observed_Response_Mean_std~ mean_temperature)
anova(model_linear)

AIC(model_linear,model_quadratic)

anova(model_linear,model_quadratic)

## best model is quadratic

## visualize trend 

temperature_ld50 <- mortality %>%
  ggplot(aes(x = mean_temperature, y = Observed_Response_Mean_std)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y~ poly(x,2), se = FALSE) +
  theme_cowplot() +
  xlab("Mean Temperature") + ylab("LD50/LC50 (ng/org)")

ggsave(temperature_ld50, filename = 'figures/temperature_ld50.jpeg')
