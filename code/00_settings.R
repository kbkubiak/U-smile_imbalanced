source('code/00_functions.R')
#
# libraries ####
library(readxl)
library(writexl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(purrr)
library(pROC) # roc(), roc.test()
library(tidyr)
library(stringr)
library(lmtest) # for lrtest()
library(gridExtra)
library(gtable)
library(ggtext)
library(cowplot) # as_grob()
library(grid)
library(sjPlot)
library(rempsyc)
library(stargazer)
library(tinytex)
library(knitr)
library(broom)
library(kableExtra)
library(ggarchery)
#
# general settings ####
Sys.setenv(LANG = 'en')

# import raw data 
ds <- read_excel('data/raw_data_rnd.xlsx')
ds <- ds %>% select(disease, age, sex, bp, chol, 
                    glu, stde, rnd_normal, strat_rnd_normal)
# divide into nonevent and event datasets 
ds_nonev <- ds[ds$disease == 0, ]
ds_event <- ds[ds$disease == 1, ]
# ref, new variables #### 
ref_vars <- c("sex", "age", "bp", "chol") 
new_vars <- c("stde", "strat_rnd_normal", "glu", "rnd_normal")

prevalence <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
prevalence_breaks <- c(1, 10, 30, 50, 70, 90, 99)

# U-smile settings
usmile_vars_pretty <- c("ST~depression", "Str~Rnd~normal", "Glucose", "Rnd~normal")
subclass_order <- c('nonev_be',
                    'nonev_wo',
                    'event_wo',
                    'event_be')
usmile_colors <- c('nonev_be' = '#0F3C78',
                   'nonev_wo' = '#0F3C78',
                   'event_wo' = '#D51424',
                   'event_be' = '#D51424')
usmile_fills  <- c('nonev_be' = '#0F3C78',
                   'nonev_wo' = '#BED2FA',
                   'event_wo' = '#FBCDB9',
                   'event_be' = '#D51424')
usmile_labels <- c('non-events with better prediction',
                   'non-events with worse prediction',
                   'events with worse prediction',
                   'events with better prediction')
labels_prevalence <- c('1%','10%','30%','50%','70%','90%','99%')
labels_model_pretty <- c("ST depression", "Str Rnd normal", "Glucose", "Rnd normal")
#