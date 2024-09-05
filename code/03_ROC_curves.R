source('code/00_settings.R')
#
raw_indices <- read.csv("results/raw_indices.txt", sep="")
#
# find iteration corresponding to mean BS_ref for every prevalence ####
results_ref_mean %>% filter(dataset == 'train') %>% select(dataset, prevalence, BS_ref) %>% select(BS_ref) -> mean_BS_ref
mean_BS_ref <- as.data.frame(mean_BS_ref)
mean_BS_ref$BS_ref

raw_results_ref_conv %>% filter(dataset == 'train', prevalence == 0.01) %>% 
  reframe(iteration = iteration, BS_ref = BS_ref, difference = abs(BS_ref - mean_BS_ref$BS_ref[1])) %>% 
  arrange(abs(difference)) # iter 628 

raw_results_ref_conv %>% filter(dataset == 'train', prevalence == 0.10) %>% 
  reframe(iteration = iteration, BS_ref = BS_ref, difference = abs(BS_ref - mean_BS_ref$BS_ref[2])) %>% 
  arrange(abs(difference)) # iter 945

raw_results_ref_conv %>% filter(dataset == 'train', prevalence == 0.30) %>% 
  reframe(iteration = iteration, BS_ref = BS_ref, difference = abs(BS_ref - mean_BS_ref$BS_ref[3])) %>% 
  arrange(abs(difference)) # iter  381 

raw_results_ref_conv %>% filter(dataset == 'train', prevalence == 0.50) %>% 
  reframe(iteration = iteration, BS_ref = BS_ref, difference = abs(BS_ref - mean_BS_ref$BS_ref[4])) %>% 
  arrange(abs(difference)) # iter  544 

raw_results_ref_conv %>% filter(dataset == 'train', prevalence == 0.70) %>% 
  reframe(iteration = iteration, BS_ref = BS_ref, difference = abs(BS_ref - mean_BS_ref$BS_ref[5])) %>% 
  arrange(abs(difference)) # iter  152 

raw_results_ref_conv %>% filter(dataset == 'train', prevalence == 0.90) %>% 
  reframe(iteration = iteration, BS_ref = BS_ref, difference = abs(BS_ref - mean_BS_ref$BS_ref[6])) %>% 
  arrange(abs(difference)) # iter  901 

raw_results_ref_conv %>% filter(dataset == 'train', prevalence == 0.99) %>% 
  reframe(iteration = iteration, BS_ref = BS_ref, difference = abs(BS_ref - mean_BS_ref$BS_ref[7])) %>% 
  arrange(abs(difference)) # iter  47 

# indices for iteration corresponding to the mean BS_ref for every prevalence 
iter_train <- c(628, 945, 381, 544, 152, 901, 47)

ind_prev01 <- raw_indices %>% filter(dataset == 'train', prevalence == 0.01, iteration == iter_train[1])
ind_prev10 <- raw_indices %>% filter(dataset == 'train', prevalence == 0.10, iteration == iter_train[2])
ind_prev30 <- raw_indices %>% filter(dataset == 'train', prevalence == 0.30, iteration == iter_train[3])
ind_prev50 <- raw_indices %>% filter(dataset == 'train', prevalence == 0.50, iteration == iter_train[4])
ind_prev70 <- raw_indices %>% filter(dataset == 'train', prevalence == 0.70, iteration == iter_train[5])
ind_prev90 <- raw_indices %>% filter(dataset == 'train', prevalence == 0.90, iteration == iter_train[6])
ind_prev99 <- raw_indices %>% filter(dataset == 'train', prevalence == 0.99, iteration == iter_train[7])

ind_nonev_prev01 <- ind_prev01[ind_prev01$class == 0,'index']
ind_nonev_prev10 <- ind_prev10[ind_prev10$class == 0,'index']
ind_nonev_prev30 <- ind_prev30[ind_prev30$class == 0,'index']
ind_nonev_prev50 <- ind_prev50[ind_prev50$class == 0,'index']
ind_nonev_prev70 <- ind_prev70[ind_prev70$class == 0,'index']
ind_nonev_prev90 <- ind_prev90[ind_prev90$class == 0,'index']
ind_nonev_prev99 <- ind_prev99[ind_prev99$class == 0,'index']

ind_event_prev01 <- ind_prev01[ind_prev01$class == 1,'index']
ind_event_prev10 <- ind_prev10[ind_prev10$class == 1,'index']
ind_event_prev30 <- ind_prev30[ind_prev30$class == 1,'index']
ind_event_prev50 <- ind_prev50[ind_prev50$class == 1,'index']
ind_event_prev70 <- ind_prev70[ind_prev70$class == 1,'index']
ind_event_prev90 <- ind_prev90[ind_prev90$class == 1,'index']
ind_event_prev99 <- ind_prev99[ind_prev99$class == 1,'index']

ds_prev01 <- rbind(ds_nonev[ind_nonev_prev01,], ds_event[ind_event_prev01,])
ds_prev10 <- rbind(ds_nonev[ind_nonev_prev10,], ds_event[ind_event_prev10,])
ds_prev30 <- rbind(ds_nonev[ind_nonev_prev30,], ds_event[ind_event_prev30,])
ds_prev50 <- rbind(ds_nonev[ind_nonev_prev50,], ds_event[ind_event_prev50,])
ds_prev70 <- rbind(ds_nonev[ind_nonev_prev70,], ds_event[ind_event_prev70,])
ds_prev90 <- rbind(ds_nonev[ind_nonev_prev90,], ds_event[ind_event_prev90,])
ds_prev99 <- rbind(ds_nonev[ind_nonev_prev99,], ds_event[ind_event_prev99,])

ds_prev01$prevalence <- rep(0.01, 300)
ds_prev10$prevalence <- rep(0.10, 300)
ds_prev30$prevalence <- rep(0.30, 300)
ds_prev50$prevalence <- rep(0.50, 300)
ds_prev70$prevalence <- rep(0.70, 300)
ds_prev90$prevalence <- rep(0.90, 300)
ds_prev99$prevalence <- rep(0.99, 300)

dataset_prev <- rbind(ds_prev01,
                      ds_prev10,
                      ds_prev30,
                      ds_prev50,
                      ds_prev70,
                      ds_prev90,
                      ds_prev99)
#
# build reference models for each prevalence: ref_models_prev ####
prevalence <- c(0.01, 0.10, 0.30, 0.50, 0.70, 0.90, 0.99)
labels_prevalence

ref_models_prev <- list()
for(i in seq_along(prevalence)){
  ref_models_prev[[i]] <- glm(disease ~ sex + age + bp + chol, family = binomial, data = dataset_prev[dataset_prev$prevalence == prevalence[i],])
}
names(ref_models_prev) <- labels_prevalence
#
# ROC curves for ref models: roc_ref_prev ####
roc_ref_prev <- lapply(ref_models_prev, function(x) roc(x$y, x$fitted.values, direction = '<'))
#
# build new models for each prevalence: new_models_XXX ####
# glucose
new_models_glu <- list()
for(i in seq_along(prevalence)){
  new_models_glu[[i]] <- glm(disease ~ sex + age + bp + chol + glu, family = binomial, data = dataset_prev[dataset_prev$prevalence == prevalence[i],])
}
names(new_models_glu) <- labels_prevalence

# stde
new_models_stde <- list()
for(i in seq_along(prevalence)){
  new_models_stde[[i]] <- glm(disease ~ sex + age + bp + chol + stde, family = binomial, data = dataset_prev[dataset_prev$prevalence == prevalence[i],])
}
names(new_models_stde) <- labels_prevalence

# rnd_normal
new_models_rnd_normal <- list()
for(i in seq_along(prevalence)){
  new_models_rnd_normal[[i]] <- glm(disease ~ sex + age + bp + chol + rnd_normal, family = binomial, data = dataset_prev[dataset_prev$prevalence == prevalence[i],])
}
names(new_models_rnd_normal) <- labels_prevalence

# strat_rnd_normal
new_models_strat_rnd_normal <- list()
for(i in seq_along(prevalence)){
  new_models_strat_rnd_normal[[i]] <- glm(disease ~ sex + age + bp + chol + strat_rnd_normal, family = binomial, data = dataset_prev[dataset_prev$prevalence == prevalence[i],])
}
names(new_models_strat_rnd_normal) <- labels_prevalence
#
# ROC curves for new models: roc_new_XXX ####
roc_new_glu <- lapply(new_models_glu, function(x) roc(x$y, x$fitted.values, direction = '<'))
roc_new_stde <- lapply(new_models_stde, function(x) roc(x$y, x$fitted.values, direction = '<'))
roc_new_rnd_normal <- lapply(new_models_rnd_normal, function(x) roc(x$y, x$fitted.values, direction = '<'))
roc_new_strat_rnd_normal <- lapply(new_models_strat_rnd_normal, function(x) roc(x$y, x$fitted.values, direction = '<'))
#
# LRT and delta_AUC: p values from the mean value of the test statistic ####
results_new_mean %>% filter(dataset == 'train') %>% 
  arrange(model) %>% 
  mutate(LRT_p = round(2 * (1 - pchisq(LRT_statistic, 1)), 3),
         AUC_p = round(2 * pnorm(AUC_statistic), 3 ) ) %>% 
  select(dataset, prevalence, model, LRT_p_value, AUC_p_value) -> p_values
#
# panel titles ####
panel_titles_glu <- paste0(labels_prevalence, 
                           ' Glucose', 
                           ifelse(filter(p_values, model == 'glu')$LRT_p_value < 0.05, " *", ""),
                           ifelse(filter(p_values, model == 'glu')$AUC_p_value < 0.05, "<sup>, #</sup>", "")    )
panel_titles_stde <- paste0(labels_prevalence, 
                            ' ST depression', 
                            ifelse(filter(p_values, model == 'stde')$LRT_p_value < 0.05, " *", ""),
                            ifelse(filter(p_values, model == 'stde')$AUC_p_value < 0.05, "<sup>, #</sup>", "")    )
panel_titles_rnd_normal <- paste0(labels_prevalence, 
                                  ' Rnd normal', 
                                  ifelse(filter(p_values, model == 'rnd_normal')$LRT_p_value < 0.05, " *", ""),
                                  ifelse(filter(p_values, model == 'rnd_normal')$AUC_p_value < 0.05, "<sup>, #</sup>", "")    )
panel_titles_strat_rnd_normal <- paste0(labels_prevalence, 
                                        ' Str Rnd normal', 
                                        ifelse(filter(p_values, model == 'strat_rnd_normal')$LRT_p_value < 0.05, " *", ""),
                                        ifelse(filter(p_values, model == 'strat_rnd_normal')$AUC_p_value < 0.05, "<sup>, #</sup>", "")    )
#
# plots of ROC curves ####
# glucose
ROC_plots_glu <- list()
for(i in seq_along(roc_ref_prev)){
  ROC_plots_glu[[i]] <- ggroc(list(roc_ref_prev[[i]], roc_new_glu[[i]]),
                              legacy.axes = TRUE) +
    scale_colour_manual(values = c('black', '#D51424')) +
    geom_abline(intercept = 0, slope = 1, color = 'grey45') +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(aspect.ratio = 1,
          legend.position = 'none',
          plot.title.position = "plot", 
          plot.title = element_markdown(size =   6),
          axis.title.x = element_markdown(size = 6), 
          axis.title.y = element_markdown(size = 6), 
          axis.text = element_text(size = 6),
          plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), 'cm')
    ) +
    ggtitle(panel_titles_glu[i])
  
}

# stde
ROC_plots_stde <- list()
for(i in seq_along(roc_ref_prev)){
  ROC_plots_stde[[i]] <- ggroc(list(roc_ref_prev[[i]], roc_new_stde[[i]]),
                               legacy.axes = TRUE) +
    scale_colour_manual(values = c('black', '#D51424')) +
    geom_abline(intercept = 0, slope = 1, color = 'grey45') +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(aspect.ratio = 1,
          legend.position = 'none',
          plot.title.position = "plot", 
          plot.title = element_markdown(size =   6),
          axis.title.x = element_markdown(size = 6), 
          axis.title.y = element_markdown(size = 6), 
          axis.text = element_text(size = 6),
          plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), 'cm')
    ) +
    ggtitle(panel_titles_stde[i])
}

# rnd_normal
ROC_plots_rnd_normal <- list()
for(i in seq_along(roc_ref_prev)){
  ROC_plots_rnd_normal[[i]] <- ggroc(list(roc_ref_prev[[i]], roc_new_rnd_normal[[i]]),
                                     legacy.axes = TRUE) +
    scale_colour_manual(values = c('black', '#D51424')) +
    geom_abline(intercept = 0, slope = 1, color = 'grey45') +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(aspect.ratio = 1,
          legend.position = 'none',
          plot.title.position = "plot", 
          plot.title = element_markdown(size =   6),
          axis.title.x = element_markdown(size = 6), 
          axis.title.y = element_markdown(size = 6), 
          axis.text = element_text(size = 6),
          plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), 'cm')
    ) +
    ggtitle(panel_titles_rnd_normal[i])
}

# strat_rnd_normal
ROC_plots_strat_rnd_normal <- list()
for(i in seq_along(roc_ref_prev)){
  ROC_plots_strat_rnd_normal[[i]] <- ggroc(list(roc_ref_prev[[i]], roc_new_strat_rnd_normal[[i]]),
                                           legacy.axes = TRUE) +
    scale_colour_manual(values = c('black', '#D51424')) +
    geom_abline(intercept = 0, slope = 1, color = 'grey45') +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(aspect.ratio = 1,
          legend.position = 'none',
          plot.title.position = "plot", 
          plot.title = element_markdown(size =   6),
          axis.title.x = element_markdown(size = 6), 
          axis.title.y = element_markdown(size = 6), 
          axis.text = element_text(size = 6) ,
          plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), 'cm')
    ) +
    ggtitle(panel_titles_strat_rnd_normal[i])
}
# 
# ROC legend and caption ####
# legend
ROC_gg <- ggroc(list(`Reference model` = roc_ref_prev[[1]], `New model` = roc_new_glu[[1]])) +
  scale_colour_manual(values = c('black', '#D51424')) +
  theme(legend.title = element_blank()) +
  guides(colour = guide_legend(nrow = 1))

tmp <- ggplot_gtable(ggplot_build(ROC_gg))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
ROC_legend <- tmp$grobs[[leg]]

# caption
roc_caption <- expression(paste("*, P value < 0.05 of the likelihood-ratio test; " ^"#" *", P value < 0.05 of the DeLong's test for two correlated ROC curves (training dataset)"))
#
# plot all ROC curves on a single plot ####
tiff('results/ROC_curves_train.tiff', res = 300, width = 0.95*2250, height = 0.95*1590, units = 'px')
grid.arrange(ROC_plots_stde[[1]],
             ROC_plots_stde[[2]],
             ROC_plots_stde[[3]],
             ROC_plots_stde[[4]],
             ROC_plots_stde[[5]],
             ROC_plots_stde[[6]],
             ROC_plots_stde[[7]],
             
             ROC_plots_strat_rnd_normal[[1]],
             ROC_plots_strat_rnd_normal[[2]],
             ROC_plots_strat_rnd_normal[[3]],
             ROC_plots_strat_rnd_normal[[4]],
             ROC_plots_strat_rnd_normal[[5]],
             ROC_plots_strat_rnd_normal[[6]],
             ROC_plots_strat_rnd_normal[[7]],
             
             ROC_plots_glu[[1]],
             ROC_plots_glu[[2]],
             ROC_plots_glu[[3]],
             ROC_plots_glu[[4]],
             ROC_plots_glu[[5]],
             ROC_plots_glu[[6]],
             ROC_plots_glu[[7]],
             
             ROC_plots_rnd_normal[[1]],
             ROC_plots_rnd_normal[[2]],
             ROC_plots_rnd_normal[[3]],
             ROC_plots_rnd_normal[[4]],
             ROC_plots_rnd_normal[[5]],
             ROC_plots_rnd_normal[[6]],
             ROC_plots_rnd_normal[[7]],
             
             ROC_legend,
             
             layout_matrix = rbind(c(1:7),
                                   c(8:14),
                                   c(15:21),
                                   c(22:28),
                                   rep(29, 7)
             ) ,
             widths = rep(1, 7),
             heights = c(1, 1, 1, 1, 0.3),
             
             top = textGrob("ROC curves (training dataset)", gp = gpar(fontsize = 12),
                            hjust = 0, x = 0),
             
             bottom = textGrob(roc_caption,
                               gp = gpar(fontsize = 7),
                               hjust = 1, x = 1)
             
)
dev.off()
#
# TEST DATASETS ####
ind_prev01_test <- raw_indices %>% filter(dataset == 'test', prevalence == 0.01, iteration == iter_train[1])
ind_prev10_test <- raw_indices %>% filter(dataset == 'test', prevalence == 0.10, iteration == iter_train[2])
ind_prev30_test <- raw_indices %>% filter(dataset == 'test', prevalence == 0.30, iteration == iter_train[3])
ind_prev50_test <- raw_indices %>% filter(dataset == 'test', prevalence == 0.50, iteration == iter_train[4])
ind_prev70_test <- raw_indices %>% filter(dataset == 'test', prevalence == 0.70, iteration == iter_train[5])
ind_prev90_test <- raw_indices %>% filter(dataset == 'test', prevalence == 0.90, iteration == iter_train[6])
ind_prev99_test <- raw_indices %>% filter(dataset == 'test', prevalence == 0.99, iteration == iter_train[7])

ind_nonev_prev01_test <- ind_prev01_test[ind_prev01_test$class == 0,'index']
ind_nonev_prev10_test <- ind_prev10_test[ind_prev10_test$class == 0,'index']
ind_nonev_prev30_test <- ind_prev30_test[ind_prev30_test$class == 0,'index']
ind_nonev_prev50_test <- ind_prev50_test[ind_prev50_test$class == 0,'index']
ind_nonev_prev70_test <- ind_prev70_test[ind_prev70_test$class == 0,'index']
ind_nonev_prev90_test <- ind_prev90_test[ind_prev90_test$class == 0,'index']
ind_nonev_prev99_test <- ind_prev99_test[ind_prev99_test$class == 0,'index']

ind_event_prev01_test <- ind_prev01_test[ind_prev01_test$class == 1,'index']
ind_event_prev10_test <- ind_prev10_test[ind_prev10_test$class == 1,'index']
ind_event_prev30_test <- ind_prev30_test[ind_prev30_test$class == 1,'index']
ind_event_prev50_test <- ind_prev50_test[ind_prev50_test$class == 1,'index']
ind_event_prev70_test <- ind_prev70_test[ind_prev70_test$class == 1,'index']
ind_event_prev90_test <- ind_prev90_test[ind_prev90_test$class == 1,'index']
ind_event_prev99_test <- ind_prev99_test[ind_prev99_test$class == 1,'index']

ds_prev01_test <- rbind(ds_nonev[ind_nonev_prev01_test,], ds_event[ind_event_prev01_test,])
ds_prev10_test <- rbind(ds_nonev[ind_nonev_prev10_test,], ds_event[ind_event_prev10_test,])
ds_prev30_test <- rbind(ds_nonev[ind_nonev_prev30_test,], ds_event[ind_event_prev30_test,])
ds_prev50_test <- rbind(ds_nonev[ind_nonev_prev50_test,], ds_event[ind_event_prev50_test,])
ds_prev70_test <- rbind(ds_nonev[ind_nonev_prev70_test,], ds_event[ind_event_prev70_test,])
ds_prev90_test <- rbind(ds_nonev[ind_nonev_prev90_test,], ds_event[ind_event_prev90_test,])
ds_prev99_test <- rbind(ds_nonev[ind_nonev_prev99_test,], ds_event[ind_event_prev99_test,])

ds_prev01_test$prevalence <- rep(0.01, 100)
ds_prev10_test$prevalence <- rep(0.10, 100)
ds_prev30_test$prevalence <- rep(0.30, 100)
ds_prev50_test$prevalence <- rep(0.50, 100)
ds_prev70_test$prevalence <- rep(0.70, 100)
ds_prev90_test$prevalence <- rep(0.90, 100)
ds_prev99_test$prevalence <- rep(0.99, 100)

dataset_prev_test <- rbind(ds_prev01_test,
                           ds_prev10_test,
                           ds_prev30_test,
                           ds_prev50_test,
                           ds_prev70_test,
                           ds_prev90_test,
                           ds_prev99_test)
## write_xlsx(dataset_prev_test, 'results/dataset_prev_test.xlsx')
# 
# TEST: build reference models for each prevalence: ref_models_prev_test ####
ref_models_prev_test <- list()
for(i in seq_along(prevalence)){
  ref_models_prev_test[[i]] <- predict(ref_models_prev[[i]],
                                       dataset_prev_test[dataset_prev_test$prevalence == prevalence[i],],
                                       type = 'r')
}
names(ref_models_prev_test) <- labels_prevalence
#
# TEST: ROC curves for ref models: roc_ref_prev_test ####
roc_ref_prev_test <- list()
for(i in seq_along(prevalence)){
  roc_ref_prev_test[[i]] <- roc(dataset_prev_test[dataset_prev_test$prevalence == prevalence[i],]$disease, 
                                ref_models_prev_test[[i]], dir = '<')
}
#
# TEST: build new models for each prevalence: new_models_XXX_test ####
# glucose
new_models_glu_test <- list()
for(i in seq_along(prevalence)){
  new_models_glu_test[[i]] <- predict(new_models_glu[[i]], 
                                      dataset_prev_test[dataset_prev_test$prevalence == prevalence[i],],
                                      type = 'r') 
}
names(new_models_glu_test) <- labels_prevalence

# stde
new_models_stde_test <- list()
for(i in seq_along(prevalence)){
  new_models_stde_test[[i]] <- predict(new_models_stde[[i]], 
                                       dataset_prev_test[dataset_prev_test$prevalence == prevalence[i],],
                                       type = 'r')
}
names(new_models_stde_test) <- labels_prevalence

# rnd_normal
new_models_rnd_normal_test <- list()
for(i in seq_along(prevalence)){
  new_models_rnd_normal_test[[i]] <- predict(new_models_rnd_normal[[i]], 
                                             dataset_prev_test[dataset_prev_test$prevalence == prevalence[i],],
                                             type = 'r')
}
names(new_models_rnd_normal_test) <- labels_prevalence

# strat_rnd_normal
new_models_strat_rnd_normal_test <- list()
for(i in seq_along(prevalence)){
  new_models_strat_rnd_normal_test[[i]] <- predict(new_models_strat_rnd_normal[[i]], 
                                                   dataset_prev_test[dataset_prev_test$prevalence == prevalence[i],],
                                                   type = 'r')}
names(new_models_strat_rnd_normal_test) <- labels_prevalence
# 
# TEST: ROC curves for new models: roc_new_XXX_test ####
roc_new_glu_test <- list()
for(i in seq_along(prevalence)){
  roc_new_glu_test[[i]] <- roc(dataset_prev_test[dataset_prev_test$prevalence == prevalence[i],]$disease,
                               new_models_glu_test[[i]], dir = '<')
}

roc_new_stde_test <- list()
for(i in seq_along(prevalence)){
  roc_new_stde_test[[i]] <- roc(dataset_prev_test[dataset_prev_test$prevalence == prevalence[i],]$disease,
                                new_models_stde_test[[i]], dir = '<')
}

roc_new_rnd_normal_test <- list()
for(i in seq_along(prevalence)){
  roc_new_rnd_normal_test[[i]] <- roc(dataset_prev_test[dataset_prev_test$prevalence == prevalence[i],]$disease,
                                      new_models_rnd_normal_test[[i]], dir = '<')
}

roc_new_strat_rnd_normal_test <- list()
for(i in seq_along(prevalence)){
  roc_new_strat_rnd_normal_test[[i]] <- roc(dataset_prev_test[dataset_prev_test$prevalence == prevalence[i],]$disease,
                                            new_models_strat_rnd_normal_test[[i]], dir = '<')
}
#
# TEST: plots of ROC curves ####
# glucose
ROC_plots_glu_test <- list()
for(i in seq_along(roc_ref_prev_test)){
  ROC_plots_glu_test[[i]] <- ggroc(list(roc_ref_prev_test[[i]], roc_new_glu_test[[i]]),
                                   legacy.axes = TRUE) +
    scale_colour_manual(values = c('black', '#D51424')) +
    geom_abline(intercept = 0, slope = 1, color = 'grey45') +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(aspect.ratio = 1,
          legend.position = 'none',
          plot.title.position = "plot", 
          plot.title = element_markdown(size =   6),
          axis.title.x = element_markdown(size = 6), 
          axis.title.y = element_markdown(size = 6), 
          axis.text = element_text(size = 6),
          plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), 'cm')
    ) +
    ggtitle(panel_titles_glu[i])
  
}

# stde
ROC_plots_stde_test <- list()
for(i in seq_along(roc_ref_prev_test)){
  ROC_plots_stde_test[[i]] <- ggroc(list(roc_ref_prev_test[[i]], roc_new_stde_test[[i]]),
                                    legacy.axes = TRUE) +
    scale_colour_manual(values = c('black', '#D51424')) +
    geom_abline(intercept = 0, slope = 1, color = 'grey45') +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(aspect.ratio = 1,
          legend.position = 'none',
          plot.title.position = "plot", 
          plot.title = element_markdown(size =   6),
          axis.title.x = element_markdown(size = 6), 
          axis.title.y = element_markdown(size = 6), 
          axis.text = element_text(size = 6),
          plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), 'cm')
    ) +
    ggtitle(panel_titles_stde[i])
}

# rnd_normal
ROC_plots_rnd_normal_test <- list()
for(i in seq_along(roc_ref_prev_test)){
  ROC_plots_rnd_normal_test[[i]] <- ggroc(list(roc_ref_prev_test[[i]], roc_new_rnd_normal_test[[i]]),
                                          legacy.axes = TRUE) +
    scale_colour_manual(values = c('black', '#D51424')) +
    geom_abline(intercept = 0, slope = 1, color = 'grey45') +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(aspect.ratio = 1,
          legend.position = 'none',
          plot.title.position = "plot", 
          plot.title = element_markdown(size =   6),
          axis.title.x = element_markdown(size = 6), 
          axis.title.y = element_markdown(size = 6), 
          axis.text = element_text(size = 6),
          plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), 'cm')
    ) +
    ggtitle(panel_titles_rnd_normal[i])
}

# strat_rnd_normal
ROC_plots_strat_rnd_normal_test <- list()
for(i in seq_along(roc_ref_prev_test)){
  ROC_plots_strat_rnd_normal_test[[i]] <- ggroc(list(roc_ref_prev_test[[i]], roc_new_strat_rnd_normal_test[[i]]),
                                                legacy.axes = TRUE) +
    scale_colour_manual(values = c('black', '#D51424')) +
    geom_abline(intercept = 0, slope = 1, color = 'grey45') +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(aspect.ratio = 1,
          legend.position = 'none',
          plot.title.position = "plot", 
          plot.title = element_markdown(size =   6),
          axis.title.x = element_markdown(size = 6), 
          axis.title.y = element_markdown(size = 6), 
          axis.text = element_text(size = 6) ,
          plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), 'cm')
    ) +
    ggtitle(panel_titles_strat_rnd_normal[i])
}
# 
# TEST: plot all ROC curves on a single plot ####

tiff('results/ROC_curves_test.tiff', res = 300, width = 0.95*2250, height = 0.95*1590, units = 'px')
grid.arrange(ROC_plots_stde_test[[1]],
             ROC_plots_stde_test[[2]],
             ROC_plots_stde_test[[3]],
             ROC_plots_stde_test[[4]],
             ROC_plots_stde_test[[5]],
             ROC_plots_stde_test[[6]],
             ROC_plots_stde_test[[7]],
             
             ROC_plots_strat_rnd_normal_test[[1]],
             ROC_plots_strat_rnd_normal_test[[2]],
             ROC_plots_strat_rnd_normal_test[[3]],
             ROC_plots_strat_rnd_normal_test[[4]],
             ROC_plots_strat_rnd_normal_test[[5]],
             ROC_plots_strat_rnd_normal_test[[6]],
             ROC_plots_strat_rnd_normal_test[[7]],
             
             ROC_plots_glu_test[[1]],
             ROC_plots_glu_test[[2]],
             ROC_plots_glu_test[[3]],
             ROC_plots_glu_test[[4]],
             ROC_plots_glu_test[[5]],
             ROC_plots_glu_test[[6]],
             ROC_plots_glu_test[[7]],
             
             ROC_plots_rnd_normal_test[[1]],
             ROC_plots_rnd_normal_test[[2]],
             ROC_plots_rnd_normal_test[[3]],
             ROC_plots_rnd_normal_test[[4]],
             ROC_plots_rnd_normal_test[[5]],
             ROC_plots_rnd_normal_test[[6]],
             ROC_plots_rnd_normal_test[[7]],
             
             ROC_legend,
             
             layout_matrix = rbind(c(1:7),
                                   c(8:14),
                                   c(15:21),
                                   c(22:28),
                                   rep(29, 7)
             ) ,
             widths = rep(1, 7),
             heights = c(1, 1, 1, 1, 0.3),
             
             top = textGrob("ROC curves (test dataset)", gp = gpar(fontsize = 12),
                            hjust = 0, x = 0),
             
             bottom = textGrob(roc_caption,
                               gp = gpar(fontsize = 7),
                               hjust = 1, x = 1)
             
)
dev.off()
#
