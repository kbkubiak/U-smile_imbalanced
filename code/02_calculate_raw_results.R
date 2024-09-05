source('code/00_settings.R')

# calculate raw results #### 
set.seed(1034765298)

raw_results_01 <- calculate.raw.results.w.test(prev = 0.01, rep = 1000, ds_nonev = ds_nonev, ds_event = ds_event, ref_vars = ref_vars, new_vars = new_vars)
raw_results_10 <- calculate.raw.results.w.test(prev = 0.10, rep = 1000, ds_nonev = ds_nonev, ds_event = ds_event, ref_vars = ref_vars, new_vars = new_vars)
raw_results_30 <- calculate.raw.results.w.test(prev = 0.30, rep = 1000, ds_nonev = ds_nonev, ds_event = ds_event, ref_vars = ref_vars, new_vars = new_vars)
raw_results_50 <- calculate.raw.results.w.test(prev = 0.50, rep = 1000, ds_nonev = ds_nonev, ds_event = ds_event, ref_vars = ref_vars, new_vars = new_vars)
raw_results_70 <- calculate.raw.results.w.test(prev = 0.70, rep = 1000, ds_nonev = ds_nonev, ds_event = ds_event, ref_vars = ref_vars, new_vars = new_vars)
raw_results_90 <- calculate.raw.results.w.test(prev = 0.90, rep = 1000, ds_nonev = ds_nonev, ds_event = ds_event, ref_vars = ref_vars, new_vars = new_vars)
raw_results_99 <- calculate.raw.results.w.test(prev = 0.99, rep = 1000, ds_nonev = ds_nonev, ds_event = ds_event, ref_vars = ref_vars, new_vars = new_vars)

raw_results_ref <- rbind(raw_results_01[[1]],
                         raw_results_10[[1]],
                         raw_results_30[[1]],
                         raw_results_50[[1]],
                         raw_results_70[[1]],
                         raw_results_90[[1]],
                         raw_results_99[[1]])

raw_results_comp <- rbind(raw_results_01[[2]],
                          raw_results_10[[2]],
                          raw_results_30[[2]],
                          raw_results_50[[2]],
                          raw_results_70[[2]],
                          raw_results_90[[2]],
                          raw_results_99[[2]])

raw_convergence <- rbind(raw_results_01[[3]],
                         raw_results_10[[3]],
                         raw_results_30[[3]],
                         raw_results_50[[3]],
                         raw_results_70[[3]],
                         raw_results_90[[3]],
                         raw_results_99[[3]])

raw_indices <- rbind(raw_results_01[[4]],
                     raw_results_10[[4]],
                     raw_results_30[[4]],
                     raw_results_50[[4]],
                     raw_results_70[[4]],
                     raw_results_90[[4]],
                     raw_results_99[[4]])

# write_xlsx(raw_results_ref,   'results/raw_results_ref.xlsx')
# write_xlsx(raw_results_comp,  'results/raw_results_comp.xlsx')
# write_xlsx(raw_convergence,   'results/raw_convergence.xlsx')
# write.table(raw_indices,      'results/raw_indices.txt')
#
# identify converged ref models: raw_results_ref_conv ####
not_converged_ref <- raw_convergence %>% filter(conv_ref == F) %>% select(iteration, prevalence)
not_converged_ref_iter_prev <- unique(not_converged_ref)
not_converged_ref_iter_prev
# 2 reference models at prev 1% and 12 ref models at prev 99% 
## iteration prevalence
##     623       0.01
##     964       0.01
##     131       0.99
##     154       0.99
##     263       0.99
##     435       0.99
##     437       0.99
##     596       0.99
##     709       0.99
##     776       0.99
##     780       0.99
##     883       0.99
##     934       0.99
##    1000       0.99

ind_not_converged_ref <- c(
  which(raw_results_ref$iteration %in% c(623, 964)  & raw_results_ref$prevalence ==  0.01),
  which(raw_results_ref$iteration %in% c( 131,
                                          154,
                                          263,
                                          435,
                                          437,
                                          596,
                                          709,
                                          776,
                                          780,
                                          883,
                                          934,
                                          1000)    & raw_results_ref$prevalence ==  0.99) )

raw_results_ref_conv <- raw_results_ref[-ind_not_converged_ref,]
raw_results_ref_conv$dataset <- factor(raw_results_ref_conv$dataset, levels = c('train', 'test'))
## write_xlsx(raw_results_ref_conv, 'results/raw_results_ref_conv.xlsx')
#
# identify converged new models: raw_results_comp_conv ####
raw_convergence %>% 
  filter(conv_ref == T & conv == T) %>% 
  group_by(prevalence, model) %>% 
  summarise(n = n()) -> how_many_models_converged

# write_xlsx(how_many_models_converged, 'results/how_many_models_converged.xlsx')

# only 1 and 99% prevalences did not converge
# identify iterations that converged for 1 and 99%

conv_iter_01_glu                <- subset(raw_convergence, prevalence == 0.01 & conv_ref == T & conv == T & model == 'glu' )$iteration
conv_iter_01_stde               <- subset(raw_convergence, prevalence == 0.01 & conv_ref == T & conv == T & model == 'stde' )$iteration
conv_iter_01_rnd_normal         <- subset(raw_convergence, prevalence == 0.01 & conv_ref == T & conv == T & model == 'rnd_normal' )$iteration
conv_iter_01_strat_rnd_normal   <- subset(raw_convergence, prevalence == 0.01 & conv_ref == T & conv == T & model == 'strat_rnd_normal' )$iteration

conv_iter_99_glu                <- subset(raw_convergence, prevalence == 0.99 & conv_ref == T & conv == T & model == 'glu' )$iteration 
conv_iter_99_stde               <- subset(raw_convergence, prevalence == 0.99 & conv_ref == T & conv == T & model == 'stde' )$iteration
conv_iter_99_rnd_normal         <- subset(raw_convergence, prevalence == 0.99 & conv_ref == T & conv == T & model == 'rnd_normal' )$iteration
conv_iter_99_strat_rnd_normal   <- subset(raw_convergence, prevalence == 0.99 & conv_ref == T & conv == T & model == 'strat_rnd_normal' )$iteration

raw_results_comp_conv <- rbind(raw_results_comp %>% filter(prevalence == 0.01 & model == 'glu'                 & iteration %in% conv_iter_01_glu),
                               raw_results_comp %>% filter(prevalence == 0.01 & model == 'stde'                & iteration %in% conv_iter_01_stde),
                               raw_results_comp %>% filter(prevalence == 0.01 & model == 'rnd_normal'          & iteration %in% conv_iter_01_rnd_normal),
                               raw_results_comp %>% filter(prevalence == 0.01 & model == 'strat_rnd_normal'    & iteration %in% conv_iter_01_strat_rnd_normal),
                               raw_results_comp %>% filter(prevalence %in% c(0.10, 0.30, 0.50, 0.70, 0.90)),
                               raw_results_comp %>% filter(prevalence == 0.99 & model == 'glu'                 & iteration %in% conv_iter_99_glu),
                               raw_results_comp %>% filter(prevalence == 0.99 & model == 'stde'                & iteration %in% conv_iter_99_stde),
                               raw_results_comp %>% filter(prevalence == 0.99 & model == 'rnd_normal'          & iteration %in% conv_iter_99_rnd_normal),
                               raw_results_comp %>% filter(prevalence == 0.99 & model == 'strat_rnd_normal'    & iteration %in% conv_iter_99_strat_rnd_normal) 
)

raw_results_comp_conv$model   <- factor(raw_results_comp_conv$model,   levels = new_vars)
raw_results_comp_conv$dataset <- factor(raw_results_comp_conv$dataset, levels = c('train', 'test'))
# write_xlsx(raw_results_comp_conv, 'results/raw_results_comp_conv.xlsx')
#
