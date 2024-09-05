source('code/00_settings.R')

# import raw results ####
raw_results_ref_conv  <- read_excel('results/raw_results_ref_conv.xlsx')
raw_results_ref_conv$dataset <- factor(raw_results_ref_conv$dataset, levels = c('train', 'test'))
raw_results_comp_conv <- read_excel('results/raw_results_comp_conv.xlsx')
raw_results_comp_conv$dataset <- factor(raw_results_comp_conv$dataset, levels = c('train', 'test'))
raw_results_comp_conv$model <- factor(raw_results_comp_conv$model, levels = new_vars)
#
which(is.na(raw_results_comp_conv$RB_overall))
ind <- c(48620, 49116, 49702)
kk <- raw_results_comp_conv[ind,]
kk %>% select(iteration, dataset, prevalence, model, starts_with('RB'))
iter_glu <- c(389, 641, 939)
# calculate NRI ####
raw_results_comp_conv$NRI <- raw_results_comp_conv$II_nonev + raw_results_comp_conv$II_event
#
# calculate summaries - means - ref model: results_ref_mean ####
raw_results_ref_conv %>% group_by(dataset, prevalence) %>% 
  summarise(AUC_ref  = mean(AUC_ref,  na.rm = T),
            BS_ref   = mean(BS_ref,   na.rm = T),
            BS_0_ref = mean(BS_0_ref, na.rm = T),
            BS_1_ref = mean(BS_1_ref, na.rm = T),
            MCC_ref  = mean(MCC_ref,  na.rm = T),
            F1_ref   = mean(F1_ref,   na.rm = T),
            TN_ref   = mean(TN_ref,   na.rm = T),
            FP_ref   = mean(FP_ref,   na.rm = T),
            FN_ref   = mean(FN_ref,   na.rm = T),
            TP_ref   = mean(TP_ref,   na.rm = T) ) %>% 
  mutate(desc_stat = 'mean') %>% relocate(desc_stat, .after = prevalence) -> results_ref_mean
write_xlsx(results_ref_mean, 'results/results_means_ref_models.xlsx')
#
# plot of strat BS and difference in strat BS between training and test ds ####
results_ref_mean %>% group_by(dataset) %>% 
  pivot_longer(cols = c(BS_0_ref, BS_1_ref), names_to = 'BS', values_to = 'value') -> str_bs
str_bs$dataset <- factor(str_bs$dataset, levels = c('train', 'test'))

ggplot(str_bs, aes(x = prevalence * 100, y = value, color = BS)) +
  geom_point(size = 1) +
  geom_smooth(se=F, size = 0.5, aes(linetype = dataset)) +
  scale_color_manual(values = c('#0F3C78',
                                '#D51424'),
                     labels = c('non-events', 'events')) +
  scale_linetype_manual(values = c('solid', 'dashed'),
                        labels = c('training', 'test')) +
  scale_x_continuous(breaks = prevalence_breaks) +
  guides(linetype = guide_legend(nrow=2, override.aes = list(color = 'black')),
         color = guide_legend(nrow=2)) +
  xlab('Percentage of the event class') +
  ylab('Stratified BS') +
  labs(color = 'class') +
  ylim(c(-0.02, 1)) +
  theme(plot.title.position = "plot", 
        legend.position = 'bottom',
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        plot.title = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        #legend.spacing.x = unit(0, "cm"),
        legend.key.width  = unit(0.5, 'cm'),
        legend.key.height = unit(0.2, 'cm'),
        aspect.ratio = 1) +
  ggtitle('A Stratified BS') -> plot_str_BS
#
# plot of the difference in strat BS between training and test datasets 
raw_results_ref_conv %>% filter(dataset == 'train') %>% 
  select(prevalence, dataset, BS_0_ref, BS_1_ref) %>% 
  mutate(prevalence = prevalence, 
         BS_0_ref_train = BS_0_ref,
         BS_1_ref_train = BS_1_ref) %>% 
  select(prevalence,
         BS_0_ref_train,
         BS_1_ref_train) -> raw_BS_train
raw_results_ref_conv %>% filter(dataset == 'test') %>% 
  select(prevalence, dataset, BS_0_ref, BS_1_ref) %>% 
  mutate(prevalence = prevalence, 
         BS_0_ref_test = BS_0_ref,
         BS_1_ref_test = BS_1_ref) %>% 
  select(BS_0_ref_test,
         BS_1_ref_test) -> raw_BS_test

raw_BS_diff <- cbind(raw_BS_train, raw_BS_test)
BS_diff <- raw_BS_diff %>% 
  mutate(diff_BS_0 = BS_0_ref_test - BS_0_ref_train,
         diff_BS_1 = BS_1_ref_test - BS_1_ref_train) %>% 
  group_by(prevalence) %>% 
  summarise(diff_BS_0 = mean(diff_BS_0),
            diff_BS_1 = mean(diff_BS_1)) %>% 
  pivot_longer(cols = c(diff_BS_0, diff_BS_1), names_to = 'diff_BS', values_to = 'value')

##
ggplot(BS_diff, aes(x = prevalence * 100, y = value,  color = diff_BS) ) +
  geom_point(size = 1) +
  geom_smooth(se=F, size = 0.5, span = 0.75) +
  scale_color_manual(values = c('#0F3C78',
                                '#D51424'),
                     labels = c('non-events', 'events')) +
  scale_x_continuous(breaks = prevalence_breaks) +
  guides(color = guide_legend(nrow=2)) +
  xlab('Percentage of the event class') +
  ylab('Stratified BS difference') +
  labs(color = 'class') +
  #ylim(c(-0.02, 1)) +
  theme(plot.title.position = "plot", 
        legend.position = 'bottom',
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        plot.title = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.width  = unit(0.5, 'cm'),
        legend.key.height = unit(0.2, 'cm'),
        aspect.ratio = 1) +
  ggtitle('B Stratified BS difference (test vs. training dataset)') -> plot_diff_BS
#
tiff('results/ref_BS_plot_smoother.tiff', res = 300, width = 0.95 * 1855, height = 950, units = 'px')
grid.arrange(plot_str_BS, plot_diff_BS, nrow = 1)
dev.off()
#
# calculate summaries - means - new models: results_new_mean ####
raw_results_comp_conv %>% group_by(dataset, prevalence, model) %>% 
  summarise(LRT_p_value = mean(LRT_p_value, na.rm = T),
            LRT_statistic = mean(LRT_statistic, na.rm = T),
            AUC_p_value = mean(AUC_p_value, na.rm = T),
            AUC_statistic = mean(AUC_statistic, na.rm = T),
            AUC         = mean(AUC,         na.rm = T),
            delta_AUC   = mean(delta_AUC,   na.rm = T),
            BS          = mean(BS,          na.rm = T),
            BS_0        = mean(BS_0,        na.rm = T),
            BS_1        = mean(BS_1,        na.rm = T),
            delta_BS    = mean(delta_BS,    na.rm = T),
            BSS         = mean(BSS,         na.rm = T),
            MCC         = mean(MCC,         na.rm = T),
            F1          = mean(F1,          na.rm = T),
            delta_MCC   = mean(delta_MCC,   na.rm = T), 
            delta_F1    = mean(delta_F1,    na.rm = T),
            NRI         = mean(NRI,         na.rm = T),
            TN          = mean(TN,          na.rm = T),
            FP          = mean(FP,          na.rm = T),
            FN          = mean(FN,          na.rm = T),
            TP          = mean(TP,          na.rm = T)      ) %>% 
  mutate(desc_stat = 'mean') %>% relocate(desc_stat, .after = model) -> results_new_mean
write_xlsx(results_new_mean, 'results/results_means_new_models.xlsx')
#
# MCC - for how many models MCC was not computed ####
raw_results_ref_conv %>% group_by(dataset, prevalence) %>% 
  summarise(n_missing_MCC_ref = sum(is.na(MCC_ref)) ) -> how_many_missing_MCC_ref

raw_results_comp_conv %>% group_by(dataset, prevalence, model) %>% 
  summarise(n_missing_MCC       = sum(is.na(MCC)),
            n_missing_delta_MCC = sum(is.na(delta_MCC))  ) -> how_many_missing_deltaMCC
write_xlsx(how_many_missing_MCC_ref , 'results/how_many_missing_MCC_ref.xlsx')
write_xlsx(how_many_missing_deltaMCC, 'results/how_many_missing_deltaMCC.xlsx')
# 
# LRT and delta_AUC: p values from the mean value of the test statistic ####
results_new_mean %>% filter(dataset == 'train') %>% 
  arrange(model) %>% 
  mutate(LRT_p = round(2 * (1 - pchisq(LRT_statistic, 1)), 3),
         AUC_p = round(2 * pnorm(AUC_statistic), 3 ) ) %>% 
  select(dataset, prevalence, model, LRT_p_value, AUC_p_value) -> p_values
write_xlsx(p_values, 'results/p_values_LRT_DeLong.xlsx')
#
# calculate BARBI ####
raw_results_comp_conv %>% group_by(dataset, prevalence, model) %>% 
  summarise(BA_nonev_be = mean(BA_nonev_be, na.rm = T),
            BA_nonev_wo = mean(BA_nonev_wo, na.rm = T),
            BA_event_wo = mean(BA_event_wo, na.rm = T),
            BA_event_be = mean(BA_event_be, na.rm = T),
            BA_nonev    = mean(BA_nonev   , na.rm = T),
            BA_event    = mean(BA_event   , na.rm = T),
            BA          = mean(BA_overall , na.rm = T),
            RB_nonev_be = mean(RB_nonev_be, na.rm = T),
            RB_nonev_wo = mean(RB_nonev_wo, na.rm = T),
            RB_event_wo = mean(RB_event_wo, na.rm = T),
            RB_event_be = mean(RB_event_be, na.rm = T),
            RB_nonev    = mean(RB_nonev   , na.rm = T),
            RB_event    = mean(RB_event   , na.rm = T),
            RB          = mean(RB_overall , na.rm = T),
            I_nonev_be  = mean(II_nonev_be, na.rm = T),
            I_nonev_wo  = mean(II_nonev_wo, na.rm = T),
            I_event_wo  = mean(II_event_wo, na.rm = T),
            I_event_be  = mean(II_event_be, na.rm = T),
            I_nonev     = mean(II_nonev   , na.rm = T),
            I_event     = mean(II_event   , na.rm = T),
            I           = mean(II_overall , na.rm = T) ) %>% 
  mutate(desc_stat = 'mean') %>% relocate(desc_stat, .after = model) -> BARBI_mean
write_xlsx(BARBI_mean, 'results/BARBI_mean.xlsx')
#
# prepare BARBIS for the U-smile plot ####
# nonev_be
BARBI_mean %>% group_by(dataset, prevalence, model) %>% 
  select(dataset, prevalence, model,
         BA_nonev_be,
         RB_nonev_be,
         I_nonev_be) %>% 
  summarise(BA = BA_nonev_be,
            RB = RB_nonev_be,
            I  = I_nonev_be,
            subclass = 'nonev_be') -> BARBI_nonev_be_for_plot
# nonev_wo
BARBI_mean %>% group_by(dataset, prevalence, model) %>% 
  select(dataset, prevalence, model,
         BA_nonev_wo,
         RB_nonev_wo,
         I_nonev_wo) %>% 
  summarise(BA = BA_nonev_wo,
            RB = RB_nonev_wo,
            I  =  I_nonev_wo,
            subclass =    'nonev_wo') -> BARBI_nonev_wo_for_plot
# event_wo
BARBI_mean %>% group_by(dataset, prevalence, model) %>% 
  select(dataset, prevalence, model,
         BA_event_wo,
         RB_event_wo,
         I_event_wo) %>% 
  summarise(BA = BA_event_wo,
            RB = RB_event_wo,
            I  =  I_event_wo,
            subclass =    'event_wo') -> BARBI_event_wo_for_plot
# event_be
BARBI_mean %>% group_by(dataset, prevalence, model) %>% 
  select(dataset, prevalence, model,
         BA_event_be,
         RB_event_be,
         I_event_be) %>% 
  summarise(BA = BA_event_be,
            RB = RB_event_be,
            I  =  I_event_be,
            subclass =    'event_be') -> BARBI_event_be_for_plot
# 
BARBI_for_plot <- rbind(BARBI_nonev_be_for_plot, BARBI_nonev_wo_for_plot, BARBI_event_wo_for_plot, BARBI_event_be_for_plot)
BARBI_for_plot$model <- factor(BARBI_for_plot$model, levels = new_vars)
BARBI_for_plot$prevalence <- factor(BARBI_for_plot$prevalence,
                                    labels = labels_prevalence)
BARBI_for_plot$model_pretty <- factor(BARBI_for_plot$model, levels = new_vars,
                                      labels = labels_model_pretty)
#
# U-smile plots ####
BARBI_for_plot %>% filter(dataset == 'train') %>% 
  ggplot(aes(x = subclass, y = BA, group = model)) +
  geom_point(aes(col = subclass, fill = subclass, size = I), shape = 21, stroke = 1) + 
  geom_line(color = 'gray30') +
  scale_x_discrete(limits = subclass_order,
                   expand = expansion(add = 1)) +
  scale_color_manual(values = usmile_colors,
                     breaks = subclass_order,
                     labels = usmile_labels) +
  scale_fill_manual(values = usmile_fills,
                    breaks = subclass_order,
                    labels = usmile_labels) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        #axis.text.y  = element_text(size = 9),
        axis.ticks.x = element_blank(),
        #axis.title.y = element_blank(),
        strip.text.x = element_text(size = 11), 
        strip.text.y = element_text(angle = -90, size = 7), 
        legend.position = 'none',
        plot.title.position = "plot") +
  ylim(c(-0.03,
         1.15*max(BARBI_for_plot$BA)) ) + 
  ggtitle('A  U-smile plots of the BA coefficient (training dataset)') +
  facet_grid(model_pretty ~ prevalence) -> plot_BA_train

BARBI_for_plot %>% filter(dataset == 'train') %>% 
  ggplot(aes(x = subclass, y = RB, group = model)) +
  geom_point(aes(col = subclass, fill = subclass, size = I), shape = 21, stroke = 1) + 
  geom_line(color = 'gray30') +
  scale_x_discrete(limits = subclass_order,
                   expand = expansion(add = 1)) +
  scale_color_manual(values = usmile_colors,
                     breaks = subclass_order,
                     labels = usmile_labels) +
  scale_fill_manual(values = usmile_fills,
                    breaks = subclass_order,
                    labels = usmile_labels) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_text(size = 9),
        axis.ticks.x = element_blank(),
        #axis.title.y = element_blank(),
        strip.text.x = element_text(size = 11), 
        strip.text.y = element_text(angle = -90, size = 7), 
        legend.position = 'bottom', legend.box="vertical", 
        plot.title.position = "plot") +
  guides(color = guide_legend(nrow = 2), 
         size  = guide_legend(nrow = 1)) + 
  coord_cartesian(ylim = c(-0.1, 1)) +
  ggtitle('B  U-smile plots of the RB coefficient (training dataset)') +
  facet_grid(model_pretty ~ prevalence) -> plot_RB_train

tiff('results/Usmile_train_grayline.tiff', res = 300, width = 0.95 * 1855, height = 0.95 * 2625, units = 'px')
plot(rbind(ggplotGrob(plot_BA_train),
           ggplotGrob(plot_RB_train), 
           size = 'first') )
dev.off()


BARBI_for_plot %>% filter(dataset == 'test') %>% 
  ggplot(aes(x = subclass, y = BA, group = model)) +
  geom_point(aes(col = subclass, fill = subclass, size = I), shape = 21, stroke = 1) + 
  geom_line(color = 'gray30') +
  scale_x_discrete(limits = subclass_order,
                   expand = expansion(add = 1)) +
  scale_color_manual(values = usmile_colors,
                     breaks = subclass_order,
                     labels = usmile_labels) +
  scale_fill_manual(values = usmile_fills,
                    breaks = subclass_order,
                    labels = usmile_labels) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        #axis.text.y  = element_text(size = 9),
        axis.ticks.x = element_blank(),
        #axis.title.y = element_blank(),
        strip.text.x = element_text(size = 11), 
        strip.text.y = element_text(angle = -90, size = 7), 
        legend.position = 'none',
        plot.title.position = "plot") +
  ylim(c(-0.03,
         1.15*max(BARBI_for_plot$BA)) ) + 
  ggtitle('A  U-smile plots of the BA coefficient (test dataset)') +
  facet_grid(model_pretty ~ prevalence) -> plot_BA_test

BARBI_for_plot %>% filter(dataset == 'test') %>% 
  ggplot(aes(x = subclass, y = RB, group = model)) +
  geom_point(aes(col = subclass, fill = subclass, size = I), shape = 21, stroke = 1) + 
  geom_line(color = 'gray30') +
  scale_x_discrete(limits = subclass_order,
                   expand = expansion(add = 1)) +
  scale_color_manual(values = usmile_colors,
                     breaks = subclass_order,
                     labels = usmile_labels) +
  scale_fill_manual(values = usmile_fills,
                    breaks = subclass_order,
                    labels = usmile_labels) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_text(size = 9),
        axis.ticks.x = element_blank(),
        #axis.title.y = element_blank(),
        strip.text.x = element_text(size = 11), 
        strip.text.y = element_text(angle = -90, size = 7), 
        legend.position = 'bottom', legend.box="vertical", 
        plot.title.position = "plot") +
  guides(color = guide_legend(nrow = 2), 
         size  = guide_legend(nrow = 1)) + 
  coord_cartesian(ylim = c(-0.1, 1)) +
  ggtitle('B  U-smile plots of the RB coefficient (test dataset)') +
  facet_grid(model_pretty ~ prevalence) -> plot_RB_test


tiff('results/Usmile_test_grayline.tiff', res = 300, width = 0.95 * 1855, height = 0.95 * 2625, units = 'px')
plot(rbind(ggplotGrob(plot_BA_test),
           ggplotGrob(plot_RB_test), 
           size = 'first') )
dev.off()
#

# Levels of BARBI - settings and functions to draw plots ####
BARBI_mean$prevalence100 <- BARBI_mean$prevalence * 100
levels_colors <- c('#15b922', 
                   '#15b922',
                   '#ab08fc', 
                   '#ab08fc'
                   
)
levels_linetypes <- c('solid', 'dashed', 'solid', 'dashed')
# extract legend
BARBI_mean %>% filter(dataset == 'train') %>% group_by(model) %>% 
  ggplot(aes(x = prevalence, y = BA_nonev_be, color = model)) +
  geom_smooth(aes(linetype = model), se = F) +
  scale_linetype_manual(values = levels_linetypes,
                        labels = labels_model_pretty) +
  scale_color_manual(values = levels_colors,
                     labels = labels_model_pretty) +
  guides(color = guide_legend(nrow=1),
         linetype = guide_legend(nrow=1)) +
  theme(legend.title = element_blank(),
        legend.key.width= unit(1, 'cm')) -> plot_lvl_legend
lvl_legend <- get_legend(plot_lvl_legend)
as_ggplot(lvl_legend)
#
# for upper plots
plot.level <- function(data, x, y, color, ymin, ymax, ylab, span = 0.75){ 
  ggplot(data, aes({{x}}, {{y}}, color = {{color}})) + 
    #ylim(c(ymin, ymax)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    geom_smooth(aes(linetype = {{color}}), se = F, size = 0.5, span = span) + 
    scale_linetype_manual(values = levels_linetypes ) +
    scale_color_manual(values = levels_colors) +
    scale_x_continuous(breaks = prevalence_breaks) +
    ylab(ylab) +
    theme(legend.position = "none",
          plot.margin = unit(c(0,0.1,0,0.1), 'cm'),   
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title.y = element_text(size = 8),
          axis.text.y  = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank() ) 
}
# for the lowest plot in group
plot.level.axis <- function(data, x, y, color, ymin, ymax, ylab, span = 0.75){ 
  ggplot(data, aes({{x}}, {{y}}, color = {{color}})) + 
    # ylim(c(ymin, ymax)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    geom_smooth(aes(linetype = {{color}}), se = F, size = 0.5, span = span) + 
    scale_linetype_manual(values = levels_linetypes ) +
    scale_color_manual(values = levels_colors) +
    scale_x_continuous(breaks = prevalence_breaks) +
    ylab(ylab) +
    theme(legend.position = "none",
          plot.margin = unit(c(0,0.1,0.2,0.1), 'cm'), 
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title.y = element_text(size = 8),
          axis.text.y  = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          axis.text.x  = element_text(size = 6)) +
    xlab('Percentage of the event class')
}
# for level 3 plots
plot.level3 <- function(data, x, y, color, ymin, ymax, ylab, span = 0.75){ 
  ggplot(data, aes({{x}}, {{y}}, color = {{color}})) + 
    #ylim(c(ymin, ymax)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    geom_smooth(aes(linetype = {{color}}), se = F, size = 0.5, span = span) + 
    scale_linetype_manual(values = levels_linetypes ) +
    scale_color_manual(values = levels_colors) +
    scale_x_continuous(breaks = prevalence_breaks) +
    ylab(ylab) +
    theme(legend.position = "none",
          plot.margin = unit(c(0,0.1,0.2,0.1), 'cm'), 
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title.y = element_text(size = 8),
          axis.text.y  = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          axis.text.x  = element_text(size = 6)) +
    xlab('Percentage of the event class')
}
#
# TRAIN:: levels plots ####
# lvl 1 BA  
min_ba_l1 <- -0.002
max_ba_l1 <- max(BARBI_mean$BA_nonev_be, BARBI_mean$BA_event_be)

p01 <-      plot.level(filter(BARBI_mean, dataset=='train'), prevalence100, BA_nonev_be, model, min_ba_l1, max_ba_l1, 'BA0+')
p02 <-      plot.level(filter(BARBI_mean, dataset=='train'), prevalence100, BA_nonev_wo, model, min_ba_l1, max_ba_l1, 'BA0-')
p03 <-      plot.level(filter(BARBI_mean, dataset=='train'), prevalence100, BA_event_wo, model, min_ba_l1, max_ba_l1, 'BA1-')
p04 <- plot.level.axis(filter(BARBI_mean, dataset=='train'), prevalence100, BA_event_be, model, min_ba_l1, max_ba_l1, 'BA1+')

# lvl 1 RB 
p05 <-      plot.level(filter(BARBI_mean, dataset=='train'), prevalence100, RB_nonev_be, model, 0, 1, 'RB0+')
p06 <-      plot.level(filter(BARBI_mean, dataset=='train'), prevalence100, RB_nonev_wo, model, 0, 1, 'RB0-', span = 0.9)
p07 <-      plot.level(filter(BARBI_mean, dataset=='train'), prevalence100, RB_event_wo, model, 0, 1, 'RB1-', span = 0.9)
p08 <- plot.level.axis(filter(BARBI_mean, dataset=='train'), prevalence100, RB_event_be, model, 0, 1, 'RB1+')

# lvl 1 I 
p09 <-      plot.level(filter(BARBI_mean, dataset=='train'), prevalence100, I_nonev_be, model, 0, 1, 'I0+')
p10 <-      plot.level(filter(BARBI_mean, dataset=='train'), prevalence100, I_nonev_wo, model, 0, 1, 'I0-')
p11 <-      plot.level(filter(BARBI_mean, dataset=='train'), prevalence100, I_event_wo, model, 0, 1, 'I1-')
p12 <- plot.level.axis(filter(BARBI_mean, dataset=='train'), prevalence100, I_event_be, model, 0, 1, 'I1+')

# lvl 2
min_ba_l2 <- -0.01
max_ba_l2 <- 1.05*max(BARBI_mean$BA_nonev, BARBI_mean$BA_event)

p13 <-      plot.level(filter(BARBI_mean, dataset=='train'), prevalence100, BA_nonev, model, min_ba_l2, max_ba_l2, 'BA0', span = 0.9)
p14 <- plot.level.axis(filter(BARBI_mean, dataset=='train'), prevalence100, BA_event, model, min_ba_l2, max_ba_l2, 'BA1', span = 0.9)

p15 <-      plot.level(filter(BARBI_mean, dataset=='train'), prevalence100, RB_nonev, model, -0.9, 0.5, 'RB0', span = 0.9)
p16 <- plot.level.axis(filter(BARBI_mean, dataset=='train'), prevalence100, RB_event, model, -0.9, 0.5, 'RB1', span = 0.9)

p17 <-      plot.level(filter(BARBI_mean, dataset=='train'), prevalence100, I_nonev, model, -0.3, 0.6, 'I0', span = 0.9)
p18 <- plot.level.axis(filter(BARBI_mean, dataset=='train'), prevalence100, I_event, model, -0.3, 0.6, 'I1', span = 0.9)

# lvl 3
p19 <- plot.level3(filter(BARBI_mean, dataset=='train'), prevalence100, BA, model, -0.002, 0.04, ylab = 'BA', span = 0.8)
p20 <- plot.level3(filter(BARBI_mean, dataset=='train'), prevalence100, RB, model, -1, 0.5, ylab = 'RB', span = 1.5)
p21 <- plot.level3(filter(BARBI_mean, dataset=='train'), prevalence100, I,  model, -0.2, 0.6, ylab = 'I',  span = 1)

# LEVELS final plot ####
txt_title_train <- textGrob("Trends in BA-RB-I coefficients (training dataset)", gp = gpar(fontisize = 12))
txt_lvl1 <- textGrob("Level 1", gp = gpar(fontisize = 10) ) #
txt_lvl2 <- textGrob("Level 2", gp = gpar(fontisize = 10) )
txt_lvl3 <- textGrob("Level 3", gp = gpar(fontisize = 10) )


tiff('results/LEVELS_trends_train.tiff', res = 300, width = 0.95 * 1855, height = 0.95 * 2625, units = 'px')
grid.arrange(txt_title_train, 
             txt_lvl1, txt_lvl2, txt_lvl3, 
             
             rbind(ggplotGrob(p01), 
                   ggplotGrob(p02), 
                   ggplotGrob(p03), 
                   ggplotGrob(p04),
                   ggplotGrob(p05), 
                   ggplotGrob(p06), 
                   ggplotGrob(p07), 
                   ggplotGrob(p08),
                   ggplotGrob(p09), 
                   ggplotGrob(p10), 
                   ggplotGrob(p11), 
                   ggplotGrob(p12), size = 'first'), 
             
             rbind(rbind(ggplotGrob(p13),
                         ggplotGrob(p14), 
                         size = 'first'), 
                   rbind(ggplotGrob(p15),
                         ggplotGrob(p16), 
                         size = 'first'),
                   rbind(ggplotGrob(p17),
                         ggplotGrob(p18), 
                         size = 'first') ), 
             
             rbind(ggplotGrob(p19), 
                   ggplotGrob(p20), 
                   ggplotGrob(p21)),
             
             lvl_legend, 
             
             #layout_matrix = rbind(c(1,2,3),
             #                      c(4,5,6),
             #                      c(7,7,7)),
             #
             #heights = c(0.4, 10, 0.4)
             
             layout_matrix = rbind(c(1,1,1),
                                   c(2,3,4),
                                   c(5,6,7),
                                   c(8,8,8)),
             heights = c(0.4, 0.3, 10, 0.4)
             
)
dev.off()
#
# 
# LEVELS FOR TEST DATASET ####
# lvl 1 BA  
p01_te <-      plot.level(filter(BARBI_mean, dataset=='test'), prevalence100, BA_nonev_be, model, min_ba_l1, max_ba_l1, 'BA0+')
p02_te <-      plot.level(filter(BARBI_mean, dataset=='test'), prevalence100, BA_nonev_wo, model, min_ba_l1, max_ba_l1, 'BA0-')
p03_te <-      plot.level(filter(BARBI_mean, dataset=='test'), prevalence100, BA_event_wo, model, min_ba_l1, max_ba_l1, 'BA1-')
p04_te <- plot.level.axis(filter(BARBI_mean, dataset=='test'), prevalence100, BA_event_be, model, min_ba_l1, max_ba_l1, 'BA1+')

# lvl 1 RB 
p05_te <-      plot.level(filter(BARBI_mean, dataset=='test'), prevalence100, RB_nonev_be, model, 0, 1, 'RB0+')
p06_te <-      plot.level(filter(BARBI_mean, dataset=='test'), prevalence100, RB_nonev_wo, model, 0, 1, 'RB0-', span = 0.9)
p07_te <-      plot.level(filter(BARBI_mean, dataset=='test'), prevalence100, RB_event_wo, model, 0, 1, 'RB1-', span = 0.9)
p08_te <- plot.level.axis(filter(BARBI_mean, dataset=='test'), prevalence100, RB_event_be, model, 0, 1, 'RB1+')

# lvl 1 I 
p09_te <-      plot.level(filter(BARBI_mean, dataset=='test'), prevalence100, I_nonev_be, model, 0, 1, 'I0+')
p10_te <-      plot.level(filter(BARBI_mean, dataset=='test'), prevalence100, I_nonev_wo, model, 0, 1, 'I0-')
p11_te <-      plot.level(filter(BARBI_mean, dataset=='test'), prevalence100, I_event_wo, model, 0, 1, 'I1-')
p12_te <- plot.level.axis(filter(BARBI_mean, dataset=='test'), prevalence100, I_event_be, model, 0, 1, 'I1+')

# lvl 2
p13_te <-      plot.level(filter(BARBI_mean, dataset=='test'), prevalence100, BA_nonev, model, min_ba_l2, max_ba_l2, 'BA0', span = 0.9)
p14_te <- plot.level.axis(filter(BARBI_mean, dataset=='test'), prevalence100, BA_event, model, min_ba_l2, max_ba_l2, 'BA1', span = 0.9)
p15_te <-      plot.level(filter(BARBI_mean, dataset=='test'), prevalence100, RB_nonev, model, -0.9, 0.5,            'RB0', span = 0.9)
p16_te <- plot.level.axis(filter(BARBI_mean, dataset=='test'), prevalence100, RB_event, model, -0.9, 0.5,            'RB1', span = 0.9)
p17_te <-      plot.level(filter(BARBI_mean, dataset=='test'), prevalence100, I_nonev,  model, -0.3, 0.6,            'I0' , span = 0.9)
p18_te <- plot.level.axis(filter(BARBI_mean, dataset=='test'), prevalence100, I_event,  model, -0.3, 0.6,            'I1' , span = 0.9)

# lvl 3
p19_te <- plot.level3(filter(BARBI_mean, dataset=='test'), prevalence100, BA, model, -0.002, 0.04, ylab = 'BA', span = 0.8)
p20_te <- plot.level3(filter(BARBI_mean, dataset=='test'), prevalence100, RB, model, -1, 0.5, ylab = 'RB', span = 1.5)
p21_te <- plot.level3(filter(BARBI_mean, dataset=='test'), prevalence100, I,  model, -0.2, 0.6, ylab = 'I' , span = 1)

# LEVELS final plot 
txt_title_test <- textGrob("Trends in BA-RB-I coefficients (test dataset)", gp = gpar(fontisize = 12))


tiff('results/LEVELS_trends_test.tiff', res = 300, width = 0.95 * 1855, height = 0.95 * 2625, units = 'px')
grid.arrange(txt_title_test,
             txt_lvl1, txt_lvl2, txt_lvl3, 
             
             rbind(ggplotGrob(p01_te), 
                   ggplotGrob(p02_te), 
                   ggplotGrob(p03_te), 
                   ggplotGrob(p04_te),
                   ggplotGrob(p05_te), 
                   ggplotGrob(p06_te), 
                   ggplotGrob(p07_te), 
                   ggplotGrob(p08_te),
                   ggplotGrob(p09_te), 
                   ggplotGrob(p10_te), 
                   ggplotGrob(p11_te), 
                   ggplotGrob(p12_te), size = 'first'), 
             
             rbind(rbind(ggplotGrob(p13_te),
                         ggplotGrob(p14_te), size = 'first'), 
                   rbind(ggplotGrob(p15_te),
                         ggplotGrob(p16_te), size = 'first'),
                   rbind(ggplotGrob(p17_te),
                         ggplotGrob(p18_te), size = 'first') ), 
             
             rbind(ggplotGrob(p19_te), 
                   ggplotGrob(p20_te), 
                   ggplotGrob(p21_te)),
             
             lvl_legend, 
             
             #layout_matrix = rbind(c(1,2,3),
             #                      c(4,5,6),
             #                      c(7,7,7)),
             #
             #heights = c(0.4, 10, 0.4)
             layout_matrix = rbind(c(1,1,1),
                                   c(2,3,4),
                                   c(5,6,7),
                                   c(8,8,8)),
             heights = c(0.4, 0.3, 10, 0.4)
             
)
dev.off()
#

