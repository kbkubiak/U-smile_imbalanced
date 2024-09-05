# examples of U-smile plot's shapes 

subclass_order <- c('nonev_dw',
                    'nonev_up',
                    'event_dw',
                    'event_up')
usmile_colors <- c('nonev_dw' = '#0F3C78',
                   'nonev_up' = '#0F3C78',
                   'event_dw' = '#D51424',
                   'event_up' = '#D51424')
usmile_fills  <- c('nonev_dw' = '#0F3C78',
                   'nonev_up' = '#BED2FA',
                   'event_dw' = '#FBCDB9',
                   'event_up' = '#D51424')
usmile_labels <- c('non-events with better prediction',
                   'non-events with worse prediction',
                   'events with worse prediction',
                   'events with better prediction')

#
panels <- data.frame(subclass = rep(c('nonev_dw',
                                      'nonev_up',
                                      'event_dw',
                                      'event_up'),9),
                     panel = factor(rep(letters[1:9], each = 4)),
                     value = c(0.103, 0.045, 0.043, 0.095, # a
                               0.102, 0.050, 0.047, 0.049, # b
                               0.046, 0.048, 0.044, 0.101, # c
                               0.045, 0.103, 0.097, 0.047, # d
                               0.053, 0.101, 0.096, 0.097, # e
                               0.098, 0.097, 0.100, 0.052, # f
                               0.051, 0.048, 0.051, 0.047, # g
                               0.097, 0.047, 0.092, 0.055, # h
                               0.045, 0.092, 0.058, 0.099  # i
                     ) )
#

tiff('results/U-smile_shapes.tiff', res = 300, 
     width = 1000, height = 1700, units = 'px')
ggplot(panels, aes(x = subclass, y = value, group = panel)) +
  geom_line() + 
  geom_point(aes(col = subclass, fill = subclass), shape = 21, size = 3.5, stroke = 1) +  # moze size = 2.5
  scale_x_discrete(limits = subclass_order) +
  scale_color_manual(values = usmile_colors,
                     breaks = subclass_order,
                     labels = usmile_labels) +
  scale_fill_manual(values = usmile_fills,
                    breaks = subclass_order,
                    labels = usmile_labels) + 
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'bottom', 
        # legend.text = element_text(size=6),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 4)  #, 
         # size  = guide_legend(nrow = 1)
  ) + 
  facet_wrap(~ panel, ncol = 3) +
  ylim(c(0, 0.13))

dev.off()
