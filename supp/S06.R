draw <- function(d) {
  library(ggplot2)
  library(patchwork)
phenotypes <- c('HbA1c', 'BMI', 'CRP', 'Albumin', 'DM_duration_years', 'DFU_duration_days')
phenotype_en <- c('HbA1c', 'BMI', 'CRP', 'Serum albumin', 'Diabetes duration', 'DFU duration')
votus <- sort(unique(d$feature[grepl('^vOTU[0-9]+$', d$feature)]))
features <- c(votus, 'Shannon', 'Total_AMG_carrier')
top <- d[order(d$raw_p, d$feature, d$comparison), ][1:10, ]
blue <- '#0072B2'; orange <- '#E69F00'
base <- theme_classic(base_family = 'Arial', base_size = 8) + theme(
  text = element_text(colour = '#222222'), axis.text = element_text(size = 7.5, colour = '#222222'),
  axis.title = element_text(size = 8), plot.title = element_text(size = 9, face = 'bold'),
  plot.title.position = 'plot', plot.margin = margin(3, 3, 3, 3, 'mm'),
  axis.line = element_line(linewidth = .25), axis.ticks = element_line(linewidth = .25),
  panel.grid = element_blank())
row_y <- setNames((length(votus)+5) - c(seq_along(votus), length(votus)+2, length(votus)+4), features)
d$x <- match(d$comparison, phenotypes); d$y <- unname(row_y[d$feature])
d$score <- sign(d$estimate) * -log10(d$raw_p)
d$star <- ifelse(d$raw_p < .05, '*', '')
limit <- ceiling(max(abs(d$score)))
row_labels <- c(votus, 'Shannon index', 'Total AMG-carrier abundance')
a <- ggplot(d, aes(x, y)) +
  geom_tile(aes(fill = score), width = .94, height = .90) +
  geom_text(aes(label = star), size = 3.1, family = 'Arial', colour = 'black', vjust = .68) +
  scale_fill_gradient2(low = blue, mid = 'white', high = orange, midpoint = 0,
    limits = c(-limit, limit), breaks = seq(-limit, limit),
    name = expression(sign(beta) %*% -log[10](P)),
    guide = guide_colourbar(display = 'gradient',
      barwidth = grid::unit(55, 'mm'), barheight = grid::unit(2.8, 'mm'),
      title.position = 'top', title.hjust = .5)) +
  scale_x_continuous(breaks = 1:6,
    labels = c('HbA1c', 'BMI', 'CRP', 'Serum\nalbumin', 'Diabetes\nduration', 'DFU\nduration'),
    limits = c(.47, 6.53), expand = c(0, 0)) +
  scale_y_continuous(breaks = unname(row_y), labels = row_labels,
    limits = c(.45, length(votus)+4.55), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = sprintf("a   Associations across %d virome endpoints", length(features)),
    caption = '* Raw P < 0.05') + base +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
    legend.position = 'bottom', legend.title = element_text(size = 9),
    legend.text = element_text(size = 7), legend.margin = margin(0, 0, 0, 0),
    plot.caption = element_text(size = 7, hjust = 0), legend.box.margin = margin(-2, 0, 0, 0))

top$y <- 11 - seq_len(nrow(top))
top$label <- paste(top$feature, phenotype_en[match(top$comparison, phenotypes)], sep = '  |  ')
top$direction <- ifelse(top$estimate < 0, 'Negative', 'Positive')
b <- ggplot(top, aes(y = y)) +
  annotate('segment', x = 0, xend = 0, y = .45, yend = 10.55,
    colour = '#777777', linetype = 'dashed', linewidth = .35) +
  geom_segment(aes(x = ci95_lower, xend = ci95_upper, yend = y, colour = direction), linewidth = .6) +
  geom_point(aes(x = estimate, colour = direction), size = 2, shape = 16) +
  geom_text(aes(x = 1.03, label = sprintf('%.4f', raw_p)), family = 'Arial', size = 2.65) +
  geom_text(aes(x = 1.47, label = sprintf('%.4f', bh_q_recalculated)), family = 'Arial', size = 2.65) +
  geom_text(aes(x = 1.83, label = n_model_effective), family = 'Arial', size = 2.65) +
  annotate('text', x = c(1.03, 1.47, 1.83), y = 11,
    label = c('Raw P', 'BH q', 'n'), family = 'Arial', fontface = 'bold', size = 2.65) +
  annotate('segment', x = -.75, xend = .82, y = .42, yend = .42, linewidth = .25) +
  scale_colour_manual(values = c(Negative = blue, Positive = orange), guide = 'none') +
  scale_y_continuous(breaks = top$y, labels = top$label, limits = c(.4, 11.5), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(-.5, 0, .5), limits = c(-.78, 1.99), expand = c(0, 0)) +
  labs(x = expression('Standardized ' * beta * ' (95% CI)'), y = NULL,
    title = 'b   Ten associations with the smallest raw P values') + base +
  theme(axis.line = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_text(hjust = .29))
figure <- a / b + plot_layout(heights = c(1.52, 1))
return(figure)
}
