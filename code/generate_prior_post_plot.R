
prior = rnorm(10000000,0,50)
likelihood  = rnorm(10000000, 100, 32)
posterior = rnorm(10000000, 65, 32)

ydata <- data.frame(prior,likelihood, posterior)

ggplot(ydata) + 
  geom_density(aes(x= prior), fill = "#f1a340", color = "#f1a340", alpha = 0.5) +
  geom_density(aes(x= likelihood), fill = "#998ec3", color = "#998ec3", alpha = 0.5) +
  geom_density(aes(x= posterior), fill = "white", color = "grey", alpha = 0.5) +
  xlab("\npitch difference in Hz (women - men)") +
  theme_classic() +
  scale_x_continuous(breaks = (c(-200,-100,0,100,200)), limits = c(-300,300)) + 
  theme(legend.position = "right",
        legend.key.height = unit(2,"line"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_blank(),
        strip.text = element_text(size = 18, face = "bold"),
        axis.line = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"),
        plot.margin = unit(c(0.2,0.1,0.2,0.1),"cm"))

ggsave(filename = "prior_like_post.pdf",
       plot = last_plot(),
       device = "pdf",
       width = 160, 
       height = 100,
       units = "mm",
       #bg = "transparent",
       dpi = 300)

likelihood = c(100, 50, 14, 150, 186)
mean(likelihood); sd(likelihood)
xdata <- data.frame(likelihood)

prior_dummy <- c(
  # define a skeptical prior for the relevant coefficients
  prior(normal(0, 50), class = "Intercept")
)


x <- brm(likelihood ~ 1, xdata)
x_prior <- brm(likelihood ~ 1, prior = prior_dummy, xdata)
