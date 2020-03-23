library(magrittr)
library(tidyr)
library(ggplot2)

df <- data.frame(
  x1 = rnorm(1e6, mean=0, sd=1/10), # alpha=10
  x2 = rnorm(1e6, mean=0, sd=1/1) # alpha=1
) %>% gather()

p <- ggplot(df, aes(x=value, fill=key)) + %>% 
  geom_density(alpha=0.6) +
  lims(x=c(-3,3)) +
  geom_vline(xintercept=0, linetype="dashed") +
  labs(x="", y="Density") +
  theme_bw() + 
  theme(
    legend.position = "top",
    axis.text = element_text(color="black"),
    axis.title = element_text(size=rel(1.2))
  )

pdf(file="/Users/ricard/thesis/ricard_thesis/Chapter2/Figs/ARD.pdf", width = 8, height = 5)
print(p)
dev.off()
