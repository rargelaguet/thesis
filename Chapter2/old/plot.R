# curve(dnorm(x, mean=0, sd=10), col="blue", lwd=2, yaxt="n", xlim=c(-35,45))
# curve(dnorm(x, mean=10, sd=10), col="red", lwd=2, yaxt="n", xlim=c(-35,45), add=T)
# 
# curve(dnorm(x, mean=0, sd=1), col="blue", lwd=2, yaxt="n", xlim=c(-5,15))
# curve(dnorm(x, mean=10, sd=1), col="red", lwd=2, yaxt="n", xlim=c(-5,15), add=T)

library(magrittr)
library(tidyr)
library(ggplot2)

df <- data.frame(
  x1 = rnorm(1e6, mean=0, sd=5),
  x2 = rnorm(1e6, mean=10, sd=5)
) %>% gather()

p1 <- ggplot(df, aes(x=value, fill=key)) +
  geom_density(alpha=0.6) +
  lims(x=c(-15,25)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_vline(xintercept=10, linetype="dashed") +
  labs(x="", y="Density") +
  theme_bw() + 
  theme(
    legend.position = "none",
    axis.text = element_text(color="black"),
    axis.title = element_text(size=rel(1.2))
  )

df <- data.frame(
  x1 = rnorm(1e6, mean=0, sd=1),
  x2 = rnorm(1e6, mean=10, sd=1)
) %>% gather()

p2 <- ggplot(df, aes(x=value, fill=key)) +
  geom_density(alpha=0.6) +
  lims(x=c(-15,25)) +
  labs(x="", y="Density") +
  theme_bw() + 
  geom_vline(xintercept=0, linetype="dashed") +
  geom_vline(xintercept=10, linetype="dashed") +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black"),
    axis.title = element_text(size=rel(1.2))
  )

p <- cowplot::plot_grid(p1,p2, align = "v", nrow = 1)

pdf(file="/Users/ricard/thesis/ricard_thesis/Chapter2/Figs/rnorm.pdf", width = 10, height = 5)
print(p)
dev.off()