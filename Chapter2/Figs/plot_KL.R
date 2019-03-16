library(magrittr)
library(tidyr)
library(ggplot2)

# Simulate observations as a mixture model of two gaussians
y = c(rnorm(1e5, mean=-1.5, sd=0.5),rnorm(1e5, mean=1.5, sd=0.5))
# Simulate fitting with KL[p||q]
kl_pq = rnorm(1e5, mean=-1.5, sd=0.3)
# Simulate fitting with KL[q||p]
kl_qp = rnorm(2e5, mean=0, sd=1.5)

df <- data.frame(y=y, kl_pq=kl_pq, kl_qp=kl_qp) %>% gather()

colors <- c(
  y="grey70",
  kl_pq = "#66c2a5",
  kl_qp = "#fc8d62"
)

p <- ggplot(df, aes(x=value, fill=key)) +
  # geom_density(aes(y=..scaled..), alpha=0.6) +
  geom_density(alpha=0.4) +
  lims(x=c(-4.5,4.5)) +
  scale_fill_manual(values=colors) +
  # geom_vline(xintercept=0, linetype="dashed") +
  labs(x="", y="Density") +
  theme_bw() + 
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(color="black"),
    axis.title = element_text(size=rel(1.2))
  )

pdf(file="/Users/ricard/thesis/ricard_thesis/Chapter2/Figs/KL.pdf", width = 8, height = 5)
print(p)
dev.off()

df.var <- data.frame(y=var(y), kl_pq=var(kl_pq), kl_qp=var(kl_qp)) %>% gather()
p <- ggplot(df.var, aes(x=key, y=value, fill=key)) +
  # geom_density(aes(y=..scaled..), alpha=0.6) +
  geom_bar(stat="identity", color="black", alpha=0.8) +
  scale_fill_manual(values=colors) +
  # geom_vline(xintercept=0, linetype="dashed") +
  labs(x="", y="Variance") +
  theme_bw() + 
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(color="black"),
    axis.title = element_text(size=rel(1.2))
  )

pdf(file="/Users/ricard/thesis/ricard_thesis/Chapter2/Figs/KL_var.pdf", width = 8, height = 5)
print(p)
dev.off()