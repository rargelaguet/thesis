library(MASS)
library(data.table)
library(purrr)
library(ggplot2)

# Simulate mean-field variational posterior with independent unobserved variables
Mu <- c(0,0)
Sigma <- matrix(c(0.5,0,0,0.5),2,2)
df1 <- mvrnorm(10000, mu = Mu, Sigma = Sigma ) %>%
  as.data.table %>% setnames(c("x","y")) %>%
  .[,key:="Mean-Field"] #%>%
  
# Simulate exact posterior with dependencies between the unobserved variables
Mu <- c(0,0)
Sigma <- matrix(c(0.75,1,1,0.75),2,2)
df2 <- mvrnorm(10000, mu = Mu, Sigma = Sigma ) %>%
  as.data.table %>% setnames(c("x","y")) %>%
  .[,key:="True Posterior"] #%>%

df <- rbind(df1,df2)

p <- ggplot(df, aes(x=x, y=y, color=key)) +
  # geom_point(alpha=0.8, size=0.25) +
  ggrastr::geom_point_rast(size=0.25) +
  # scale_fill_manual(values=colors) +
  coord_cartesian(ylim=c(-4,4), xlim=c(-4,4)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x="Variable 1", y="Variable 2") +
  theme_bw() + 
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(color="black", size=rel(1.2)),
    axis.title = element_text(size=rel(1.3))
  )

pdf(file="/Users/ricard/thesis/ricard_thesis/Chapter2/Figs/mean_field2.pdf", width = 8, height = 5)
print(p)
dev.off()