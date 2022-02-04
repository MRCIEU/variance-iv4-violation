library("dplyr")
library("broom")
library("varGWASR")
library("ggplot2")
source("funs.R")
set.seed(123)

n_sim <- 200
n_obs <- 100

# simulate MR estimate with IV-exp interaction effect
results <- data.frame()
for (i in 1:n_sim){
    data <- data.frame()
    for (pb in seq(0,1,0.2)){
        c <- rnorm(n_obs)
        z <- get_simulated_genotypes(0.25, n_obs)
        u <- rbinom(n_obs, 1, pb)
        x <- z*u + c + rnorm(n_obs)
        y <- x + c + rnorm(n_obs)
        data <- rbind(data, data.frame(
            x, z, u, x, y
        ))
  }
    exp_fit <- lm(x ~ z, data)
    b_exp <- exp_fit %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
    f_exp <- summary(exp_fit)$fstatistic[['value']]
    b_out <- lm(y ~ z, data) %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
    b_mr <- b_out/b_exp
    b_obs <- lm(y ~ x, data) %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::select("estimate") %>% as.numeric
    fit <- varGWASR::model(data.frame(z,x), "z", "x")
    fit$b_mr <- b_mr
    fit$b_obs <- b_obs
    fit$pb <- pb
    fit$lin_var <- var(x[z==2]) - var(x[z==0])
    fit$f_exp <- f_exp
    results <- rbind(results, fit)
}
results$iv1 <- results$f_exp > 10

# plot variance test P vs MR point estimate
pdf("plot.pdf")
ggplot(results, aes(x=-log10(phi_p), y=b_mr, color=iv1)) +
  geom_point(alpha=0.7) +
  scale_colour_grey(start = 0.8,end = 0.2) +
  labs(x="Variance test -log10(P)", y="MR point estimate", color="IV1 (F>10)") +
  theme_classic() +
  scale_x_continuous(breaks=scales::pretty_breaks()) +
  geom_hline(yintercept = mean(b_mr), color="grey", linetype="dashed") +
  theme(legend.position="bottom", legend.box.background = element_rect(colour = "black"))
dev.off()