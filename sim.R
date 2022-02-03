library("dplyr")
library("broom")
library("varGWASR")
library("ggplot2")
set.seed(123)

n_sim <- 200
n_obs <- 1000

#' Function to simulate genotypes in HWE
#' @param q Recessive/alternative allele frequency
#' @param n_obs Number of observations to return
get_simulated_genotypes <- function(q, n_obs){
  p <- 1 - q
  x <- sample(c(0, 1, 2), n_obs, prob=c(p^2, 2 * p * q, q^2), replace=T)
  return(x)
}

# simulate MR estimate with IV-exp interaction effect
results <- data.frame()
for (pb in seq(0.2, 1, 0.2)){
  for (i in 1:n_sim){
      c <- rnorm(n_obs)
      z <- get_simulated_genotypes(0.25, n_obs)
      u <- rbinom(n_obs, 2, pb)
      x <- z*u + c + rnorm(n_obs)
      y <- x + c + rnorm(n_obs)
      y <- y / sd(y)
      b_exp <- lm(x ~ z) %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
      b_out <- lm(y ~ z) %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
      b_mr <- b_out/b_exp
      b_obs <- lm(y ~ x) %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::select("estimate") %>% as.numeric
      fit <- varGWASR::model(data.frame(z,x), "z", "x")
      fit$b_mr <- b_mr
      fit$b_obs <- b_obs
      fit$pb <- pb
      results <- rbind(results, fit)
  }
}

# plot
# -log10(P)
# fix the X axes
# only present the IV-exp when the F stat > 10 & use log scale
# present as a cube with variance effect size rather than P value
# apply systematically to opengwas
# change Y axis to SD
# surface plot
# plot coef 1 vs 2 and then use MR point estimate SD
pdf("plot.pdf")
ggplot(results, aes(x=-log10(phi_p), y=b_mr)) +
  geom_point(alpha=0.7) +
  facet_grid(~ pb) +
  labs(x="Variance test -log10(P)", y="MR point estimate") +
  theme_classic()
dev.off()

# plot effect size vs SD
k <- kmeans(results %>% dplyr::select(phi_x1, phi_x2) %>% as.matrix, 20)
results$cluster <- k$cluster
k_sd <- results %>% dplyr::group_by(cluster) %>% dplyr::summarize(sd=sd(b_mr))
k_sd <- cbind(k_sd, k$centers)
ggplot(k_sd, aes(x=phi_x1, y=phi_x2, size=sd)) +
  geom_point()

# simulate MR estimate with IV-exp interaction effect
results <- data.frame()
for (pb in seq(0, 1, 0.2)){
  for (i in 1:n_sim){
      c <- rnorm(n_obs)
      z <- get_simulated_genotypes(0.25, n_obs)
      u <- rbinom(n_obs, 2, pb)
      x <- z + z*u + c + rnorm(n_obs)
      y <- x + c + rnorm(n_obs)
      y <- y / sd(y)
      b_exp <- lm(x ~ z) %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
      b_out <- lm(y ~ z) %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
      b_mr <- b_out/b_exp
      b_obs <- lm(y ~ x) %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::select("estimate") %>% as.numeric
      fit <- data.frame(
        b_var=var(x[z==2]) - var(x[z==0]),
        b_mr, b_obs, pb
      )
      results <- rbind(results, fit)
  }
}

# plot effect size vs SD
ggplot(results, aes(x=b_var, y=b_mr)) +
  geom_point()