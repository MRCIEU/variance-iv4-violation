library("dplyr")
library("broom")
library("varGWASR")

#' Function to simulate genotypes in HWE
#' @param q Recessive/alternative allele frequency
#' @param n_obs Number of observations to return
get_simulated_genotypes <- function(q, n_obs){
  p <- 1 - q
  x <- sample(c(0, 1, 2), n_obs, prob=c(p^2, 2 * p * q, q^2), replace=T)
  return(x)
}

get_est <- function(z, x, y){
    exp_fit <- lm(x ~ z)
    b_exp <- exp_fit %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
    se_exp <- exp_fit %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("std.error") %>% as.numeric
    f_exp <- summary(exp_fit)$fstatistic[['value']]
    out_fit <- lm(y ~ z)
    b_out <- out_fit %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
    se_out <- out_fit %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("std.error") %>% as.numeric
    mr <- TwoSampleMR::mr_wald_ratio(b_exp, b_out, se_exp, se_out, NULL)
    b_mr <- mr$b
    se_mr <- mr$se
    p_mr <- mr$pval
    b_obs <- lm(y ~ x) %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::select("estimate") %>% as.numeric
    fit1 <- varGWASR::model(data.frame(z,x), "z", "x")
    names(fit1) <- paste0(names(fit1), ".zx")
    fit2 <- varGWASR::model(data.frame(z,y), "z", "y")
    names(fit2) <- paste0(names(fit2), ".zy")
    fit <- cbind(fit1, fit2)
    fit$b_mr <- b_mr
    fit$se_mr <- se_mr
    fit$b_obs <- b_obs
    fit$p_mr <- p_mr
    fit$lin_var <- var(x[z==2]) - var(x[z==0])
    fit$f_exp <- f_exp
    return(fit)
}

get_plot <- function(results){
  results$iv1 <- results$f_exp > 10
  p <- ggplot(results, aes(x=-log10(phi_p), y=b_mr, color=iv1)) +
    geom_point(alpha=0.7) +
    facet_grid(~ pb, scale="free_x") +
    scale_colour_grey(start = 0.8,end = 0.2) +
    labs(x="Variance test -log10(P)", y="MR point estimate", color="IV1 (F>10)") +
    theme_classic() +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    geom_hline(yintercept = mean(b_mr), color="grey", linetype="dashed") +
    theme(legend.position="bottom", legend.box.background = element_rect(colour = "black"))
    return(p)
}