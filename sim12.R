library("dplyr")
library("ggplot2")
library("grid")
library("gtable")
library("TwoSampleMR")
library("viridis")
source("funs.R")
set.seed(123)

mr <- function(b_exp, b_out, se_exp, se_out, phi_p, q){
    # define P-threshold using quantile
    p_var <- quantile(phi_p, q)

    # Drop top q SNPs with IV-exp variance effect
    b_exp <- b_exp[phi_p >= p_var]
    b_out <- b_out[phi_p >= p_var]
    se_exp <- se_exp[phi_p >= p_var]
    se_out <- se_out[phi_p >= p_var]

    # IVW
    ivw <- mr_ivw(b_exp, b_out, se_exp, se_out, NULL)
    ivw$q <- q

    return(ivw)
}

n_obs <- 100000
n_sim <- 500
r2_u <- 0.2 # U-Y main effect (small Cohen d)
r2_x <- 0.2 # X-Y main effect (small Cohen d)
r2_xu <- r2_x * 0.5 # X-Y interaction effect half the size of the main effect
n_snps <- 6 # number of SNPs for use as IVs in MR
n_isnps <- 3 # number of IVs with interaction effects on X
phi <- 0.5 # interaction effect size relative to main effect
x_b <- sqrt(r2_x)
u_b <- sqrt(r2_u)
xu_b <- sqrt(r2_xu)

results <- data.frame()
for (i in 1:n_sim){
    # select SNP betas
    z_b <- runif(n_snps, min = 0.02, max = 0.06)
    
    # set size of interaction effect relative to main effect
    zu_b <- z_b * phi

    # simulate variables
    z <- sapply(1:n_snps, function(j) get_simulated_genotypes(0.25, n_obs)); z <- scale(z) # instruments
    u <- rnorm(n_obs) # modifier
    x <- sapply(1:n_obs, function(j) sum(z[j,]*z_b)) + u*u_b + sapply(1:n_obs, function(j) sum(z[j,1:n_isnps]*u[j]*zu_b[1:n_isnps])) + rnorm(n_obs, sd=sqrt(1-(sum(z_b^2) + u_b^2 + sum(zu_b[1:n_isnps]^2)))) # exposure
    y <- x*x_b + u*u_b + x*u*xu_b + rnorm(n_obs, sd=sqrt(1-(r2_x+r2_u+r2_xu))) # outcome
    varX <- var(x)
    varY <- var(y)
    r2 <- summary(lm(x~z))$r.squared

    # SNP estimates
    b_exp <- sapply(1:n_snps, function(i) lm(x~z[,i]) %>% tidy %>% dplyr::filter(term!="(Intercept)") %>% dplyr::pull(estimate))
    se_exp <- sapply(1:n_snps, function(i) lm(x~z[,i]) %>% tidy %>% dplyr::filter(term!="(Intercept)") %>% dplyr::pull(std.error))
    b_out <- sapply(1:n_snps, function(i) lm(y~z[,i]) %>% tidy %>% dplyr::filter(term!="(Intercept)") %>% dplyr::pull(estimate))
    se_out <- sapply(1:n_snps, function(i) lm(y~z[,i]) %>% tidy %>% dplyr::filter(term!="(Intercept)") %>% dplyr::pull(std.error))
    f_stat <- sapply(1:n_snps, function(i) summary(lm(x~z[,i]))$fstatistic[1] %>% as.numeric) # main effect only
    phi_p <- sapply(1:n_snps, function(i) model(data.frame(z, x), paste0("X", i), "x")$phi_p)

    # estimate IVW with subsets of instruemnts
    oracle <- mr_ivw(b_exp[4:6], b_out[4:6], se_exp[4:6], se_out[4:6], NULL)
    oracle$q <- NA
    results <- rbind(results, mr(b_exp, b_out, se_exp, se_out, phi_p, 0.75))
    results <- rbind(results, mr(b_exp, b_out, se_exp, se_out, phi_p, 0.5))
    results <- rbind(results, mr(b_exp, b_out, se_exp, se_out, phi_p, 0.25))
    results <- rbind(results, mr(b_exp, b_out, se_exp, se_out, phi_p, 0.1))
    results <- rbind(results, mr(b_exp, b_out, se_exp, se_out, phi_p, 0.05))
    results <- rbind(results, mr(b_exp, b_out, se_exp, se_out, phi_p, 0))
    results <- rbind(results, oracle)
}

# plot bias and efficiency
bias <- results %>% dplyr::group_by(q) %>% dplyr::summarize(t.test(b) %>% tidy) %>% dplyr::mutate(estimand="b")
bias$yintercept <- x_b
efficiency <- results %>% dplyr::group_by(q) %>% dplyr::summarize(t.test(se) %>% tidy) %>% dplyr::mutate(estimand="se")
efficiency$yintercept <- NA
s <- rbind(bias, efficiency)
s$q <- as.character(s$q)
s$q[is.na(s$q)] <- "Oracle"
s$q <- factor(s$q)
s$q <- factor(s$q, levels = rev(levels(s$q)))
pdf("sim12.pdf")
ggplot(s, aes(x=q, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() +
    geom_errorbar(width=0.3) +
    coord_flip() +
    geom_hline(data=s, aes(yintercept = yintercept), linetype="dashed", color="grey") +
    facet_grid(. ~ estimand, scales="free", switch = 'x', labeller = as_labeller(c(b="SD (95% CI)", se="SE (95% CI)"))) +
    labs(x="Proportion of top instrument-variance effects removed") +
    theme_classic() +
    theme(
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.placement = "outside",
        legend.position = "bottom",
        legend.background = element_blank(),
        axis.title.x=element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.spacing = unit(1, "lines")
    )
dev.off()