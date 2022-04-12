library("dplyr")
library("ggplot2")
library("grid")
library("gtable")
library("TwoSampleMR")
library("viridis")
source("funs.R")
set.seed(123)

n_obs <- 100000
n_sim <- 30
r2_u <- 0.2 # U-Y main effect (small Cohen d)
r2_x <- 0.2 # X-Y main effect (small Cohen d)
r2_xu <- r2_x * 0.5 # X-Y interaction effect half the size of the main effect
n_snps <- 5 # number of SNPs for use as IVs in MR
n_isnps <- 1 # number of IVs with interaction effects on X
r2_z <- 0.05 # combined variance explained by the IV
r2_zu <- 0.005 # interaction effect size relative to main effect
p_thresh <- 0.05

# betas
z_b <- sqrt(r2_z / n_snps)
zu_b <- sqrt(r2_zu)
x_b <- sqrt(r2_x)
u_b <- sqrt(r2_u)
xu_b <- sqrt(r2_xu)

results <- data.frame()
for (i in 1:n_sim){
    # simulate variables
    q <- runif(n_snps, min = 0.05, max = 0.5) # alternative allele frequency
    q_zu <- sample(c(rep(1, n_isnps), rep(0, n_snps-n_isnps))) # IVs that have interactions
    z <- sapply(q, function(q) get_simulated_genotypes(q, n_obs)); z <- scale(z) # instrument
    u <- rnorm(n_obs) # modifier
    x <- sapply(1:n_obs, function(i) sum(z[i,]*z_b)) + u*u_b + sapply(1:n_obs, function(i) sum(z[i,]*q_zu*u[i]*zu_b)) + rnorm(n_obs, sd=sqrt(1-(r2_z+r2_u+r2_zu))) # exposure
    y <- x*x_b + u*u_b + x*u*xu_b + rnorm(n_obs, sd=sqrt(1-(r2_x+r2_u+r2_xu)))

    # MR effect
    iv_exp <- lm(x~z)
    iv_out <- lm(y~z)
    b_exp <- iv_exp %>% tidy %>% dplyr::filter(term!="(Intercept)") %>% dplyr::pull(estimate)
    se_exp <- iv_exp %>% tidy %>% dplyr::filter(term!="(Intercept)") %>% dplyr::pull(std.error)
    f_stat <- summary(lm(x~z))$fstatistic[1] %>% as.numeric
    b_out <- iv_out %>% tidy %>% dplyr::filter(term!="(Intercept)") %>% dplyr::pull(estimate)
    se_out <- iv_out %>% tidy %>% dplyr::filter(term!="(Intercept)") %>% dplyr::pull(std.error)

    # IV-exp variance effect
    phi_p <- sapply(1:n_snps, function(i) model(data.frame(z, x), paste0("X", i), "x")$phi_p)

    # IVW with all SNPs
    mr_all <- mr_ivw(b_exp, b_out, se_exp, se_out, NULL)$b

    # IVW with homogeneous subset
    b_exp <- b_exp[phi_p > p_thresh]
    b_out <- b_out[phi_p > p_thresh]
    se_exp <- se_exp[phi_p > p_thresh]
    se_out <- se_out[phi_p > p_thresh]

    if (length(b_exp) > 1){
        mr_homo <- mr_ivw(b_exp, b_out, se_exp, se_out, NULL)$b
    } else if (length(b_exp) == 0){
        mr_homo <- NA
    } else {
        mr_homo <- mr_wald_ratio(b_exp, b_out, se_exp, se_out, NULL)$b
    }

    # store result
    results <- rbind(data.frame(mr_all, mr_homo), results)
}

est <- rbind(
    t.test(results$mr_all) %>% tidy %>% dplyr::mutate(iv="All"),
    t.test(results$mr_homo) %>% tidy %>% dplyr::mutate(iv="Brown-Forsythe test P > 0.05")
)
est <- as.data.frame(est)
est$iv <- as.factor(est$iv)
levels(est$iv) <- rev(levels(est$iv))

# plot
pdf("sim10.pdf", height=3)
ggplot(est, aes(x=iv, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point(size=2, position=position_dodge(width = .5)) +
    geom_errorbar(width = .3, position=position_dodge(width = .5)) +
    coord_flip() +
    labs(x="Instruments", y="Estimate (SD, 95% CI)", shape="Outcome") + 
    geom_hline(yintercept=x_b, color="grey", linetype = "dashed") +
    theme_classic() +
    theme(
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")
    )
dev.off()