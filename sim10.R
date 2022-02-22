library("dplyr")
library("ggplot2")
library("TwoSampleMR")
library("viridis")
source("funs.R")
set.seed(123)

n_sim <- 30
n_obs <- 10000
n_snps <- 10

# variance explained
r2_u <- 0.05
r2_z <- 0.05
r2_x <- 0.05

results <- data.frame()
for (r2_zu in seq(0, 0.1, 0.02)){
    for (r2_xu in seq(0, 0.1, 0.02)){
        for (i in 1:n_sim){
            # betas
            u_b <- sqrt(r2_u)
            z_b <- sqrt(r2_z/n_snps)
            zu_b <- sqrt(r2_zu)
            xu_b <- sqrt(r2_xu)
            x_b <- sqrt(r2_x)

            # simulate variables
            z <- sapply(runif(n_snps, min=0.05, max=0.5), function(q) get_simulated_genotypes(q, n_obs)); z <- scale(z)
            u <- rbinom(n_obs, 1, 0.5); u <- scale(u)
            x <- rowSums(z*z_b) + u*u_b + z[,1]*u*zu_b + rnorm(n_obs, sd=sqrt(1-(r2_z+r2_u+r2_zu))); x <- scale(x)
            y <- x*x_b + u*u_b + x*u*xu_b + rnorm(n_obs, sd=sqrt(1-(r2_x+r2_u+r2_xu))); y <- scale(y) # variance explained by U is greater than set due to effect on X

            # estimate MR effect
            b_exp <- sapply(1:n_snps, function(n) lm(x ~ z[,n]) %>% tidy %>% dplyr::filter(term == "z[, n]") %>% dplyr::pull(estimate))
            se_exp <- sapply(1:n_snps, function(n) lm(x ~ z[,n]) %>% tidy %>% dplyr::filter(term == "z[, n]") %>% dplyr::pull(std.error))
            b_out <- sapply(1:n_snps, function(n) lm(y ~ z[,n]) %>% tidy %>% dplyr::filter(term == "z[, n]") %>% dplyr::pull(estimate))
            se_out <- sapply(1:n_snps, function(n) lm(y ~ z[,n]) %>% tidy %>% dplyr::filter(term == "z[, n]") %>% dplyr::pull(std.error))
            exp_df <- data.frame(SNP=paste0("rs", seq(1, n_snps)), beta=b_exp, se=se_exp, effect_allele=rep("A", n_snps))
            exp_df <- format_data(exp_df, type="exposure")
            out_df <- data.frame(SNP=paste0("rs", seq(1, n_snps)), beta=b_out, se=se_out, effect_allele=rep("A", n_snps))
            out_df <- format_data(out_df, type="outcome")
            dat <- harmonise_data(exp_df, out_df, action=1)
            het_p <- mr_heterogeneity(dat) %>% dplyr::filter(method == "Inverse variance weighted") %>% dplyr::pull(Q_pval)
            ivw <- TwoSampleMR::mr_ivw(b_exp, b_out, se_exp, se_out)
            wald <- TwoSampleMR::mr_wald_ratio(b_exp[1], b_out[1], se_exp[1], se_out[1])
            
            # variance effect
            result <- model(data.frame(z1=z[,1],x), "z1", "x")
            
            # store results
            result$het_p <- het_p
            result$b_ivw <- ivw$b
            result$b_wald <- wald$b
            result$r2_zu <- r2_zu
            result$r2_xu <- r2_xu
            result$x_b <- x_b
            results <- rbind(results, result)
        }
    }
}

# TODO
estimates <- results %>%
    dplyr::group_by(r2_zu, r2_xu) %>% 
    dplyr::summarize(cbind(t.test(b_mr) %>% tidy, phi_p.x=mean(phi_p.x)))

pdf("bias.pdf")
ggplot(estimates, aes(x=r2_zu, y=estimate, ymin=conf.low, ymax=conf.high, color=-log10(phi_p.x))) +
    geom_point() + theme_classic() + geom_hline(yintercept=results$b[1], linetype="dashed", color="grey") +
    labs(y="Wald estimate (95% CI)",x="Z-X interaction explained variance") + facet_grid(~r2_xu) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
    geom_errorbar() +
    labs(color="-log10(P_variance) Z-X") +
    scale_color_viridis(direction = 1) +
    theme(
            strip.background = element_blank(),
            strip.text.y = element_text(angle = 0),
            legend.position = "bottom",
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            panel.spacing = unit(1, "lines"), 
            plot.title = element_text(hjust = 0.5, size=11)
        ) +
        ggtitle("X-Y interaction explained variance")
dev.off()