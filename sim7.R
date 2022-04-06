library("dplyr")
library("ggplot2")
library("TwoSampleMR")
library("viridis")
source("funs.R")
set.seed(123)

n_sim <- 30
n_obs <- 10000

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
            z_b <- sqrt(r2_z)
            zu_b <- sqrt(r2_zu)
            xu_b <- sqrt(r2_xu)
            x_b <- sqrt(r2_x)

            # simulate variables
            z <- get_simulated_genotypes(0.25, n_obs); z <- scale(z)
            u <- rnorm(n_obs)
            x <- z*z_b + u*u_b + z*u*zu_b + rnorm(n_obs, sd=sqrt(1-(r2_z+r2_u+r2_zu)))
            y <- x*x_b + u*u_b + x*u*xu_b + rnorm(n_obs, sd=sqrt(1-(r2_x+r2_u+r2_xu)))

            # estimate MR effect
            b_exp <- lm(x~z) %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(estimate)
            se_exp <- lm(x~z) %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(std.error)
            b_out <- lm(y~z) %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(estimate)
            se_out <- lm(y~z) %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(std.error)
            wald <- mr_wald_ratio(b_exp, b_out, se_exp, se_out, NULL)
            v_z <- model(data.frame(z, x), "z", "x")
            names(v_z) <- paste0(names(v_z), ".x")
            v_y <- model(data.frame(z, y), "z", "y")
            names(v_y) <- paste0(names(v_y), ".y")
            result <- cbind(v_z, v_y)
            result$r2_zu <- r2_zu
            result$r2_xu <- r2_xu
            result$b_mr <- wald$b
            result$se_mr <- wald$se_mr
            result$b <- sqrt(r2_x)

            results <- rbind(results, result)
        }
    }
}

estimates <- results %>%
    dplyr::group_by(r2_zu, r2_xu) %>% 
    dplyr::summarize(cbind(t.test(b_mr) %>% tidy, phi_p.x=mean(phi_p.x)))

estimates$r2_zu <- as.factor(estimates$r2_zu)
pdf("bias.pdf")
ggplot(estimates, aes(x=r2_zu, y=estimate, ymin=conf.low, ymax=conf.high, color=-log10(phi_p.x))) +
    geom_point() + theme_classic() + geom_hline(yintercept=results$b[1], linetype="dashed", color="grey") +
    labs(y="Wald estimate (95% CI)",x="Z-X interaction explained variance") + facet_grid(~r2_xu) +
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
            plot.title = element_text(hjust = 0.5, size=11),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        ) +
        ggtitle("X-Y interaction explained variance")
dev.off()