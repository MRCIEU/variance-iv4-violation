library("dplyr")
library("ggplot2")
library("TwoSampleMR")
library("viridis")
source("funs.R")
set.seed(123)

n_sim <- 50
n_obs <- 100000
r2_z <- 0.1
r2_x <- 0.1

results <- data.frame()
for (h in seq(-1, 1, .5)){
    for (s in seq(0, 1, 0.2)){
        if (h == 0 & s == 1){
            next
        }
        # betas
        b0 <- 0.5
        b1 <- b0-h*s
        b2 <- b0-s
        for (i in 1:n_sim){
            # genotype
            z <- get_simulated_genotypes(0.25, n_obs)
            u0 <- z==0
            u1 <- z==1
            u2 <- z==2
            z_ev <- var(z*u0*b0 + z*u1*b1 + z*u2*b2)
            x_tv <- z_ev / r2_z

            # exposure           
            x <- z*u0*b0 + z*u1*b1 + z*u2*b2 + rnorm(n_obs, sd=sqrt(x_tv-z_ev))

            # check IV1 holds
            f <- summary(lm(x~z))$fstatistic[1] %>% as.numeric
            if (f < 10){
                warning(paste(h, s, f))
                next
            }

            # outcome
            x_ev <- var(x)
            y_tv <- x_ev / r2_x
            y <- x + rnorm(n_obs, sd=sqrt(y_tv-x_ev))

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
            result$b_mr <- wald$b
            result$b0 <- b0
            result$b1 <- b1
            result$b2 <- b2
            result$se_mr <- wald$se_mr
            result$b <- sqrt(r2_x)
            result$h <- h
            result$s <- s
            result$f <- f

            results <- rbind(results, result)
        }
    }
}

estimates <- results %>%
    dplyr::group_by(h, s) %>% 
    dplyr::summarize(cbind(t.test(b_mr) %>% tidy, phi_p.x=mean(phi_p.x)))

pdf("monotonicity-bias.pdf")
ggplot(estimates, aes(x=b1, y=estimate, ymin=conf.low, ymax=conf.high, color=-log10(phi_p.x))) +
    geom_point() + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color="grey") +
    labs(y="Wald estimate (95% CI)",x="Heterozygote effect") + facet_grid(~b2) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=5)) +
    geom_errorbar() +
    labs(color="-log10(P_variance) Z-X") +
    scale_color_viridis(direction = 1) +
    theme(
            strip.background = element_blank(),
            axis.text.x = element_text(angle = 90),
            strip.text.y = element_text(angle = 0),
            legend.position = "bottom",
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            panel.spacing = unit(1, "lines"), 
            plot.title = element_text(hjust = 0.5, size=11)
        ) +
        ggtitle("Recessive allele effect")
dev.off()