library("dplyr")
library("broom")
library("ggplot2")
library("TwoSampleMR")
source("funs.R")
set.seed(123)

n_sim <- 1000
n_obs <- 10000
r2_z <- 0.05 # explained variance of main effect of Z on X
r2_x <- 0.05 # explained variance of main effect of X on Y
r2_u <- 0.05 # explained variance of main effect of U on X/Y
r2_c <- 0.05 # explained variance of main effect of C on X/Y

results <- data.frame()
for (r2_zu in 0){ #seq(0, 0.5, 0.1)
    for (r2_xu in 0){ #seq(0, 0.5, 0.1)
        for (i in 1:n_sim){
            # simulate variables
            z_q <- runif(1, min=0.05, max=0.5)
            z <- get_simulated_genotypes(z_q, n_obs)
            c <- rnorm(n_obs)
            u <- rbinom(n_obs, 1, 0.5)
            
            # calculate betas
            u_b <- sqrt(r2_u/0.5^2)
            c_b <- sqrt(r2_c/1)
            x_b <- sqrt(r2_x/1)
            z_p <- 1-z_q
            v <- 2*z_p*z_q
            z_b <- sqrt(r2_z/v)
            v_zu <- var(z*u) + 2*cov(z, z*u) + 2*cov(u, z*u) # variance of interaction term # TODO expectation
            zu_b <- sqrt(r2_zu/v_zu)

            # simulate exposure
            evx <- r2_z + r2_u + r2_zu + r2_c
            stopifnot(!evx > 1)
            x <- z*z_b + u*u_b + z*u*zu_b + c*c_b + rnorm(n_obs, sd=sqrt(1 - evx))

            # simulate outcome
            v_xu <- var(x*u) + 2*cov(x, x*u) + 2*cov(u, x*u) # variance of interaction term # TODO expectation
            xu_b <- sqrt(r2_xu/v_xu)
            y <- x*x_b + u*u_b + x*u*xu_b + c*c_b
            evy <- var(y)
            stopifnot(!evy > 1)
            y <- y + rnorm(n_obs, sd=sqrt(1 - evy))

            # check effect sizes
            results <- rbind(results, data.frame(
                var_x=var(x),
                var_y=var(y),
                zx_r2=summary(lm(x~z))$r.squared,
                ux_r2=summary(lm(x~u))$r.squared,
                zux_r2=summary(lm(x~z*u))$r.squared,
                cx_r2=summary(lm(x~c))$r.squared,
                xy_r2=summary(lm(y~x+c+u))$r.squared,
                uy_r2=summary(lm(y~u))$r.squared,
                uy_r2=summary(lm(y~u))$r.squared,
                xuy_r2=summary(lm(y~x*u))$r.squared,
                xuy_r2=summary(lm(y~x*u))$r.squared,
                cy_r2=summary(lm(y~c))$r.squared,
                r2_zu=r2_zu,
                r2_xu=r2_xu
            ))

            # test causal effect & IV variance effects
            #res <- get_est(z, x, y)

            # store result
            #res$r2_zu <- r2_zu
            #res$r2_xu <- r2_xu
            #res$x_b <- x_b
            #results <- rbind(results, res)
        }
    }
}

# check values
results %>% dplyr::select(-r2_zu, -r2_xu) %>% apply(., 2, function(x) t.test(x) %>% tidy)

# MR causal estimate mean
mr_mean <- results %>% 
    dplyr::group_by(NOSH,r2_zu,r2_xu) %>% 
    dplyr::summarize(t.test(b_mr) %>% tidy)

# Obs estimate mean
results %>% 
    dplyr::group_by(NOSH,r2_zu,r2_xu) %>% 
    dplyr::summarize(t.test(b_obs) %>% tidy)

# plot
pdf("bias.pdf")
ggplot(mr_mean, aes(x=r2_zu, y=estimate, ymin=conf.low, ymax=conf.high, group=NOSH)) +
    geom_point() + geom_errorbar(width=.05) + theme_classic() + geom_hline(yintercept=results$x_b[1], linetype="dashed", color="grey") +
    labs(y="Mean MR estimate (95% CI [Monte Carlo])",x="IV-X interaction variance explained") + facet_grid(r2_xu~NOSH)
dev.off()