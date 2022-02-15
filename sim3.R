library("dplyr")
library("broom")
library("TwoSampleMR")
source("funs.R")
set.seed(123)

n_sim <- 50
n_obs <- 10000
r2_z <- 0.05 # explained variance of main effect of Z on X
r2_u <- 0.05 # explained variance of main effect of U on X
r2_c <- 0.05 # explained variance of main effect of C on X
r2_x <- 0.05 # explained variance of main effect of X on Y
NOSH <- TRUE # no simultaneous heterogeneity of Z-X and X-Y

# NOSH violation
results <- data.frame()
for (r2_zu in seq(0, 0.2, 0.05)){
    for (r2_xu in seq(0, 0.2, 0.05)){
        for (i in 1:n_sim){
            # simulate variables
            z_q <- runif(1, min=0.05, max=0.5)
            z <- get_simulated_genotypes(z_q, n_obs)
            c <- rnorm(n_obs)
            u1 <- rbinom(n_obs, 1, 0.5)
            u2 <- rbinom(n_obs, 1, 0.5)
            
            # calculate betas
            u_b <- sqrt(r2_u/0.5^2)
            c_b <- sqrt(r2_c/1)
            x_b <- sqrt(r2_x/1)
            z_p <- 1-z_q; v <- 2*z_p*z_q; z_b <- sqrt(r2_z/v)
            v_zu1 <- var(z*u1) + 2*cov(z, z*u1) + 2*cov(u1, z*u1) # variance of interaction term # TODO expectation
            zu1_b <- sqrt(r2_zu/v_zu1)

            # simulate exposure
            evx <- r2_z + r2_u + r2_zu + r2_c
            stopifnot(!evx > 1)
            x <- z*z_b + u1*u_b + z*u1*zu1_b + c*c_b + rnorm(n_obs, sd=sqrt(1 - evx))

            # simulate outcome
            v_xu1 <- var(x*u1) + 2*cov(x, x*u1) + 2*cov(u1, x*u1) # variance of interaction term # TODO expectation
            v_xu2 <- var(x*u2) + 2*cov(x, x*u2) + 2*cov(u2, x*u2) # variance of interaction term # TODO expectation
            xu1_b <- sqrt(r2_xu/v_xu1)
            xu2_b <- sqrt(r2_xu/v_xu2)
            evy <- r2_x + r2_u + r2_xu + r2_c
            stopifnot(!evy > 1)
            if (NOSH){
                y <- x*x_b + u2*u_b + x*u2*xu2_b + c*c_b + rnorm(n_obs, sd=sqrt(1 - evy))
            } else {
                y <- x*x_b + u1*u_b + x*u1*xu1_b + c*c_b + rnorm(n_obs, sd=sqrt(1 - evy))
            }

            # check effect sizes
            #results <- rbind(results, data.frame(
            #    var_x=var(x),
            #    var_y=var(y),
            #    zx_r2=summary(lm(x~z))$r.squared,
            #    ux_r2=summary(lm(x~u1))$r.squared,
            #    zux_r2=summary(lm(x~z*u1))$r.squared,
            #    cx_r2=summary(lm(x~c))$r.squared,
            #    xy_r2=summary(lm(y~x+c+u1+u2))$r.squared,
            #    u1y_r2=summary(lm(y~u1))$r.squared,
            #    u2y_r2=summary(lm(y~u2))$r.squared,
            #    xu1y_r2=summary(lm(y~x*u1))$r.squared,
            #    xu2y_r2=summary(lm(y~x*u2))$r.squared,
            #    cy_r2=summary(lm(y~c))$r.squared,
            #    r2_zu=r2_zu,
            #    r2_xu=r2_xu,
            #    NOSH
            #))

            # test causal effect & IV variance effects
            res <- get_est(z, x, y)

            # store result
            res$r2_zu <- r2_zu
            res$r2_xu <- r2_xu
            res$x_b <- x_b
            res$NOSH <- NOSH
            results <- rbind(results, res)
        }
    }
}

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
ggplot(mr_mean, aes(x=phi_zx, y=estimate, ymin=conf.low, ymax=conf.high, group=NOSH)) +
    geom_point() + geom_errorbar(width=.05) + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color="grey") +
    labs(y="Mean MR estimate (95% CI [Monte Carlo])",x="IV-X interaction variance explained") + facet_grid(~phi_xy)
dev.off()