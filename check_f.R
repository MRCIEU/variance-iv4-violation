library("broom")
source("funs.R")
set.seed(123)

n_obs <- 1000
n_sim <- 1000
r2 <- 0.05
maf <- 0.25

# single SNP
# check estimated r2 and F match true values
results <- data.frame()
for (i in 1:n_sim){
    x <- get_simulated_genotypes(maf, n_obs); x <- scale(x)
    y <- x*sqrt(r2) + rnorm(n_obs, sd=sqrt(1-r2))
    fit <- lm(y ~ x)
    er2 <- summary(fit)$r.squared
    tfit <- tidy(fit)
    b <- tfit$estimate[2]
    se <- tfit$std.error[2]
    f <- summary(fit)$fstatistic[1] %>% as.numeric
    results <- rbind(results, data.frame(r2, f, ef=get_f(r2, n_obs, 1), er2, eer2=get_r2(b, se, 1-maf, n_obs), n_obs))
}

# 10 SNPs
# check estimated r2 and F match true values
r2 <- runif(10, max=.1)
results <- data.frame()
for (i in 1:n_sim){
    x <- sapply(1:10, function(x) get_simulated_genotypes(maf, n_obs)); x <- scale(x)
    y <- rowSums(t(t(x)*sqrt(r2))) + rnorm(n_obs, sd=sqrt(1-sum(r2)))
    fit <- lm(y ~ x)
    rss <- sapply(1:10, function(k) summary(lm(y ~ x[,k]))$r.squared)
    fss <- sapply(1:10, function(k) summary(lm(y ~ x[,k]))$fstatistic[1] %>% as.numeric)
    er2 <- summary(fit)$r.squared
    tfit <- tidy(fit)
    b <- tfit$estimate[-1]
    se <- tfit$std.error[-1]
    rs <- sapply(1:10, function(j) get_r2(b[j], se[j], maf, n_obs))
    fs <- sapply(1:10, function(j) get_f(r2[j], n_obs, 1))
    fs10 <- summary(lm(y ~ x))$fstatistic[1] %>% as.numeric
    results <- rbind(results, data.frame(b, se))
    #f <- summary(fit)$fstatistic[1] %>% as.numeric
    #results <- rbind(results, data.frame(r2, f, ef=get_f(r2, n_obs, 1), er2, eer2=get_r2(b, se, 1-maf, n_obs), n_obs))
}