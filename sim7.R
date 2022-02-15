library("dplyr")
library("TwoSampleMR")
source("funs.R")
set.seed(123)

#' Test for effect of SNP on outcome variance using the LAD-BF model
#'
#' @param data Dataframe of observations
#' @param x Name of SNP dosage
#' @param y Name of outcome
#' @param covar1 Optional vector of covariates to include in the first-stage model
#' @param covar2 Optional vector of covariates to include in the second-stage model
#' @return Dataframe containing variance effect for SNP=1 (phi_x1) and SNP=2 (phi_x2) with SE and p and F-stat
#' @export
model <- function(data, x, y, covar1=NULL, covar2=NULL){
    if (any(is.na(data))) stop("Dataframe contains NA values")
    if (!x %in% names(data)) stop(paste0(x, " was not in dataframe"))
    if (!y %in% names(data)) stop(paste0(y, " was not in dataframe"))
    if (!all(covar1 %in% names(data))) stop(paste0("Dataframe is missing some of: ", paste0(covar1, collapse=" ")))
    if (!all(covar2 %in% names(data))) stop(paste0("Dataframe is missing some of: ", paste0(covar2, collapse=" ")))
    if (!is.numeric(data[[x]])) stop("Genotype must be numeric")

    # prepare first-stage fit matrix
    if (!is.null(covar1)){
        X <- data %>% dplyr::select(!!x, !!covar1) %>% as.matrix
    } else {
        X <- data %>% dplyr::select(!!x) %>% as.matrix
    }
    # first-stage fit
    fit <- cqrReg::qrfit(X=X, y=data[[y]], tau=.5, method="mm")
    b <- rbind(fit$b, fit$beta)
    # predicted
    X <- cbind(rep(1, nrow(X)), X)
    fitted <- X %*% b
    # residual
    d <- data[[y]] - fitted
    # abs residual
    d <- abs(as.vector(d))
    # second-stage model
    data[[x]] <- as.factor(data[[x]])
    if (!is.null(covar2)){
        X <- data %>% dplyr::select(!!x, !!covar2)
        fit2 <- lm(d ~ ., data=X)
        X <- data %>% dplyr::select(!!covar2)
        fit_null <- lm(d ~ ., data=X)
        p <- anova(fit_null, fit2) %>% broom::tidy(.) %>% dplyr::pull(p.value) %>% dplyr::nth(2)
        f <- anova(fit_null, fit2) %>% broom::tidy(.) %>% dplyr::pull(statistic) %>% dplyr::nth(2)
    } else {
        X <- data %>% dplyr::select(!!x)
        fit2 <- lm(d ~ ., data=X)
        fit_null <- lm(d ~ 1, data=data)
        p <- anova(fit_null, fit2) %>% broom::tidy(.) %>% dplyr::pull(p.value) %>% dplyr::nth(2)
        f <- anova(fit_null, fit2) %>% broom::tidy(.) %>% dplyr::pull(statistic) %>% dplyr::nth(2)
    }

    # deltamethod
    v1 <- car::deltaMethod(fit2, "(2*b0*b1+b1^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"), vcov=sandwich::vcovHC(fit2, type = "HC0"))
    v2 <- car::deltaMethod(fit2, "(2*b0*b2+b2^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"), vcov=sandwich::vcovHC(fit2, type = "HC0"))

    res <- data.frame(
        phi_x1=v1$Estimate,
        se_x1=v1$SE,
        phi_x2=v2$Estimate,
        se_x2=v2$SE,
        phi_f=f,
        phi_p=p
    )
    
    return(res)
}

n_sim <- 50
n_obs <- 10000

# variance explained
r2_u <- 0.05
r2_z <- 0.05
r2_x <- 0.05

results <- data.frame()
for (r2_zu in seq(0, 0.25, 0.05)){
    for (r2_xu in seq(0, 0.25, 0.05)){
        for (i in 1:n_sim){
            # betas
            u_b <- sqrt(r2_u)
            z_b <- sqrt(r2_z)
            zu_b <- sqrt(r2_zu)
            xu_b <- sqrt(r2_xu)
            x_b <- sqrt(r2_x)

            # simulate variables
            z <- get_simulated_genotypes(0.25, n_obs); z <- scale(z)
            u <- rbinom(n_obs, 1, 0.5); u <- scale(u)
            x <- z*z_b + u*u_b + z*u*zu_b + rnorm(n_obs, sd=sqrt(1-(r2_z+r2_u+r2_zu))); x <- scale(x)
            y <- x*x_b + u*u_b + x*u*xu_b + rnorm(n_obs, sd=sqrt(1-(r2_x+r2_u+r2_xu))); y <- scale(y) # variance explained by U is greater than set due to effect on X

            # estimate MR effect
            b_exp <- lm(x~z) %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(estimate)
            se_exp <- lm(x~z) %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(std.error)
            b_out <- lm(y~z) %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(estimate)
            se_out <- lm(y~z) %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(std.error)
            wald <- mr_wald_ratio(b_exp, b_out, se_exp, se_out, NULL)
            v_z <- model(data.frame(z, x), "z", "x")
            names(v_z) <- paste0(names(v_z), ".z")
            v_x <- model(data.frame(z, x), "z", "x")
            names(v_x) <- paste0(names(v_x), ".x")
            result <- cbind(v_z, v_x)
            result$r2_zu <- r2_zu
            result$r2_xu <- r2_xu
            result$b_mr <- wald$b
            result$se_mr <- wald$se_mr
            result$b <- sqrt(r2_x)

            results <- rbind(results, result)
        }
    }
}