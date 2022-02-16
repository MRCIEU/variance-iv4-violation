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