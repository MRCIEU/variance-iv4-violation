library("TwoSampleMR")
library("ggplot2")
source("funs.R")
options(ieugwasr_api="http://64.227.44.193:8006/")
set.seed(123)

get_ivw <- function(exp_id, out_id, vgwas, q){
    # Get instruments
    exposure_dat <- extract_instruments(exp_id)

    # Merge instruments with vQTL P
    exposure_dat <- merge(
        exposure_dat,
        vgwas %>% dplyr::select(rsid, phi_p),
        by.x="SNP",
        by.y="rsid"
    )

    # define P-threshold using quantile
    p_var <- quantile(exposure_dat$phi_p, q)

    # Select instruments if they do not have strong vQTL effects
    exposure_dat <- exposure_dat %>% dplyr::filter(phi_p > !!p_var)

    # Get effects of instruments on outcome
    outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=out_id, proxies = F)

    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure_dat, outcome_dat)
    stopifnot(all(dat$effect_allele.exposure == dat$effect_allele.outcome))
    stopifnot(all(dat$other_allele.exposure == dat$other_allele.outcome))

    # Perform MR
    res <- mr(dat, method_list=c("mr_ivw"))
    res <- cbind(res, mr_heterogeneity(dat, method_list=c("mr_ivw"))[6:8])
    res$p_var <- p_var
    res$q <- q

    return(res)
}

# LDL -> CVD
ldl_vgwas <- get_variants("ldl_direct.30780.0.0")
ldl_mr <- rbind(
    get_ivw("ukb-d-30780_irnt", "ieu-a-7", ldl_vgwas, 0.25),
    get_ivw("ukb-d-30780_irnt", "ieu-a-7", ldl_vgwas, 0.10),
    get_ivw("ukb-d-30780_irnt", "ieu-a-7", ldl_vgwas, 0.05),
    get_ivw("ukb-d-30780_irnt", "ieu-a-7", ldl_vgwas, 0.01)
)