library("TwoSampleMR")
library("ggplot2")
source("funs.R")
set.seed(123)

get_ivw <- function(exp_id, out_id, vgwas, p_var=5e-8){
    # Get instruments
    exposure_dat <- extract_instruments(exp_id)

    # Get effects of instruments on outcome
    outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=out_id)

    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure_dat, outcome_dat)
    stopifnot(all(dat$effect_allele.exposure == dat$effect_allele.outcome))
    stopifnot(all(dat$other_allele.exposure == dat$other_allele.outcome))

    # Perform MR
    res <- mr(dat)
}

# LDL -> CVD
ldl_vgwas <- 
ldl_mr <- get_ivw("ukb-d-30780_irnt", "ieu-a-7", )