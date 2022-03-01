library("TwoSampleMR")
library("ieugwasr")
library("ggplot2")
library("viridis")
library("meta")
source("funs.R")
options(ieugwasr_api="http://64.227.44.193:8006/")
set.seed(123)

get_meta <- function(exp_id, out_id, vgwas, p_threshold=0.05){
    # Get instruments
    exposure_dat <- extract_instruments(exp_id)

    # Merge instruments with vQTL P
    exposure_dat <- merge(
        exposure_dat,
        vgwas %>% dplyr::select(rsid, phi_p),
        by.x="SNP",
        by.y="rsid"
    )

    # split instruments by p_threshold
    exposure_dat_h0 <- exposure_dat %>% dplyr::filter(phi_p > p_threshold)
    exposure_dat_h1 <- exposure_dat %>% dplyr::filter(phi_p <= p_threshold)

    # Get effects of instruments on outcome
    outcome_dat_h0 <- extract_outcome_data(snps=exposure_dat_h0$SNP, outcomes=out_id, proxies = F)
    outcome_dat_h1 <- extract_outcome_data(snps=exposure_dat_h1$SNP, outcomes=out_id, proxies = F)

    # Harmonise the exposure and outcome data
    dat_h0 <- harmonise_data(exposure_dat_h0, outcome_dat_h0)
    dat_h1 <- harmonise_data(exposure_dat_h1, outcome_dat_h1)

    # Perform MR
    res_h0 <- mr(dat_h0, method_list=c("mr_ivw"))
    res_h1 <- mr(dat_h1, method_list=c("mr_ivw"))

    # random effect meta-analysis of IVW effects
    res <- metagen(c(res_h0$b, res_h1$b), c(res_h0$se, res_h1$se), studlab=c("H0", "H1"), comb.fixed = TRUE, comb.random = TRUE, prediction=TRUE, sm="SMD")

    return(res)
}

# load vGWAS summary stats
cigarettes_per_day <- read.table("data/number_of_cigarettes_previously_smoked_daily.2887.0.0_iv_variance.txt")
ldl_vgwas <- get_variants("ldl_direct.30780.0.0")
glucose_vgwas <- get_variants("glucose.30740.0.0")
urate_vgwas <- get_variants("urate.30880.0.0")

# estimate meta-analysis of IVW estimates for IV-exp h0 vs h1
fev1 <- get_meta("ieu-b-25", "ukb-b-19657", cigarettes_per_day)
fvc <- get_meta("ieu-b-25", "ukb-b-7953", cigarettes_per_day)
ldl <- get_meta("ukb-d-30780_irnt", "ieu-a-7", ldl_vgwas)
glucose <- get_meta("ukb-d-30740_irnt", "ieu-a-24", glucose_vgwas)
urate <- get_meta("ukb-d-30880_irnt", "ieu-a-1055", ldl_vgwas)

# plots
