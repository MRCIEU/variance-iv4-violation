library("TwoSampleMR")
library("ieugwasr")
library("ggplot2")
library("viridis")
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

    # Drop top q SNPs with IV-exp variance effect
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

wrapper <- function(trait, exp_id, out_id, label){
    vgwas <- get_variants(trait)

    mr1 <- get_ivw(exp_id, out_id, vgwas, 0.75)
    mr2 <- get_ivw(exp_id, out_id, vgwas, 0.5)
    mr3 <- get_ivw(exp_id, out_id, vgwas, 0.25)
    mr4 <- get_ivw(exp_id, out_id, vgwas, 0.1)
    mr5 <- get_ivw(exp_id, out_id, vgwas, 0.05)
    mr6 <- get_ivw(exp_id, out_id, vgwas, 0.0)
    mr <- rbind(mr1, mr2, mr3, mr4, mr5, mr6)
    mr$lci <- exp(mr$b - (mr$se * 1.96))
    mr$uci <- exp(mr$b + (mr$se * 1.96))
    mr$b <- exp(mr$b)
    mr$trait <- label
    mr$q <- as.factor(mr$q)

    return(mr)
}

ldl <- wrapper("ldl_direct.30780.0.0", "ukb-d-30780_irnt", "ieu-a-7", "LDL-CVD")
hba1c <- wrapper("glycated_haemoglobin.30750.0.0", "ukb-d-30750_irnt", "ieu-a-24", "HbA1C-T2DM")
urate <- wrapper("urate.30880.0.0", "ukb-d-30880_irnt", "ieu-a-1055", "Urate-Gout")
results <- rbind(ldl, hba1c, urate)

# plot effects
pdf("plot.pdf")
ggplot(results, aes(x=q, y=b, ymin=lci, ymax=uci, color=-log10(Q_pval))) +
    geom_point() +
    geom_errorbar() +
    coord_flip() +
    labs(y="IVW estimate (95% CI)", x="Top instrument variance effects removed") +
    scale_color_viridis(direction = 1) +
    theme_classic() +
    facet_grid(trait~)
dev.off()