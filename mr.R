library("TwoSampleMR")
library("ieugwasr")
library("ggplot2")
library("viridis")
source("funs.R")
options(ieugwasr_api="http://104.248.164.99/")
set.seed(123)

get_mr <- function(exp_id, out_id, vgwas, q){
    # Get instruments
    exposure_dat <- extract_instruments(exp_id)

    # add R^2
    exposure_dat$r2.exposure <- sapply(1:nrow(exposure_dat), function(k) get_r2(exposure_dat$beta.exposure[k], exposure_dat$se.exposure[k], exposure_dat$eaf.exposure[k], exposure_dat$samplesize.exposure[k]))

    # add F-stat
    exposure_dat$F.exposure <- sapply(1:nrow(exposure_dat), function(k) get_f(exposure_dat$r2[k], exposure_dat$samplesize.exposure[k], 1))

    # Merge instruments with vQTL P
    exposure_dat <- merge(
        exposure_dat,
        vgwas %>% dplyr::select(rsid, phi_p),
        by.x="SNP",
        by.y="rsid"
    )

    # define P-threshold using quantile
    p_var <- quantile(exposure_dat$phi_p, q)

    # Drop top q SNPs with IV-exp variance effect & record IV-exp P-dist to test for enrichment
    f_keep <- exposure_dat %>% dplyr::filter(phi_p >= !!p_var) %>% dplyr::pull(F.exposure)
    f_drop <- exposure_dat %>% dplyr::filter(phi_p < !!p_var) %>% dplyr::pull(F.exposure)
    f_all <- exposure_dat$F.exposure
    exposure_dat <- exposure_dat %>% dplyr::filter(phi_p >= !!p_var)
    f_mean <- wilcox.test(exposure_dat$F.exposure, conf.int=T, exact = FALSE) %>% tidy

    # Test for enrichment of IV-exp assocition using Mann Whitney U test
    w_all_diff <- wilcox.test(x=f_all, y=f_keep, exact = FALSE) %>% tidy %>% dplyr::pull(p.value)
    if (q == 0){
        w_diff <- NA
    } else {
        w_diff <- wilcox.test(x=f_drop, y=f_keep, exact = FALSE) %>% tidy %>% dplyr::pull(p.value)
    }

    # Get effects of instruments on outcome
    outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=out_id, proxies = F)

    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure_dat, outcome_dat)
    stopifnot(all(dat$effect_allele.exposure == dat$effect_allele.outcome))
    stopifnot(all(dat$other_allele.exposure == dat$other_allele.outcome))

    # Perform MR
    res <- mr(dat, method_list=c("mr_ivw"))
    
    # store results
    res$f_mean <- f_mean$estimate
    res$f_mean_l <- f_mean$conf.low
    res$f_mean_h <- f_mean$conf.high
    res$w_all_diff <- w_all_diff
    res$w_diff <- w_diff
    res$p_var <- p_var
    res$q <- q

    # save instruments to file for paper
    if (q == 0){
        write.table(dat, file=paste0(exp_id, "_", out_id, ".txt"), sep="\t", quote=F)
    }

    return(res)
}

wrapper <- function(vgwas, exp_id, out_id, label, or=T){
    mr1 <- get_mr(exp_id, out_id, vgwas, 0.75)
    mr2 <- get_mr(exp_id, out_id, vgwas, 0.5)
    mr3 <- get_mr(exp_id, out_id, vgwas, 0.25)
    mr4 <- get_mr(exp_id, out_id, vgwas, 0.1)
    mr5 <- get_mr(exp_id, out_id, vgwas, 0.05)
    mr6 <- get_mr(exp_id, out_id, vgwas, 0.0)

    mr <- rbind(mr1, mr2, mr3, mr4, mr5, mr6)

    mr$lci <- mr$b - (mr$se * 1.96)
    mr$uci <- mr$b + (mr$se * 1.96)

    if (or){
        mr$b <- exp(mr$b)
        mr$lci <- exp(mr$lci)
        mr$uci <- exp(mr$uci)
    }
    
    mr$trait <- label
    mr$q <- as.factor(mr$q)

    return(mr)
}

ldl_vgwas <- get_variants("ldl_direct.30780.0.0")
glucose_vgwas <- get_variants("glucose.30740.0.0")
urate_vgwas <- get_variants("urate.30880.0.0")

ldl <- wrapper(ldl_vgwas, "ukb-d-30780_irnt", "ieu-a-7", "LDL-CAD")
glucose <- wrapper(glucose_vgwas, "ukb-d-30740_irnt", "ieu-a-24", "Glucose-T2DM")
urate <- wrapper(urate_vgwas, "ukb-d-30880_irnt", "ieu-a-1055", "Urate-Gout")
results <- rbind(ldl, glucose, urate)

# plot binary outcomes
pdf("binary.pdf")
results$q <- factor(results$q, levels = rev(levels(results$q)))
ggplot(results, aes(x=q, y=b, ymin=lci, ymax=uci, color=-log10(w_all_diff))) +
    geom_point() +
    geom_errorbar(width=0.3) +
    coord_flip() +
    geom_hline(yintercept=1, linetype="dashed", color="grey") +
    labs(y="OR (95% CI)", x="Proportion of top instrument-variance effects removed", color="Mann-Whitney test -log10(P) for weak instrument bias") +
    theme_classic() +
    facet_grid(~trait, scales="free_x") +
    theme(
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.spacing = unit(1, "lines")
    )
dev.off()

save.image(file = "data/mr.RData")