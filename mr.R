library("TwoSampleMR")
library("ieugwasr")
library("ggplot2")
library("viridis")
source("funs.R")
options(ieugwasr_api="http://64.227.44.193:8006/")
set.seed(123)

# Taken from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5446088/#E1
Isq = function(y,s){
    k = length(y)
    w = 1/s^2; sum.w = sum(w)
    mu.hat = sum(y*w)/sum.w
    Q = sum(w*(y-mu.hat)^2)
    Isq = (Q - (k-1))/Q
    Isq = max(0,Isq)
    return(Isq)
}

get_mr <- function(exp_id, out_id, vgwas, q){
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
    exposure_dat <- exposure_dat %>% dplyr::filter(phi_p >= !!p_var)

    # Get effects of instruments on outcome
    outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=out_id, proxies = F)

    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure_dat, outcome_dat)
    stopifnot(all(dat$effect_allele.exposure == dat$effect_allele.outcome))
    stopifnot(all(dat$other_allele.exposure == dat$other_allele.outcome))

    # Perform MR
    res <- mr(dat, method_list=c("mr_egger_regression"))
    res$Isq <- Isq(dat$beta.exposure/dat$se.outcome,dat$se.exposure/dat$se.outcome)
    res$p_var <- p_var
    res$q <- q

    return(res)
}

wrapper <- function(trait, exp_id, out_id, label, or=T){
    vgwas <- get_variants(trait)

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

ldl <- wrapper("ldl_direct.30780.0.0", "ukb-d-30780_irnt", "ieu-a-7", "LDL-CAD")
hba1c <- wrapper("glycated_haemoglobin.30750.0.0", "ukb-d-30750_irnt", "ieu-a-24", "HbA1c-T2DM")
urate <- wrapper("urate.30880.0.0", "ukb-d-30880_irnt", "ieu-a-1055", "Urate-Gout")
results <- rbind(ldl, hba1c, urate)

# plot effects
pdf("plot.pdf")
results$q <- factor(results$q, levels = rev(levels(results$q)))
ggplot(results, aes(x=q, y=b, ymin=lci, ymax=uci, color=Isq)) +
    geom_point() +
    geom_errorbar(width=0.3) +
    coord_flip() +
    geom_hline(yintercept=1, linetype="dashed", color="grey") +
    labs(y="OR (95% CI)", x="Proportion of top instrument-variance effects removed", color="I^2GX") +
    scale_color_viridis(direction = 1) +
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