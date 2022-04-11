library("TwoSampleMR")
library("ggplot2")
set.seed(212)

# Get instruments
exposure_dat <- extract_instruments("ieu-b-25")
chrna3_exp <- exposure_dat %>% dplyr::filter(SNP=="rs8034191")
nochrna3_exp <- exposure_dat %>% dplyr::filter(SNP!="rs8034191")

# Get effects of instruments on outcome
chrna3_out <- extract_outcome_data(snps=chrna3_exp$SNP, outcomes=c("ukb-b-19657", "ukb-b-7953", "ukb-b-12019"))
nochrna3_out <- extract_outcome_data(snps=nochrna3_exp$SNP, outcomes=c("ukb-b-19657", "ukb-b-7953", "ukb-b-12019"))

# Harmonise the exposure and outcome data
chrna3_dat <- harmonise_data(chrna3_exp, chrna3_out)
nochrna3_dat <- harmonise_data(nochrna3_exp, nochrna3_out)

# Perform MR
chrna3_res <- mr(chrna3_dat)
nochrna3_res <- mr(nochrna3_dat)

# merge tables
ivw <- nochrna3_res %>% dplyr::filter(method == "Inverse variance weighted")
wald <- chrna3_res
dat <- rbind(ivw, wald)
dat$lci <- dat$b - (1.96 * dat$se)
dat$uci <- dat$b + (1.96 * dat$se)
dat$outcome[dat$outcome == "Peak expiratory flow (PEF) || id:ukb-b-12019"] <- "Peak expiratory flow"
dat$outcome[dat$outcome == "Forced vital capacity (FVC) || id:ukb-b-7953"] <- "Forced vital capacity"
dat$outcome[dat$outcome == "Forced expiratory volume in 1-second (FEV1) || id:ukb-b-19657"] <- "Forced expiratory volume"
dat$exposure[dat$exposure == "Cigarettes per Day || id:ieu-b-25"] <- "Cigarettes per day"

# Plot
#pdf("smokingh_fev1.pdf")
ggplot(dat, aes(x=method, y=b, ymin=lci, ymax=uci, shape=outcome, group=outcome)) +
    geom_point(size=2, position=position_dodge(width = .5)) +
    geom_errorbar(width = .3, position=position_dodge(width = .5)) +
    coord_flip() +
    labs(x="Method", y="Estimate (SD, 95% CI)", shape="Outcome") + 
    geom_hline(yintercept=0, color="grey", linetype = "dashed") +
    theme_classic() +
    theme(
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")
    )
#dev.off()