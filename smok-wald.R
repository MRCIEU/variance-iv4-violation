load("data/pheno.RData")

library("ieugwasr")
library("dplyr")
library("ggplot2")
library("TwoSampleMR")
library("broom")
library("rbgen")
source("funs.R")
options(ieugwasr_api="http://64.227.44.193:8006/")
set.seed(123)

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")

# load covariates
pc <- get_genetic_principal_components()

# merge data
dat <- merge(linker, pc, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")

# SD scale outcomes
dat$forced_vital_capacity_best_measure.20151.0.0 <- dat$forced_vital_capacity_best_measure.20151.0.0 / sd(dat$forced_vital_capacity_best_measure.20151.0.0, na.rm=T)
dat$forced_expiratory_volume_best_measure.20150.0.0 <- dat$forced_expiratory_volume_best_measure.20150.0.0 / sd(dat$forced_expiratory_volume_best_measure.20150.0.0, na.rm=T)
dat$number_of_cigarettes_previously_smoked_daily.2887.0.0 <- dat$number_of_cigarettes_previously_smoked_daily.2887.0.0 / sd(dat$number_of_cigarettes_previously_smoked_daily.2887.0.0, na.rm=T)

# select SNPs and extract from UKBB
iv <- ieugwasr::tophits("ieu-b-25")
for (i in 1:nrow(iv)){
    dosage <- tryCatch(
        {
            extract_variant_from_bgen(iv$chr[i], as.double(iv$position[i]), iv$nea[i], iv$ea[i])
        },
        error=function(cond) {
            return(NULL)
        }
    )
    if (!is.null(dosage)){
        dat <- merge(dat, dosage, "appieu")
    }
}

# test for SNP effects
exp_df <- data.frame()
out_df <- data.frame()
for (i in 1:nrow(iv)){
    snp <- paste0("chr", iv$chr[i], "_", as.double(iv$position[i]), "_", iv$nea[i], "_", iv$ea[i])

    if (!snp %in% names(dat)){
        next
    }

    # IV-exp effect
    fit <- lm(as.formula(paste0("number_of_cigarettes_previously_smoked_daily.2887.0.0 ~ ", snp, " + sex.31.0.0 + age_at_recruitment.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=dat)
    b_exp <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(estimate)
    se_exp <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(std.error)
    covar <- c("sex.31.0.0", "age_at_recruitment.21022.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
    #v_exp <- model(dat %>% dplyr::select(number_of_cigarettes_previously_smoked_daily.2887.0.0, !!snp, sex.31.0.0, age_at_recruitment.21022.0.0, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% tidyr::drop_na(.), snp, "number_of_cigarettes_previously_smoked_daily.2887.0.0", covar1 = covar, covar2 = covar)

    # interaction effect
    fit <- lm(as.formula(paste0("number_of_cigarettes_previously_smoked_daily.2887.0.0 ~ smoking_status.20116.0.0 *", snp, " + sex.31.0.0 + age_at_recruitment.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=dat)
    p_int <- fit %>% tidy %>% dplyr::filter(grepl(":",term)) %>% dplyr::pull(p.value)

    # IV-outcome
    fev_fit <- lm(as.formula(paste0("forced_expiratory_volume_best_measure.20150.0.0 ~ ", snp, " + sex.31.0.0 + age_at_recruitment.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=dat)
    fev_b <- fev_fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(estimate)
    fev_se <- fev_fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(std.error)
    fvc_fit <- lm(as.formula(paste0("forced_vital_capacity_best_measure.20151.0.0 ~ ", snp, " + sex.31.0.0 + age_at_recruitment.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=dat)
    fvc_b <- fvc_fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(estimate)
    fvc_se <- fvc_fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(std.error)

    exp_df <- rbind(exp_df, data.frame(
        Phenotype="number_of_cigarettes_previously_smoked_daily.2887.0.0",
        SNP=snp, 
        beta=b_exp, 
        se=se_exp,
        effect_allele=iv$ea[i],
        other_allele=iv$nea[i],
        units="SD"
    ))
    out_df <- rbind(out_df, data.frame(
        Phenotype=c("forced_expiratory_volume_best_measure.20150.0.0", "forced_vital_capacity_best_measure.20151.0.0"),
        SNP=c(snp, snp), 
        beta=c(fev_b, fvc_b),
        se=c(fev_se, fvc_se),
        effect_allele=c(iv$ea[i], iv$ea[i]),
        other_allele=c(iv$nea[i], iv$nea[i]),
        units=c("SD", "SD")
    ))
}

# check for X-Y effect modification
lm(forced_expiratory_volume_best_measure.20150.0.0 ~ number_of_cigarettes_previously_smoked_daily.2887.0.0 * smoking_status.20116.0.0 + sex.31.0.0 + age_at_recruitment.21022.0.0, data=dat) %>% tidy
lm(forced_vital_capacity_best_measure.20151.0.0 ~ number_of_cigarettes_previously_smoked_daily.2887.0.0 * smoking_status.20116.0.0 + sex.31.0.0 + age_at_recruitment.21022.0.0, data=dat) %>% tidy

# Wald effect of CHRNA3 on lung function
exposure_dat <- format_data(exp_df %>% dplyr::filter(SNP == "chr15_78806023_T_C"), type="exposure")
outcome_dat <- format_data(out_df %>% dplyr::filter(SNP == "chr15_78806023_T_C"), type="outcome")
mr_dat <- harmonise_data(exposure_dat, outcome_dat, action=1)
wald_res <- mr(mr_dat)

# IVW effect wo CHRNA3 on lung function
exposure_dat <- format_data(exp_df %>% dplyr::filter(SNP != "chr15_78806023_T_C"), type="exposure")
outcome_dat <- format_data(out_df %>% dplyr::filter(SNP != "chr15_78806023_T_C"), type="outcome")
mr_dat <- harmonise_data(exposure_dat, outcome_dat, action=1)
ivw_res <- mr(mr_dat)

# combine MR estimates
mr_res <- rbind(wald_res, ivw_res)
mr_res$method <- as.character(mr_res$method)
mr_res$outcome <- as.character(mr_res$outcome)
mr_res$lci <- mr_res$b - (1.96 * mr_res$se)
mr_res$uci <- mr_res$b + (1.96 * mr_res$se)
mr_res$method[mr_res$method == "Wald ratio"] <- "CHRNA3 (Wald ratio)"
mr_res$outcome[mr_res$outcome == "forced_expiratory_volume_best_measure.20150.0.0"] <- "FEV1"
mr_res$outcome[mr_res$outcome == "forced_vital_capacity_best_measure.20151.0.0"] <- "FVC"
mr_res$method <- factor(mr_res$method, levels = rev(c("CHRNA3 (Wald ratio)", "Inverse variance weighted", "Weighted median", "MR Egger", "Simple mode", "Weighted mode")))

# plot
pdf("smokingh_fev1.pdf")
ggplot(mr_res, aes(x=method, y=b, ymin=lci, ymax=uci, shape=outcome, group=outcome)) +
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
dev.off()