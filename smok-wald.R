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
    f_exp <- summary(fit)$fstatistic[1] %>% as.numeric
    b_exp <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(estimate)
    se_exp <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(std.error)
    covar <- c("sex.31.0.0", "age_at_recruitment.21022.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
    v_exp <- model(dat %>% dplyr::select(number_of_cigarettes_previously_smoked_daily.2887.0.0, !!snp, sex.31.0.0, age_at_recruitment.21022.0.0, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% tidyr::drop_na(.), snp, "number_of_cigarettes_previously_smoked_daily.2887.0.0", covar1 = covar, covar2 = covar)

    # interaction effect
    fit <- lm(as.formula(paste0("number_of_cigarettes_previously_smoked_daily.2887.0.0 ~ smoking_status.20116.0.0 *", snp, " + sex.31.0.0 + age_at_recruitment.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=dat)
    p_int <- fit %>% tidy %>% dplyr::filter(grepl(":",term)) %>% dplyr::pull(p.value)

    # IV-outcome
    fit <- lm(as.formula(paste0("forced_expiratory_volume_best_measure.20150.0.0 ~ ", snp, " + sex.31.0.0 + age_at_recruitment.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=dat)
    b_out <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(estimate)
    se_out <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(std.error)

    exp_df <- rbind(exp_df, data.frame(
        Phenotype="number_of_cigarettes_previously_smoked_daily.2887.0.0",
        SNP=snp, 
        beta=b_exp, 
        se=se_exp,
        effect_allele=iv$ea[i],
        other_allele=iv$nea[i],
        units="SD",
        f_exp, p_int, phi_p=v_exp$phi_p
    ))
    out_df <- rbind(out_df, data.frame(
        Phenotype="forced_expiratory_volume_best_measure.20150.0.0",
        SNP=snp, 
        beta=b_out, 
        se=se_out,
        effect_allele=iv$ea[i],
        other_allele=iv$nea[i],
        units="SD"
    ))
}

# check for X-Y effect modification
lm(forced_expiratory_volume_best_measure.20150.0.0 ~ number_of_cigarettes_previously_smoked_daily.2887.0.0 * smoking_status.20116.0.0 + sex.31.0.0 + age_at_recruitment.21022.0.0, data=dat,family="binomial") %>% summary

# Wald effect of CHRNA3 on lung function
exposure_dat <- format_data(exp_df %>% dplyr::filter(SNP == "chr15_78806023_t_c"), type="exposure")
outcome_dat <- format_data(out_df %>% dplyr::filter(SNP == "chr15_78806023_t_c"), type="outcome")
mr_dat <- harmonise_data(exposure_dat, outcome_dat, action=1)
wald_res <- mr(mr_dat)

# IVW effect wo CHRNA3 on lung function
exposure_dat <- format_data(exp_df %>% dplyr::filter(SNP != "chr15_78806023_t_c"), type="exposure")
outcome_dat <- format_data(out_df %>% dplyr::filter(SNP != "chr15_78806023_t_c"), type="outcome")
mr_dat <- harmonise_data(exposure_dat, outcome_dat, action=1)
ivw_res <- mr(mr_dat)

# plot
res_loo$lci <- res_loo$b - (res_loo$se * 1.96)
res_loo$uci <- res_loo$b + (res_loo$se * 1.96)
res_loo$delta <- res_loo$b - dplyr::filter(method == "Inverse variance weighted") %>% dplyr::pull(b)
res_loo$label <- NA
res_loo$label[res_loo$SNP=="chr15_78806023_t_c"] <- "CHRNA3"
pdf("smokingh_fev1.pdf")
ggplot(res_loo, aes(x=-log10(phi_p), y=b, ymin=lci, ymax=uci, label=label)) +
    geom_point() +
    labs(x="Instrument-exposure variance test -log10(P)", y="IVW estimate (SD, 95% CI)") +
    coord_flip() +
    geom_text(nudge_y = 0.025) +
    geom_hline(yintercept=res %>% dplyr::filter(method == "Inverse variance weighted") %>% dplyr::pull(b)) +
    theme_classic() +
    geom_segment(mapping=aes(x=-log10(phi_p), y=res %>% dplyr::filter(method == "Inverse variance weighted") %>% dplyr::pull(b), xend=-log10(phi_p), yend=b)) +
    geom_rect(
        inherit.aes = F,
        alpha = 0.01,
        fill="#bdbdbd",
        aes(
            xmin=-Inf, 
            xmax=Inf, 
            ymin=res %>% dplyr::mutate(lci=b-(1.96*se)) %>% dplyr::filter(method == "Inverse variance weighted") %>% dplyr::pull(lci),
            ymax=res %>% dplyr::mutate(uci=b+(1.96*se)) %>% dplyr::filter(method == "Inverse variance weighted") %>% dplyr::pull(uci)
        )
    )
dev.off()