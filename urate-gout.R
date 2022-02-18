load("/user/home/ml18692/projects/varGWAS-ukbb-biomarkers/data/pheno.RData")

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
dat$urate.30880.0.0 <- dat$urate.30880.0.0 / sd(dat$urate.30880.0.0, na.rm=T)

# select SNPs and extract from UKBB
iv <- ieugwasr::tophits("ieu-a-1055")
for (i in 1:nrow(iv)){
    dosage <- extract_variant_from_bgen(iv$chr[i], as.double(iv$position[i]), iv$nea[i], iv$ea[i])
    dat <- merge(dat, dosage, "appieu")
}

# test for SNP-by-sex effects
exp_df <- data.frame()
out_df <- data.frame()
for (i in 1:nrow(iv)){
    snp <- paste0("chr", iv$chr[i], "_", as.double(iv$position[i]), "_", iv$nea[i], "_", iv$ea[i])

    # IV-exp effect
    fit <- lm(as.formula(paste0("urate.30880.0.0 ~ sex.31.0.0 + ", snp, " + age_at_recruitment.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=dat)
    f_exp <- summary(fit)$fstatistic[1] %>% as.numeric
    b_exp <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(estimate)
    se_exp <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(std.error)

    # interaction effect
    fit <- lm(as.formula(paste0("urate.30880.0.0 ~ sex.31.0.0 * ", snp, " + age_at_recruitment.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=dat)
    p_int <- fit %>% tidy %>% dplyr::filter(grepl(":",term)) %>% dplyr::pull(p.value)

    # IV-outcome
    fit <- glm(as.formula(paste0("gout ~ sex.31.0.0 + ", snp, " + age_at_recruitment.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=dat, family="binomial")
    b_out <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(estimate)
    se_out <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(std.error)

    exp_df <- rbind(exp_df, data.frame(
        Phenotype="Urate",
        SNP=snp, 
        beta=b_exp, 
        se=se_exp,
        effect_allele=iv$ea[i],
        other_allele=iv$nea[i],
        units="SD",
        f_exp, p_int
    ))
    out_df <- rbind(out_df, data.frame(
        Phenotype="Gout",
        SNP=snp, 
        beta=b_out, 
        se=se_out,
        effect_allele=iv$ea[i],
        other_allele=iv$nea[i],
        units="Log odds"
    ))
}

# IV-exp
exposure_dat <- format_data(exp_df, type="exposure")

# IV-out
outcome_dat <- format_data(out_df, type="outcome")

# Harmonise the exposure and outcome data
mr_dat <- harmonise_data(exposure_dat, outcome_dat, action=1)

# Perform MR
res <- mr(mr_dat)

# MR effect of urate on gout
res_loo <- mr_leaveoneout(mr_dat)

# merge IV x sex P
res_loo <- merge(res_loo, exp_df %>% dplyr::select(SNP, p_int) %>% dplyr::mutate(SNP=tolower(SNP)))

# test for heterogenity


# sex-stratified analysis
stratified <- function(males){
    exp_df <- data.frame()
    out_df <- data.frame()
    for (i in 1:nrow(iv)){
        snp <- paste0("chr", iv$chr[i], "_", as.double(iv$position[i]), "_", iv$nea[i], "_", iv$ea[i])

        # IV-exp effect
        fit <- lm(as.formula(paste0("urate.30880.0.0 ~ ", snp, " + age_at_recruitment.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=dat %>% dplyr::filter(sex.31.0.0 == !!males))
        f_exp <- summary(fit)$fstatistic[1] %>% as.numeric
        b_exp <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(estimate)
        se_exp <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(std.error)

        # IV-outcome
        fit <- glm(as.formula(paste0("gout ~ ", snp, " + age_at_recruitment.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=dat %>% dplyr::filter(sex.31.0.0 == !!males), family="binomial")
        b_out <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(estimate)
        se_out <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(std.error)

        exp_df <- rbind(exp_df, data.frame(
            Phenotype="Urate",
            SNP=snp, 
            beta=b_exp, 
            se=se_exp,
            effect_allele=iv$ea[i],
            other_allele=iv$nea[i],
            units="SD"
        ))
        out_df <- rbind(out_df, data.frame(
            Phenotype="Gout",
            SNP=snp, 
            beta=b_out, 
            se=se_out,
            effect_allele=iv$ea[i],
            other_allele=iv$nea[i],
            units="Log odds"
        ))
    }
    exposure_dat <- format_data(exp_df, type="exposure")

    # IV-out
    outcome_dat <- format_data(out_df, type="outcome")

    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure_dat, outcome_dat, action=1)

    # Perform MR
    res_loo <- mr_leaveoneout(dat)

    return(res_loo)
}