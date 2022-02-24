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

# get sex-stratified urate GWAS from Neale et al
#download.file("https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30880_irnt.gwas.imputed_v3.female.varorder.tsv.bgz", "data/30880_irnt.gwas.imputed_v3.female.varorder.txt.gz")
female <- fread("data/30880_irnt.gwas.imputed_v3.female.varorder.txt.gz")
#download.file("https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30880_irnt.gwas.imputed_v3.male.varorder.tsv.bgz", "data/30880_irnt.gwas.imputed_v3.male.varorder.txt.gz")
male <- fread("data/30880_irnt.gwas.imputed_v3.male.varorder.txt.gz")
#download.file("https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz", "data/variants.txt.gz")
anno <- fread("data/variants.txt.gz", select=c("variant", "rsid", "chr", "pos", "ref", "alt", "AF"))
female <- merge(female, anno, "variant")
male <- merge(male, anno, "variant")

# extract tophits
female <- female %>% 
    dplyr::mutate(key = paste0(chr, ":", pos)) %>% 
    dplyr::group_by(key) %>% 
    dplyr::filter(n() == 1) %>% 
    dplyr::filter(pval < 5e-8)
male <- male %>% 
    dplyr::mutate(key = paste0(chr, ":", pos)) %>% 
    dplyr::group_by(key) %>% 
    dplyr::filter(n() == 1) %>% 
    dplyr::filter(pval < 5e-8)

# map to tsmr
female_tsmr <- female %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(Phenotype="Urate",units="SD") %>% 
    dplyr::select(Phenotype, rsid, beta, se, alt, ref, AF, pval, units, n_complete_samples, chr, pos) %>% 
    dplyr::rename(SNP="rsid", effect_allele="alt", other_allele="ref", eaf="AF", samplesize="n_complete_samples")
female_tsmr <- format_data(female_tsmr, type="exposure")
male_tsmr <- male %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(Phenotype="Urate",units="SD") %>% 
    dplyr::select(Phenotype, rsid, beta, se, alt, ref, AF, pval, units, n_complete_samples, chr, pos) %>% 
    dplyr::rename(SNP="rsid", effect_allele="alt", other_allele="ref", eaf="AF", samplesize="n_complete_samples")
male_tsmr <- format_data(male_tsmr, type="exposure")

# clump exposure SNPs
female_tsmr <- clump_data(female_tsmr)
male_tsmr <- clump_data(male_tsmr)

# select SNPs and extract from UKBB
iv <- unique(rbind(
    female_tsmr %>% dplyr::select(chr.exposure, pos.exposure, other_allele.exposure, effect_allele.exposure) %>% dplyr::rename(chr="chr.exposure", position="pos.exposure", nea="other_allele.exposure", ea="effect_allele.exposure"),
    male_tsmr %>% dplyr::select(chr.exposure, pos.exposure, other_allele.exposure, effect_allele.exposure) %>% dplyr::rename(chr="chr.exposure", position="pos.exposure", nea="other_allele.exposure", ea="effect_allele.exposure")
))
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

get_estimates <- function(iv, sex){
    # test for SNP effects
    exp_df <- data.frame()
    out_df <- data.frame()
    # stratify by sex
    datstrat <- dat %>% dplyr::filter(sex.31.0.0 == !!sex)
    # SD scale outcomes
    datstrat$urate.30880.0.0 <- datstrat$urate.30880.0.0 / sd(datstrat$urate.30880.0.0, na.rm=T)
    for (i in 1:nrow(iv)){
        snp <- paste0("chr", iv$chr.exposure[i], "_", as.double(iv$pos.exposure[i]), "_", iv$other_allele.exposure[i], "_", iv$effect_allele.exposure[i])

        if (!snp %in% names(dat)){
            next
        }

        # IV-exp mean
        fit <- lm(as.formula(paste0("urate.30880.0.0 ~ ", snp, " + age_at_recruitment.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=datstrat)
        f_exp <- summary(fit)$fstatistic[1] %>% as.numeric
        b_exp <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(estimate)
        se_exp <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(std.error)

        # IV-exp variance
        covar <- c("age_at_recruitment.21022.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
        v_exp <- model(datstrat %>% dplyr::select(urate.30880.0.0, !!snp, age_at_recruitment.21022.0.0, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% tidyr::drop_na(.), snp, "urate.30880.0.0", covar1 = covar, covar2 = covar)

        # IV-outcome
        fit <- glm(as.formula(paste0("gout ~ ", snp, " + age_at_recruitment.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=datstrat, family="binomial")
        b_out <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(estimate)
        se_out <- fit %>% tidy %>% dplyr::filter(grepl("chr",term)) %>% dplyr::pull(std.error)

        exp_df <- rbind(exp_df, data.frame(
            Phenotype="urate.30880.0.0",
            SNP=snp, 
            beta=b_exp, 
            se=se_exp,
            effect_allele=iv$ea[i],
            other_allele=iv$nea[i],
            units="SD",
            f_exp, phi_p=v_exp$phi_p
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
}