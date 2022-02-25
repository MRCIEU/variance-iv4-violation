load("/user/home/ml18692/projects/varGWAS-ukbb-biomarkers/data/pheno.RData")

library("TwoSampleMR")
library("ieugwasr")
library("ggplot2")
library("viridis")
library("varGWASR")
source("funs.R")
options(ieugwasr_api="http://64.227.44.193:8006/")
set.seed(123)

# load urate IV
exposure_dat <- extract_instruments("ukb-d-30880_irnt")

# load vGWAS
vgwas <- get_variants("urate.30880.0.0")

# Merge instruments with vQTL P
exposure_dat <- merge(
    exposure_dat,
    vgwas %>% dplyr::select(rsid, phi_p),
    by.x="SNP",
    by.y="rsid"
)

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")

# load covariates
pc <- get_genetic_principal_components()

# merge data
dat <- merge(linker, pc, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")
dat$urate.30880.0.0 <- dat$urate.30880.0.0 / sd(dat$urate.30880.0.0, na.rm=T)

# select SNPs and extract from UKBB
for (i in 1:nrow(exposure_dat)){
    dosage <- tryCatch(
        {
            extract_variant_from_bgen(exposure_dat$chr.exposure[i], as.double(exposure_dat$pos.exposure[i]), exposure_dat$other_allele.exposure[i], exposure_dat$effect_allele.exposure[i])
        },
        error=function(cond) {
            return(NULL)
        }
    )
    if (!is.null(dosage)){
        dat <- merge(dat, dosage, "appieu")
    }
}

# test for IV effect on urate variance w/wo adjusting for by-sex effect
results <- data.frame()
for (i in 1:nrow(exposure_dat)){
    # define SNP & int term
    snp <- paste0("chr", exposure_dat$chr.exposure[i], "_", as.double(exposure_dat$pos.exposure[i]), "_", exposure_dat$other_allele.exposure[i], "_", exposure_dat$effect_allele.exposure[i])
    dat$gxs <- dat[[snp]] * dat$sex.31.0.0
    # fit wo by-sex effect
    covar <- c("age_at_recruitment.21022.0.0", "sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
    fit0 <- varGWASR::model(dat %>% dplyr::select(!!covar, snp, "urate.30880.0.0") %>% tidyr::drop_na(.), snp, "urate.30880.0.0", covar1 = covar, covar2 = covar)
    # fit w by-sex effect
    covar <- c(covar, "gxs")
    fit1 <- varGWASR::model(dat %>% dplyr::select(!!covar, snp, "urate.30880.0.0") %>% tidyr::drop_na(.), snp, "urate.30880.0.0", covar1 = covar, covar2 = covar)
    # store result
    names(fit0) <- paste0(names(fit0), ".m")
    names(fit1) <- paste0(names(fit1), ".i")
    result <- cbind(fit0, fit1)
    result$snp <- snp
    results <- rbind(results, result)
}

# plot -log10(P) variance effect w/wo adjusting for by-sex effect
