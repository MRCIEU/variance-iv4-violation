library("TwoSampleMR")
library("ggplot2")
set.seed(123)

# could use testosterone on CVD outcomes

# List available GWASs
ao <- available_outcomes()

# Get instruments
exposure_dat <- extract_instruments("ieu-a-835")

# Get effects of instruments on outcome
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-38")

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)
stopifnot(all(dat$effect_allele.exposure == dat$effect_allele.outcome))
stopifnot(all(dat$other_allele.exposure == dat$other_allele.outcome))

# Perform MR
res <- mr(dat)
res_single <- mr_singlesnp(dat)

# estimate delta
res_single$egger.delta <- res_single$b - res %>% dplyr::filter(method == "MR Egger") %>% dplyr::pull(b)
res_single$wmedian.delta <- res_single$b - res %>% dplyr::filter(method == "Weighted median") %>% dplyr::pull(b)
res_single$ivw.delta <- res_single$b - res %>% dplyr::filter(method == "Inverse variance weighted") %>% dplyr::pull(b)
res_single$smode.delta <- res_single$b - res %>% dplyr::filter(method == "Simple mode") %>% dplyr::pull(b)
res_single$wmode.delta <- res_single$b - res %>% dplyr::filter(method == "Weighted mode") %>% dplyr::pull(b)

# drop non-SNP
res_single <- res_single %>% dplyr::filter(grepl("^rs", res_single$SNP))

# append alleles
res_single <- merge(res_single, dat %>% dplyr::select(SNP, effect_allele.exposure, other_allele.exposure) %>% dplyr::rename(effect_allele="effect_allele.exposure", other_allele = "other_allele.exposure"), "SNP")

# add key
res_single$key <- paste0(res_single$SNP, "_", res_single$other_allele, "_", res_single$effect_allele)

# IV variance effects
ve <- fread("/Users/ml18692/projects/variance-iv4-violation/bmi_iv_var.txt")
ve$key <- paste0(ve$term, "_", ve$nea, "_", ve$ea)
ve$term <- NULL
ve$nea <- NULL
ve$ea <- NULL
BMI <- ve %>% dplyr::filter(trait == "body_mass_index.21001.0.0") %>% dplyr::select(!trait)
names(BMI) <- paste0(names(BMI), ".BMI")
SBP <- ve %>% dplyr::filter(trait == "systolic_blood_pressure_automated_reading.4080.0.0") %>% dplyr::select(!trait)
names(SBP) <- paste0(names(SBP), ".SBP")

# merge
res_single <- merge(res_single, BMI, by.x="key", by.y="key.BMI")
res_single <- merge(res_single, SBP, by.x="key", by.y="key.SBP")

# wide to long
long <- res_single %>% dplyr::select(-egger.delta, -wmedian.delta, -ivw.delta, -smode.delta, -wmode.delta)

res_single_long <- data.frame()

long$delta <- res_single$egger.delta
long$method <- "MR Egger"
res_single_long <- rbind(res_single_long, long)

long$delta <- res_single$wmedian.delta
long$method <- "Weighted Median"
res_single_long <- rbind(res_single_long, long)

long$delta <- res_single$ivw.delta
long$method <- "IVW"
res_single_long <- rbind(res_single_long, long)

long$delta <- res_single$smode.delta
long$method <- "Simple Mode"
res_single_long <- rbind(res_single_long, long)

long$delta <- res_single$wmode.delta
long$method <- "Weighted Mode"
res_single_long <- rbind(res_single_long, long)

# plot
ggplot(res_single_long, aes(x=-log10(phi_p.BMI), y=abs(delta))) +
    geom_point() +
    facet_grid(~method) +
    theme_classic()
ggplot(res_single_long, aes(x=phi_x1.BMI, y=phi_x2.BMI, size=abs(delta))) +
    geom_point() +
    facet_grid(~method) +
    theme_classic()
ggplot(res_single_long, aes(x=-log10(phi_p.SBP), y=delta)) +
    geom_point() +
    facet_grid(~method) +
    theme_classic()
ggplot(res_single_long, aes(x=phi_x1.SBP, y=phi_x2.SBP, size=abs(delta))) +
    geom_point() +
    facet_grid(~method) +
    theme_classic()