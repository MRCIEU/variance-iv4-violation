library("TwoSampleMR")
library("ggplot2")
set.seed(123)

mr_scatter_var_plot <- function(mr_results, dat) {
    # restict to IVW
    mr_results <- mr_results %>% dplyr::filter(method=="Inverse variance weighted")

    # add intercept
    mr_results$a <- 0
    
    # drop failing SNPs
    dat <- subset(dat, mr_keep)

    # flip betas
    index <- dat$beta.exposure < 0
    dat$beta.exposure[index] <- dat$beta.exposure[index] * -1
    dat$beta.outcome[index] <- dat$beta.outcome[index] * -1
    
    # plot
    p <- ggplot2::ggplot(data=dat, ggplot2::aes(x=beta.exposure, y=beta.outcome, size=-log10(phi_p), color=gene, shape=inv)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(data=mr_results, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE) +
        ggplot2::scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
        ggplot2::labs(colour="MR Test", x=paste("SNP effect on", dat$exposure[1]), y=paste("SNP effect on", dat$outcome[1])) +
        ggplot2::theme(legend.position="top", legend.direction="vertical") +
        ggplot2::guides(colour=ggplot2::guide_legend(ncol=2))
    
    return(p)
}

# List available GWASs
ao <- available_outcomes()

# Get instruments
exposure_dat <- extract_instruments("ukb-d-30880_irnt")

# Get effects of instruments on outcome
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST001790")

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

# load IV-LDL variance effects
ve <- data.frame()
for (file in list.files("/Users/ml18692/projects/varGWAS-ukbb-biomarkers/data/urate")){
    message(paste0("working on:", file))
    ve <- rbind(ve, fread(paste0("/Users/ml18692/projects/varGWAS-ukbb-biomarkers/data/urate/", file)))
}
dat <- merge(dat, ve %>% dplyr::select(rsid, phi_x1, phi_x2, phi_p) %>% dplyr::rename(SNP="rsid"), "SNP")

# append SLC2A9 status
dat <- dat %>% dplyr::mutate(gene = case_when(chr.exposure == "4" & pos.exposure >= 9772777 - 500000 & pos.exposure <= 10056560 + 500000 ~ TRUE, TRUE ~ FALSE))

# scatter plot
dat$inv <- dat$phi_x2 < 0
mr_scatter_var_plot(res, dat)

# BMI -> SBP
# GIANT -> 