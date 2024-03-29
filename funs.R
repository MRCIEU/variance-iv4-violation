library("dplyr")
library("broom")
#library("varGWASR")

#' Function to simulate genotypes in HWE
#' @param q Recessive/alternative allele frequency
#' @param n_obs Number of observations to return
get_simulated_genotypes <- function(q, n_obs){
  p <- 1 - q
  x <- sample(c(0, 1, 2), n_obs, prob=c(p^2, 2 * p * q, q^2), replace=T)
  return(x)
}

get_est <- function(z, x, y){
    exp_fit <- lm(x ~ z)
    b_exp <- exp_fit %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
    se_exp <- exp_fit %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("std.error") %>% as.numeric
    f_exp <- summary(exp_fit)$fstatistic[['value']]
    out_fit <- lm(y ~ z)
    b_out <- out_fit %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
    se_out <- out_fit %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("std.error") %>% as.numeric
    mr <- TwoSampleMR::mr_wald_ratio(b_exp, b_out, se_exp, se_out, NULL)
    b_mr <- mr$b
    se_mr <- mr$se
    p_mr <- mr$pval
    b_obs <- lm(y ~ x) %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::select("estimate") %>% as.numeric
    fit1 <- varGWASR::model(data.frame(z,x), "z", "x")
    names(fit1) <- paste0(names(fit1), ".zx")
    fit2 <- varGWASR::model(data.frame(z,y), "z", "y")
    names(fit2) <- paste0(names(fit2), ".zy")
    fit <- cbind(fit1, fit2)
    fit$b_mr <- b_mr
    fit$se_mr <- se_mr
    fit$b_obs <- b_obs
    fit$p_mr <- p_mr
    fit$lin_var <- var(x[z==2]) - var(x[z==0])
    fit$f_exp <- f_exp
    return(fit)
}

get_plot <- function(results){
  results$iv1 <- results$f_exp > 10
  p <- ggplot(results, aes(x=-log10(phi_p), y=b_mr, color=iv1)) +
    geom_point(alpha=0.7) +
    facet_grid(~ pb, scale="free_x") +
    scale_colour_grey(start = 0.8,end = 0.2) +
    labs(x="Variance test -log10(P)", y="MR point estimate", color="IV1 (F>10)") +
    theme_classic() +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    geom_hline(yintercept = mean(b_mr), color="grey", linetype="dashed") +
    theme(legend.position="bottom", legend.box.background = element_rect(colour = "black"))
    return(p)
}

#' Test for effect of SNP on outcome variance using the LAD-BF model
#'
#' @param data Dataframe of observations
#' @param x Name of SNP dosage
#' @param y Name of outcome
#' @param covar1 Optional vector of covariates to include in the first-stage model
#' @param covar2 Optional vector of covariates to include in the second-stage model
#' @return Dataframe containing variance effect for SNP=1 (phi_x1) and SNP=2 (phi_x2) with SE and p and F-stat
#' @export
model <- function(data, x, y, covar1=NULL, covar2=NULL){
    if (any(is.na(data))) stop("Dataframe contains NA values")
    if (!x %in% names(data)) stop(paste0(x, " was not in dataframe"))
    if (!y %in% names(data)) stop(paste0(y, " was not in dataframe"))
    if (!all(covar1 %in% names(data))) stop(paste0("Dataframe is missing some of: ", paste0(covar1, collapse=" ")))
    if (!all(covar2 %in% names(data))) stop(paste0("Dataframe is missing some of: ", paste0(covar2, collapse=" ")))
    if (!is.numeric(data[[x]])) stop("Genotype must be numeric")

    # prepare first-stage fit matrix
    if (!is.null(covar1)){
        X <- data %>% dplyr::select(!!x, !!covar1) %>% as.matrix
    } else {
        X <- data %>% dplyr::select(!!x) %>% as.matrix
    }
    # first-stage fit
    fit <- cqrReg::qrfit(X=X, y=data[[y]], tau=.5, method="mm")
    b <- rbind(fit$b, fit$beta)
    # predicted
    X <- cbind(rep(1, nrow(X)), X)
    fitted <- X %*% b
    # residual
    d <- data[[y]] - fitted
    # abs residual
    d <- abs(as.vector(d))
    # second-stage model
    data[[x]] <- as.factor(data[[x]])
    if (!is.null(covar2)){
        X <- data %>% dplyr::select(!!x, !!covar2)
        fit2 <- lm(d ~ ., data=X)
        X <- data %>% dplyr::select(!!covar2)
        fit_null <- lm(d ~ ., data=X)
        p <- anova(fit_null, fit2) %>% broom::tidy(.) %>% dplyr::pull(p.value) %>% dplyr::nth(2)
        f <- anova(fit_null, fit2) %>% broom::tidy(.) %>% dplyr::pull(statistic) %>% dplyr::nth(2)
    } else {
        X <- data %>% dplyr::select(!!x)
        fit2 <- lm(d ~ ., data=X)
        fit_null <- lm(d ~ 1, data=data)
        p <- anova(fit_null, fit2) %>% broom::tidy(.) %>% dplyr::pull(p.value) %>% dplyr::nth(2)
        f <- anova(fit_null, fit2) %>% broom::tidy(.) %>% dplyr::pull(statistic) %>% dplyr::nth(2)
    }

    # deltamethod
    v1 <- car::deltaMethod(fit2, "(2*b0*b1+b1^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"), vcov=sandwich::vcovHC(fit2, type = "HC0"))
    v2 <- car::deltaMethod(fit2, "(2*b0*b2+b2^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"), vcov=sandwich::vcovHC(fit2, type = "HC0"))

    res <- data.frame(
        phi_x1=v1$Estimate,
        se_x1=v1$SE,
        phi_x2=v2$Estimate,
        se_x2=v2$SE,
        phi_f=f,
        phi_p=p
    )
    
    return(res)
}

get_filtered_linker <- function(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="16729") {
    message("Preparing IEU linker")

    if (application == "16729"){
        # load IEU identifiers to application 16729 linker
        linker=fread("/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/16729/2019-04-29/data/ieu_id_linker.csv")
        withdrawn=read.table(pipe( 'ssh bc3 "cat /projects/MRC-IEU/research/data/ukbiobank/phenotypic/applications/16729/withdrawals/meta.withdrawn.20200820.csv"' ), stringsAsFactors=F)$V1
    } else if (application == "15825"){
        # load IEU identifiers to application 15825 linker
        linker <- fread("/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/linker_app15825.csv", col.names=c("appieu", "app15825"))
        withdrawn <- read.table(pipe( 'cat /mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/withdrawals/*' ), stringsAsFactors=F)$V1
    } else {
        abort(paste0("No method for application: ", application))
    }

    # load exclusions
    standard=read.table("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/standard_exclusions/data.combined_recommended.qctools.txt", stringsAsFactors=F)$V1
    minrelated=read.table("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/related/relateds_exclusions/data.minimal_relateds.qctools.txt", stringsAsFactors=F)$V1
    highrelated=read.table("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/related/relateds_exclusions/data.highly_relateds.qctools.txt", stringsAsFactors=F)$V1
    nonwhite=read.table("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/ancestry/data.non_white_british.qctools.txt", stringsAsFactors=F)$V1

    # exclude samples
    linker=linker[which(!(get(paste0("app", application), linker) %in% withdrawn)),]
    
    if (drop_non_white_british){
        linker=linker[which(!(linker$appieu %in% nonwhite)),]
    }  

    if (drop_standard_excl){
        linker=linker[which(!(linker$appieu %in% standard)),]
    }

    if (drop_related){
        linker=linker[which(!(linker$appieu %in% minrelated)),]
        linker=linker[which(!(linker$appieu %in% highrelated)),]
    }

    message(paste0("Left ", nrow(linker), " samples after exclusions."))

    return(linker)
}

get_genetic_principal_components = function(){
    message("Extracting genetic PCs")
    d = fread("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/principal_components/data.pca1-40.qctools.txt", header=FALSE, fill = FALSE, sep=" ")
    names(d) = c("appieu", paste0("PC", seq(1, 40)))
    return(d)
}

get_covariates <- function() {
    cov <- fread("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/standard_covariates/data.covariates.qctools.txt", header=FALSE, fill = FALSE, sep=" ")
    names(cov) = c("appieu", "sex", "chip")
    return(cov)
}

# Extract dosage for variant from BGEN file. Assumes all variants are diploid.
#'
#' @param chrom A single chromosome
#' @param pos Basepair location
#' @param ref Reference allele aka non-effect allele or major allele
#' @param alt Alternative allele aka effect allele or minor allele
#'
#' @return Dataframe of dosages for supplied variant
extract_variant_from_bgen <- function(chrom, pos, ref, alt){
    
    # check inputs
    stopifnot(length(chrom) == 1 && typeof(chrom) == "character")
    stopifnot(length(pos) == 1 && typeof(pos) == "double")
    stopifnot(length(ref) == 1 && typeof(ref) == "character")
    stopifnot(length(alt) == 1 && typeof(alt) == "character")
    # format chrom as expected
    chrom <- gsub("^chr", "", chrom)
    chrom_n <- as.numeric(chrom)
    if (!is.na(chrom_n) && chrom_n < 10){
        chrom <- paste0("0", chrom_n)
    }
    message(paste0("Extracting variant: ", chrom, ":", pos, ref, ">", alt, " from BGEN file"))
    
    # extract variant from bgen
    variant <- bgen.load(paste0("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr", chrom, ".bgen"),
        data.frame(chromosome = chrom, start = pos, end = pos)
    )
    # check we have a single variant at the correct locus
    if (nrow(variant$variants) == 0){
        stop("Variant not found")
    }
    if (nrow(variant$variants) > 1){
        stop("Variant is multiallelic and not currently supported")
    }
    if (variant$variants$chromosome != chrom){
        stop("Wrong variant returned")
    }
    if (variant$variants$position != pos){
        stop("Wrong variant returned")
    }
    # harmonise
    if (variant$variants$allele0 == ref && variant$variants$allele1 == alt){
        # convert dosage for each genotype to copies of alt allele
        dosage <- as.data.frame(
            apply(variant$data, 1, function(data) { return(data[,1]*0 + data[,2]*1 + data[,3]*2) })
        )
    } else if (variant$variants$allele1 == ref && variant$variants$allele0 == alt){
        # convert dosage for each genotype to copies of alt allele
        dosage <- as.data.frame(
            apply(variant$data, 1, function(data) { return(data[,1]*2 + data[,2]*1 + data[,3]*0) })
        )
    } else {
        stop("Alleles do not match input")
    }
    # tidy up
    names(dosage) <- paste0("chr", gsub("^0", "", chrom), "_", pos, "_", ref, "_", alt)
    dosage$appieu <- row.names(dosage)
    row.names(dosage) <- NULL
    return(dosage)
}

get_variants <- function(trait, chrs=seq(1,22)){
    # load vGWAS & SNP stats; QC loci
    data <- data.frame()

    for (chr in chrs){
        message(paste0("loading chr", chr))
        if (chr < 10){
            gwas <- fread(paste0("/user/home/ml18692/projects/varGWAS-ukbb-biomarkers/data/", trait, ".vgwas.chr0", chr, ".txt"))
            snp_stats <- fread(paste0("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/snp-stats/data.chr0", chr, ".snp-stats"), skip=15)
        } else {
            gwas <- fread(paste0("/user/home/ml18692/projects/varGWAS-ukbb-biomarkers/data/", trait, ".vgwas.chr", chr, ".txt"))
            snp_stats <- fread(paste0("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/snp-stats/data.chr", chr, ".snp-stats"), skip=15)
        }

        # drop multiallelics by rsid
        counts <- table(snp_stats$rsid)
        ma <- as.data.frame(counts[which(counts>1)])
        snp_stats <- snp_stats[!snp_stats$rsid %in% ma$Var1]

        # drop multiallelics by position
        counts <- table(snp_stats$position)
        ma <- as.data.frame(counts[which(counts>1)])
        snp_stats <- snp_stats[!snp_stats$position %in% ma$Var1]
        
        # exclude MAF < q
        #snp_stats <- snp_stats[which(snp_stats$minor_allele_frequency > 0.05)]
        
        # exclude HWE violations
        snp_stats <- snp_stats[which(snp_stats$HW_exact_p_value > 1e-5)]

        # exclude high missingness
        snp_stats <- snp_stats[which(snp_stats$missing_proportion < 0.05)]

        # exclude low imputation quality
        snp_stats <- snp_stats[which(snp_stats$info > 0.3)]

        # drop HLA region with 10Mb pad
        snp_stats <- snp_stats[!(snp_stats$chromosome == 6 & snp_stats$position >= (28477797-5000000) & snp_stats$position <= 33448354+5000000),]

        # drop vGWAS failed rows
        gwas <- gwas %>% dplyr::filter(phi_p != -1)

        # merge/filter vGWAS
        snp_stats$key <- paste0(snp_stats$chromosome, "_", snp_stats$position, "_", snp_stats$alleleA, "_", snp_stats$alleleB)
        snp_stats$rsid <- NULL
        gwas$key <- paste0(gwas$chr, "_", gwas$pos, "_", gwas$oa, "_", gwas$ea)
        gwas <- merge(gwas, snp_stats, "key")

        # store results
        data <- rbind(data, gwas)
    }

    return(data)
}

get_r2 <- function(b, se, eaf, n){
    r2 <- 2*b^2*(eaf)*(1-eaf)/ (2*b^2*(eaf)*(1-eaf) + (se)^2*2*n*eaf*(1-eaf))
    return(r2)
}

# n_snps = number of snps in the model. 1=univariable; n=multivarible
get_f <- function(r2, n, n_snps){
    f <- r2*(n-1-n_snps)/((1- r2)*n_snps)
    return(f)
}