library("epigraphdb")
library("data.table")
library("ReactomeContentService4R")
library("broom")
library("TwoSampleMR")
library("ieugwasr")
library("ggplot2")
source("funs.R")
set.seed(124)

# Get instruments for LDL-c
#exposure_dat <- extract_instruments("ukb-d-30780_irnt")
#write.csv(exposure_dat, file="ldl.csv")
exposure_dat <- fread("ldl.csv")

# Get variance effect for LDL-c instruments
vgwas <- get_variants("ldl_direct.30780.0.0")

# Merge instruments with vQTL P
exposure_dat <- merge(
    exposure_dat,
    vgwas %>% dplyr::select(rsid, phi_p),
    by.x="SNP",
    by.y="rsid"
)

# define P-threshold using quantile
p_var <- quantile(exposure_dat$phi_p, .75)

# annotate IVs with gene using eqtlgen data
#genes <- phewas(exposure_dat$SNP, batch=c("eqtl-a"))
#write.csv(genes, file="genes.csv")
genes <- fread("genes.csv")

# Test for enrichment of pathways

# Map genes to proteins
proteins <- query_epigraphdb(
  route = "/mappings/gene-to-protein",
  params = list(by_gene_id = TRUE, gene_id_list = genes %>% dplyr::pull(trait) %>% unique),
  mode = "table",
  method = "POST"
)

# Map proteins to pathways
proteins_uniprot_ids <- proteins$protein.uniprot_id %>% I()
pathway_df <- query_epigraphdb(
  route = "/protein/in-pathway",
  params = list(uniprot_id_list = proteins_uniprot_ids),
  mode = "table",
  method = "POST"
)

# restructure pathways
pathway_df <- pathway_df %>% 
  dplyr::group_by(uniprot_id) %>% 
  dplyr::summarize(pathway_reactome_id=unlist(pathway_reactome_id)) %>%
  dplyr::ungroup()

# merge data
pp <- merge(proteins, pathway_df, by.x="protein.uniprot_id", by.y="uniprot_id")

# unique pathways
pathways <- unlist(pp$pathway_reactome_id) %>% unique

# Test for enrichment of MR evidence for each pathway
results <- data.frame()
for (reactome_id in pathways){
  print(paste0("Working on pathway:", reactome_id))
  
  # Gene list for proteins in pathway
  pathway_genes <- pp %>% dplyr::filter(pathway_reactome_id == !!reactome_id) %>% dplyr::pull(gene.name) %>% unique

  if (length(pathway_genes) == 0){
    next
  }

  # check if genes are in pathway
  cis_mr$gene_in_pathway <- cis_mr$HGNC.symbol.protein %in% pathway_genes

  # 2x2 table for genes in pathway by outcome
  res <- cis_mr %>%
    dplyr::filter(outcome == !!outcome) %>%
    dplyr::summarize(
        associated_pathway=sum(pval < 0.05 & HGNC.symbol.protein %in% pathway_genes),
        associated_not_pathway=sum(pval < 0.05 & !(HGNC.symbol.protein %in% pathway_genes)),
        not_associated_pathway=sum(pval >= 0.05 & HGNC.symbol.protein %in% pathway_genes),
        not_associated_not_pathway=sum(pval >= 0.05 & !(HGNC.symbol.protein %in% pathway_genes)),
        n_tested_genes=sum(HGNC.symbol.protein %in% pathway_genes)
      ) %>%
    dplyr::summarize(cbind(n_tested_genes, fisher.test(matrix(c(associated_pathway, associated_not_pathway, not_associated_pathway, not_associated_not_pathway), nrow = 2), alternative = "greater") %>% tidy)) %>%
    dplyr::mutate(reactome_id=!!reactome_id, n_pathway_genes=length(!!pathway_genes))

  # store results
  results <- rbind(results, res)
}

# annotate with reactome name and plot
e$name <- apply(e, 1, function(x) query(id = x[["reactome_id"]])$displayName)
e$name <- factor(e$name)
e$name <- factor(e$name, levels=rev(levels(e$name)))
e$outcome <- factor(e$outcome, levels=c("AS", "AIS", "LAS", "CES", "SVS"))
pdf("epathways.pdf", width=10)
ggplot(e, aes(x=name, y=-log10(p.value), fill=outcome)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position="bottom") + 
  labs(x="Pathway")
dev.off()