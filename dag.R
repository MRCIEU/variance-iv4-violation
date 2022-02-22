library("dagitty")
library("ggdag")
library("ggpubr")
set.seed(123)

# positive control DAG: smoking heaviness-by-status on lung function
dag <- dagify(Y ~ U + XU, XU ~ X + U, ZU ~ Z + U, X ~ Z + ZU, exposure="X", outcome="Y", labels=c("Z" = "Variant", "Y" = "Lung\nfunction", "U" = "Smoking\nstatus", "XU" = "Smoking\nheaviness-by-status", "X" = "Smoking\nheaviness", "ZU" = "Genotype-by-status"))
dag <- dag_paths(dag)

# positive
pos <- ggdag(dag, use_labels = "label", layout = "circle") + theme_dag_blank()

ggsave("dag.pdf", pos)