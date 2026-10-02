gene_names <- c("BRCA1", "TP53", "EGFR", "MYC", "PTEN", "KRAS")
control <- c(5.2, 7.8, 3.1, 9.4, 6.0, 4.5)
treated <- c(8.9, 7.6, 6.7, 12.1, 2.3, 9.8)
chromosome <- c("17", "17", "7", "8", "10", "12")

names(control) <- gene_names
names(treated) <- gene_names

fold_change <- treated / control

classification <- ifelse(fold_change > 1.2, "upregulated",
                         ifelse(fold_change < 0.8, "downregulated", "stable"))
classification_factor <- factor(classification, 
                                levels = c("upregulated", "downregulated", "stable"))
gene_data <- list(
  names = gene_names,
  stats = list(
    control = control,
    treated = treated
  ),
  source = "Example gene expression dataset"
)

gene_data$names

gene_data$date_analyzed <- "2026-10-02"

gene_table <- data.frame(
  gene_names = gene_names,
  control = control,
  treated = treated,
  chromosome = chromosome,
  fold_change = fold_change,
  classification = classification
)

gene_table[gene_table$chromosome == "17", ]

gene_table$log2_fc <- log2(gene_table$fold_change)
gene_table <- gene_table[order(gene_table$fold_change, decreasing = TRUE), ]
write.csv(gene_table, "gene_table_sorted.csv", row.names = FALSE)


png("expression_plot.png", width = 800, height = 600)
plot(
  gene_table$control,
  gene_table$treated,
  xlab = "Control Expression",
  ylab = "Treated Expression",
  main = "Gene Expression: Control vs Treated",
  col = ifelse(
    gene_table$classification == "upregulated", "pink",
    ifelse(
      gene_table$classification == "downregulated", "purple",
      "gray"
    )
  ),
  pch = 19
)

dev.off()

