library(GEOquery)
library(ggplot2)
library(limma)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(patchwork)
library(ggrepel)


#STEP 1 — Load GEO dataset
gse <- getGEO("GSE20916", GSEMatrix = TRUE)

expr <- gse[[1]]


#STEP 2 — Extract expression matrix
expr_matrix <- exprs(expr)

metadata <- pData(expr)

#STEP 3 — Create biological groups
group <- ifelse(grepl("normal", metadata$title, ignore.case = TRUE),
                "Normal",
                ifelse(grepl("carcinoma", metadata$title, ignore.case = TRUE),
                       "Cancer", "Adenoma"))

names(group) <- colnames(expr_matrix)


#STEP 4 — Differential Expression
expr_matrix <- na.omit(expr_matrix)

#Filtracji low-expression genes
keep <- rowMeans(expr_matrix) > 5
expr_matrix <- expr_matrix[keep, ]

#Differential expression analysis using LIMMA
design <- model.matrix(~0 + group)

colnames(design) <- levels(factor(group))

fit <- lmFit(expr_matrix, design)

contrast.matrix <- makeContrasts(
  Adenoma_vs_Normal = Adenoma - Normal,
  Cancer_vs_Normal = Cancer - Normal,
  Adenoma_vs_Cancer = Adenoma - Cancer,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Results for different groups
results_CN <- topTable(fit2, coef = "Cancer_vs_Normal", number = Inf)

results_AN <- topTable(fit2, coef = "Adenoma_vs_Normal", number = Inf)

results_AC <- topTable(fit2, coef = "Adenoma_vs_Cancer", number = Inf)


#STEP 5 - Probe-to-gene annotation
mapping <- fData(expr)[, c("ID", "Gene Symbol")]
colnames(mapping) <- c("probe_id", "gene_symbol")

expr_df <- as.data.frame(expr_matrix)
expr_df$probe_id <- rownames(expr_df)

expr_annot <- merge(mapping, expr_df, by = "probe_id")

expr_annot <- expr_annot[expr_annot$gene_symbol != "" & !is.na(expr_annot$gene_symbol), ]

expr_values <- expr_annot[, !(names(expr_annot) %in% c("probe_id", "gene_symbol"))]

expr_annot$variance <- apply(expr_values, 1, var, na.rm = TRUE)

expr_annot <- expr_annot[order(expr_annot$gene_symbol, -expr_annot$variance),]

expr_annot <- expr_annot[!duplicated(expr_annot$gene_symbol),]

expr_annot$variance <- NULL

#Gene matrix
rownames(expr_annot) <- expr_annot$gene_symbol

expr_matrix_genes <- expr_annot[, !(names(expr_annot) %in% c("probe_id", "gene_symbol"))]

#Mapping LIMMA results with genes
results_CN$Gene <- mapping$gene_symbol[
  match(rownames(results_CN), mapping$probe_id)
]

results_AN$Gene <- mapping$gene_symbol[
  match(rownames(results_AN), mapping$probe_id)
]

results_AC$Gene <- mapping$gene_symbol[
  match(rownames(results_AC), mapping$probe_id)
]


#STEP 6 — Volcano plot
#Cancer vs Normal
results_CN$status <- "Non-significant"

results_CN$status[results_CN$adj.P.Val < 0.05 &
                                  results_CN$logFC > 1] <- "Upregulated"

results_CN$status[results_CN$adj.P.Val < 0.05 &
                                  results_CN$logFC < -1] <- "Downregulated"

top_labels_CN <- results_CN[order(results_CN$adj.P.Val),][1:15, ]

p1.1 <- ggplot(results_CN,
               aes(logFC, -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.7, size = 1.5) + 
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", 
                                "Non-significant" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + geom_text_repel(
    data = top_labels_CN, aes(label = Gene), size = 3, max.overlaps = 20) +
  theme_minimal() + ggtitle("Cancer vs Normal")

#Adenoma vs Normal
results_AN$status <- "Non-significant"

results_AN$status[results_AN$adj.P.Val < 0.05 &
                                   results_AN$logFC > 1] <- "Upregulated"

results_AN$status[results_AN$adj.P.Val < 0.05 &
                                   results_AN$logFC < -1] <- "Downregulated"

top_labels_AN <- results_AN[order(results_AN$adj.P.Val),][1:15, ]

p1.2 <- ggplot(results_AN,
               aes(logFC, -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.7, size = 1.5) + 
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue",
                                "Non-significant" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + geom_text_repel(
    data = top_labels_AN, aes(label = Gene), size = 3, max.overlaps = 20) +
  theme_minimal() + ggtitle("Adenoma vs Normal")

#Adenoma vs Cancer
results_AC$status <- "Non-significant"

results_AC$status[results_AC$adj.P.Val < 0.05 &
                                   results_AC$logFC > 1] <- "Upregulated"

results_AC$status[results_AC$adj.P.Val < 0.05 &
                                   results_AC$logFC < -1] <- "Downregulated"

top_labels_AC <- results_AC[order(results_AC$adj.P.Val),][1:15, ]

p1.3 <- ggplot(results_AC,
               aes(logFC, -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.7, size = 1.5) + 
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", 
                                "Non-significant" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + geom_text_repel(
    data = top_labels_AC, aes(label = Gene), size = 3, max.overlaps = 20) +
  theme_minimal() + ggtitle("Adenoma vs Cancer")


#STEP 7 - PCA
pca <- prcomp(t(expr_matrix), scale. = TRUE)

pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Group = group)

p2 <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal()


#STEP 8 - Top genes
top_genes <- results_CN[order(results_CN$adj.P.Val, -abs(results_CN$logFC)),]

top_genes <- na.omit(top_genes$Gene)

top_genes <- top_genes[1:min(30, length(top_genes))]

top_genes <- intersect(top_genes, rownames(expr_matrix_genes))


#STEP 9 - Heatmap
heatmap_data <- expr_matrix_genes[top_genes, , drop = FALSE]

heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd, na.rm = TRUE) > 0, ]

heatmap_scaled <- t(scale(t(heatmap_data)))

heatmap_scaled[!is.finite(heatmap_scaled)] <- 0

names(group) <- colnames(expr_matrix_genes)

ordered_samples <- names(group)[order(factor(group, levels = c("Normal", "Adenoma", "Cancer")))]

heatmap_scaled <- heatmap_scaled[, ordered_samples, drop = FALSE]

annotation_df <- data.frame(Group = group[ordered_samples])

rownames(annotation_df) <- ordered_samples

#Plotting
my_colors <- colorRampPalette(rev(brewer.pal(7, "RdBu")))(100)

p3.2 <- pheatmap(heatmap_scaled,color = my_colors, annotation_col = annotation_df, 
                 cluster_cols = FALSE, cluster_rows = TRUE, show_colnames = FALSE, 
                 show_rownames = TRUE, fontsize_row = 8, border_color = NA, 
                 main = "Gene-level Heatmap: Normal → Adenoma → Cancer")


#STEP 10 - Functional Enrichment Analysis
#Data preparation

sig_genes_CN <- results_CN[results_CN$adj.P.Val < 0.01 &
                             abs(results_CN$logFC) > 2,]
sig_genes_AN <- results_AN[results_AN$adj.P.Val < 0.01 &
                             abs(results_AN$logFC) > 2,]
sig_genes_AC <- results_AC[results_AC$adj.P.Val < 0.01 & 
                             abs(results_AC$logFC) > 2,]

genes_CN <- sig_genes_CN$Gene
genes_AN <- sig_genes_AN$Gene
genes_AC <- sig_genes_AC$Gene


genes_CN <- unique(genes_CN)
genes_CN <- genes_CN[!is.na(genes_CN) & genes_CN != ""]

genes_AN <- unique(genes_AN)
genes_AN <- genes_AN[!is.na(genes_AN) & genes_AN != ""]

genes_AC <- unique(genes_AC)
genes_AC <- genes_AC[!is.na(genes_AC) & genes_AC != ""]


gene_entrez_CN <- bitr(genes_CN, fromType = "SYMBOL", toType = "ENTREZID", 
                       OrgDb = org.Hs.eg.db)
gene_entrez_AN <- bitr(genes_AN, fromType = "SYMBOL", toType = "ENTREZID", 
                       OrgDb = org.Hs.eg.db)
gene_entrez_AC <- bitr(genes_AC, fromType = "SYMBOL", toType = "ENTREZID", 
                       OrgDb = org.Hs.eg.db)

#GO enrichment analysis
go_results_CN <- enrichGO(
  gene = gene_entrez_CN$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)
p4.1 <- dotplot(go_results_CN, showCategory = 10) +
  ggtitle("GO Biological Processes: Cancer vs Normal")

go_results_AN <- enrichGO(
  gene = gene_entrez_AN$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)
p4.2 <-  dotplot(go_results_AN, showCategory = 10) +
  ggtitle("GO Biological Processes: Adenoma vs Normal")

go_results_AC <- enrichGO(
  gene = gene_entrez_AC$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)
p4.3 <- dotplot(go_results_AC, showCategory = 10) +
  ggtitle("GO Biological Processes: Adenoma vs Cancer")

#KEGG enrichment
kegg_results_CN <- enrichKEGG(
  gene = gene_entrez_CN$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH"
)
p5.1 <- dotplot(kegg_results_CN, showCategory = 10) +
  ggtitle("KEGG Pathway Enrichment: Cancer vs Normal")

kegg_results_AN <- enrichKEGG(
  gene = gene_entrez_AN$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH"
)
p5.2 <- dotplot(kegg_results_AN, showCategory = 10) +
  ggtitle("KEGG Pathway Enrichment: Adenoma vs Normal")

kegg_results_AC <- enrichKEGG(
  gene = gene_entrez_AC$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH"
)
p5.3 <- dotplot(kegg_results_AC, showCategory = 10) +
  ggtitle("KEGG Pathway Enrichment: Adenoma vs Cancer")


#STEP 11 - Network plot
go_results_CN <- simplify(go_results_CN)
go_results_CN <- pairwise_termsim(go_results_CN)

go_results_AN <- simplify(go_results_AN)
go_results_AN <- pairwise_termsim(go_results_AN)

go_results_AC <- simplify(go_results_AC)
go_results_AC <- pairwise_termsim(go_results_AC)

p6.1 <- emapplot(go_results_CN,
                       showCategory = 10) + 
  ggtitle("GO Biological Processes: Cancer vs Normal")
p6.2 <- emapplot(go_results_AN, showCategory = 10) +
  ggtitle("GO Biological Processes: Adenoma vs Normal")
p6.3 <- emapplot(go_results_AC, showCategory = 10) +
  ggtitle("GO Biological Processes: Adenoma vs Cancer")


#STEP 12 - Saving figures
ggsave("figures/volcano_plots.png", p1.1 + p1.2 + p1.3, width = 14, height = 5)
ggsave("figures/pca_plot.png", p2, width = 6, height = 5)
ggsave("figures/go_dotplots.png", p4.1 + p4.2 + p4.3, width = 14, height = 5)
ggsave("figures/kegg_dotplots.png", p5.1 + p5.2 + p5.3, width = 14, height = 5)
ggsave("figures/go_network.png", p6.1 + p6.2 + p6.3, width = 14, height = 5)
png("figures/heatmap.png", width = 1000, height = 900)

pheatmap(heatmap_scaled, color = my_colors, annotation_col = annotation_df,
         cluster_cols = FALSE, cluster_rows = TRUE, show_colnames = FALSE,
         show_rownames = TRUE, fontsize_row = 8, border_color = NA,
         main = "Gene-level Heatmap: Normal → Adenoma → Cancer")
dev.off()
