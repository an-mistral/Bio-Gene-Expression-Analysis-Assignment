# Chunk 1
#_______________________________________________

# Specify the path to the archive
downloads_path = "/Users/anastasiiadeviataieva/Desktop/UCD Autumn/BioPrin and CellOrg/Gene expression analysis/Gene Expression Analysis"
file_path = paste(downloads_path, "brca_tcga_pan_can_atlas_2018.tar", sep = "/")

# Unpack the archive into the same folder
untar(file_path)

# Check the contents of the folder after unpacking
# list.files(path = downloads_path)

# Path to the unarchived file
folder_path = paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/" )

# Paths to files
data_patient_path = paste(folder_path,"data_clinical_patient.txt", sep = "/")
path_RNA = paste(folder_path,"data_mrna_seq_v2_rsem.txt", sep = "/")
path_cna = paste(folder_path,"data_cna.txt", sep = "/")

# Reading files
data_Rnaseq = read.delim(path_RNA)
data_patient = read.delim(data_patient_path)
data_cna = read.delim(path_cna)

data_patient = data_patient[5:dim(data_patient)[1],] # we will skip 5 rows of column descriptions. 


# Chunk 2
#_______________________________________________
# Finding patient identification data
cna_ids <- substr(colnames(data_cna)[3:ncol(data_cna)], 1, 12)
patient_ids <- substr(gsub("-", "\\.", data_patient[, 1]), 1, 12)
# Matching patient identification data
common_ids <- intersect(cna_ids, patient_ids)

# Cols 1 and 2 are gene names
assay <- round(as.matrix(data_Rnaseq[, -c(1, 2)])) 
rownames(assay) <- data_Rnaseq[, 1]

# Build metadata.
metadata <- matrix(0, ncol(assay), 2)

# Creating metadata with HER2 status
erbb2_row <- which(data_cna$Hugo_Symbol == "ERBB2")  # Finding the line ERBB2
for (i in 1:dim(assay)[2]) {
  pat_barcode <- substr(colnames(assay)[i], 1, 12)
  # Check if the id is in common_ids
  if (pat_barcode %in% common_ids) {
    idx <- match(pat_barcode, cna_ids)
    if (!is.na(idx)) {
      metadata[i, 1] <- pat_barcode
      metadata[i, 2] <- ifelse(as.numeric(data_cna[erbb2_row, idx + 2]) > 0, 1, 0)
    }
  }else {
    metadata[i, 1] <- pat_barcode
    metadata[i, 2] <- NA
  } 
}
colnames(metadata) <- c("Patient_ID", "HER2_Status") # Setting column names
# Changing metadata types
metadata <- as.data.frame(metadata) 
metadata$HER2_Status <- as.factor(metadata$HER2_Status)

#anyNA(metadata) # Checking for missing values
#sum(is.na(metadata)) # Calculating amount of missing values

metadata <- metadata[complete.cases(metadata), ] # Remove lines containing NA


# Chunk 3
#_______________________________________________
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# Install DeSeq2
#if (!require("DESeq2", quietly = TRUE))
#    BiocManager::install("DESeq2")


# Chunk 4
#_______________________________________________
assay[assay<0] = 0

# Findind identifiers that are not in metadata
ids_to_remove <- setdiff(substr(colnames(assay), 1, 12), metadata[, 1])
# Finding the indexes of the columns corresponding to these identifiers
cols_to_remove <- which(substr(colnames(assay), 1, 12) %in% ids_to_remove)
# Removing columns from assay
filtered_assay <- assay[, -cols_to_remove]

# Check: order must be identical to colnames(assay)
# all(colnames(filtered_assay) == metadata$`Patient ID`)

# Filtering out genes with too many missing values
smallestGroupSize = 3
keep = rowSums(filtered_assay >= 10) >= smallestGroupSize
filtered_assay = filtered_assay[keep,]


# Chunk 5
#_______________________________________________
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = filtered_assay,
  colData = metadata,
  design = ~ HER2_Status
)
dds <- DESeq(dds)  # Normalization and analysis


# Chunk 6
#_______________________________________________
resultsNames(dds)
res = results(dds)
print(res)
#Top 10 most differentially expressed genes with the smallest Fold Change 
res[order(res$log2FoldChange)[1:10], ]
# Top 10 most differentially expressed genes with the largest Fold Change
res[order(res$log2FoldChange, decreasing = TRUE)[1:10], ]


# Chunk 7
#_______________________________________________
#if (!require("EnhancedVolcano", quietly = TRUE))
#  BiocManager::install("EnhancedVolcano")
#if (!require("ggplot2", quietly = TRUE))
#  install.packages("ggplot2")
#if (!require("ggrepel", quietly = TRUE))
#  install.packages("ggrepel")


# Chunk 8
#_______________________________________________
library(EnhancedVolcano)
pdf("VolcanoPlotHER2.pdf", width = 10, height = 5.8)
# Building a volcano plot
EnhancedVolcano(
  res,  # DESeq2 results
  lab = rownames(res),  # Genes names
  x = 'log2FoldChange',  # Axis X: log2FoldChange
  y = 'padj',  # Axis Y: Adjusted p-values 
  pCutoff = 0.05,  # p-value significance threshold
  FCcutoff = 1.5,  # Change threshold ($fold change = 2^{1.5} \\approx 2.8$)
  title = 'Volcano Plot: HER2 Amplified vs Not Amplified',
  subtitle = 'Differentially Expressed Genes',
  caption = 'log2FoldChange vs Adjusted p-value',
  pointSize = 2,  # Dot size
  labSize = 3,  # Text size for genes
  xlim = c(-5, 5),  # X-axis limits
  ylim = c(0, -log10(min(res$padj, na.rm = TRUE))),  # Y-axis limits
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'), # Dot colors
  colAlpha = 0.6,  # Transparency of dots
  legendPosition = 'top',  # Legend Position
  legendLabSize = 10
)
dev.off()


# Chunk 9
#_______________________________________________
vsd = vst(dds)
pdf("PCA.pdf", width = 8, height = 6)
plotPCA(vsd, intgroup=c("HER2_Status"))
dev.off()


# Chunk 10
#_______________________________________________
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("fgsea")
#if (!requireNamespace("clusterProfiler", quietly = TRUE))
#  BiocManager::install("clusterProfiler")
#if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
#  BiocManager::install("org.Hs.eg.db")
#if (!requireNamespace("enrichplot", quietly = TRUE))
#  install.packages("enrichplot")


# Chunk 11
#_______________________________________________
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# get subset of differentially expressed genes.
res_sig = res[res$padj<0.05,] # Filtering genes with adjusted p-value $(padj) < 0.05$
head(res_sig)

# separate into over and under expressed using log2foldchange
DE_over = rownames(res_sig[res_sig$log2FoldChange>0,]) # List of overexpressed genes 
DE_under = rownames(res_sig[res_sig$log2FoldChange<0,]) 

# Enrichment analysis for overexpressed genes
go_results_over = enrichGO(
  gene          = DE_over,      # Vector of overexpressed genes
  OrgDb         = org.Hs.eg.db, # Human Gene Annotation Database
  keyType       = "SYMBOL",     # Gene identifier type (gene symbols)
  ont           = "BP",         # Ontology: Biological Process
  pAdjustMethod = "BH",         # p-value adjustment method (Benjamini-Hochberg)
  pvalueCutoff  = 0.05,         # p-value threshold 
  qvalueCutoff  = 0.05          # q-value threshold 
)
# print and plot results
print(head(go_results_over))
dotplot(go_results_over, showCategory=10) + ggtitle("Gene Ontology Enrichment over Expressed")

# Enrichment analysis for underexpressed genes
go_results_under = enrichGO(
  gene          = DE_under,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
# print and plot results
print(head(go_results_under))
dotplot(go_results_under, showCategory=10) + ggtitle("Gene Ontology Enrichment Under Expressed")


# Chunk 12
#_______________________________________________
#if (!requireNamespace("pathview", quietly = TRUE))
#  BiocManager::install("pathview")
#if (!requireNamespace("ReactomePA", quietly = TRUE))
#  BiocManager::install("ReactomePA")


# Chunk 13
#_______________________________________________
library(ReactomePA)
library(pathview)

# Convert gene symbols (SYMBOL) to their ENTREZ identifiers to use them in Reactome and KEGG analysis
gene_entrez_over <- bitr(
  DE_over,
  fromType = "SYMBOL",    # Type of source data identifier (gene symbols)
  toType   = "ENTREZID",  # Target ID type (ENTREZ ID)
  OrgDb    = org.Hs.eg.db
)
gene_entrez_under <- bitr(
  DE_under,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

#KEGG
kegg_results_over =  enrichKEGG(
  gene          = gene_entrez_over[,2],  # Column with ENTREZ ID from transformed data
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
print(head(kegg_results_over))
dotplot(kegg_results_over, showCategory=10) + ggtitle("Kegg Pathway Enrichment Over Expressed")

kegg_results_under =  enrichKEGG(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
print(head(kegg_results_under))
dotplot(kegg_results_under, showCategory=10) + ggtitle("Kegg Pathway Enrichment Under Expressed")


# Chunk 14
#_______________________________________________
# Reactome
pdf("ReactomeOver.pdf", width = 10, height = 5.8)
reactome_results_over =  enrichPathway(
  gene          = gene_entrez_over[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
)
print(head(reactome_results_over))
dotplot(reactome_results_over, showCategory=10) + ggtitle("Reactome Pathway Enrichment Over Expressed")
dev.off()

pdf("ReactomeUnder.pdf", width = 10, height = 5.8)
reactome_results_under =  enrichPathway(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
)
print(head(reactome_results_under))
dotplot(reactome_results_under, showCategory=10) + ggtitle("Reactome Pathway Enrichment Under Expressed")
dev.off()


# Chunk 15
#_______________________________________________
# GO
# go_results_over_pw = pairwise_termsim(go_results_over)
# treeplot(go_results_over_pw) + ggtitle("GO Enrichment Over Expressed")
# go_results_under_pw = pairwise_termsim(go_results_under)
# treeplot(go_results_under_pw) + ggtitle("Go Enrichment Under Expressed")
# KEGG
# kegg_results_over_pw = pairwise_termsim(kegg_results_over)
# treeplot(kegg_results_over_pw) + ggtitle("KEGG Enrichment Over Expressed")
# kegg_results_under_pw = pairwise_termsim(kegg_results_under)
# treeplot(kegg_results_under_pw) + ggtitle("KEGG Enrichment Under Expressed")
# REACTOME
# reactome_results_over_pw = pairwise_termsim(reactome_results_over)
# treeplot(reactome_results_over_pw) + ggtitle("Reactome Enrichment Over Expressed")
# reactome_results_under_pw = pairwise_termsim(reactome_results_under)
# treeplot(reactome_results_under_pw) + ggtitle("Reactome Enrichment Under Expressed")

reactome_results_over_pw <- pairwise_termsim(reactome_results_over)
treeplot(reactome_results_over_pw, main = "Reactome Enrichment Over Expressed")
reactome_results_under_pw <- pairwise_termsim(reactome_results_under)
treeplot(reactome_results_under_pw, main = "Reactome Enrichment Under Expressed")


# Chunk 16
#_______________________________________________
# install packages for nicer heatmap than R's base one. 
# if (!requireNamespace("pheatmap", quietly = TRUE))
#  install.packages("pheatmap")
library(pheatmap)


# Chunk 17
#_______________________________________________
# Subset the dataset on differentially expressed gene 
top_DE = order(res$padj)
# To simplify visualization, the top 20 genes with the most significant expression are selected
vsd_DE = assay(vsd)[top_DE[1:20],]
# Creation annotations for columns (samples) based on metadata.
annotation_col = data.frame(HER2_Status = as.matrix(metadata[,2]))
rownames(annotation_col) = colnames(vsd) # # Set row names to match column names in vsd.
# Change colors for HER2_Status annotation
annotation_colors <- list(
  HER2_Status = c("0" = "springgreen",  
                  "1" = "blue") 
)
pdf("Heatmap.pdf", width = 10, height = 6)
pheatmap(
  vsd_DE,                   # Top 20 gene expression data
  cluster_rows = TRUE,      # Clustering of strings (genes)
  cluster_cols = TRUE,      # Clustering of columns (samples)
  scale = 'row',            # Row-wise scaling (highlights relative differences)
  show_colnames = FALSE,    # Hiding column names (samples)
  show_rownames = TRUE,     # Displaying names of genes
  annotation_col = annotation_col,           # Adding annotations to columns
  annotation_colors = annotation_colors)     # Applying specified colors to annotations
dev.off()


# Chunk 18
#_______________________________________________
# if (!requireNamespace("survival", quietly = TRUE)) {
#  install.packages("survival")}
# if (!requireNamespace("glmnet", quietly = TRUE)) {
#  install.packages("glmnet")}
library(survival)
library(glmnet)
library(survminer)


# Chunk 19
#_______________________________________________
de_genes <- rownames(res[res$padj < 0.05, ]) # DE gene selection
x <- t(assay(vsd)[de_genes, ])  # VST values of DE genes, transposed (patients x genes)

# Conversion of identifier formats and preparation of survival data
surv <- data_patient[, c('X.Patient.Identifier', 'Overall.Survival.Status', 'Overall.Survival..Months.')]
surv$X.Patient.Identifier <- substr(gsub("-", "\\.", surv$X.Patient.Identifier), 1, 12)
rownames(x) <- substr(rownames(x), 1, 12)                  # Transforming patient ids
surv <- surv[surv$X.Patient.Identifier %in% rownames(x), ] # Leave only matching ids

# Transforming survival data
time <- as.numeric(surv$Overall.Survival..Months.)  # Survival Time
status <- as.numeric(ifelse(surv$Overall.Survival.Status == "1:DECEASED", 1, 0))  # 1 = dead, 0 = alive

# Remove lines with invalid survival time
set.seed(123) # For reproducibility
valid_indices <- time > 0 # Check that time > 0
time <- time[valid_indices]
status <- status[valid_indices]
x <- x[valid_indices, ]  # Removing rows with invalid data

y <- Surv(time, status) # Creating a Survival Object

set.seed(123) # For reproducibility
fit <- glmnet(x, y, family = "cox")  # Building a Lasso-regularized Cox model
# plot(fit) 

# Selecting the optimal value of λ through cross-validation
set.seed(123) # For reproducibility
cv_fit <- cv.glmnet(x, y, family = "cox", type.measure = "deviance")
# plot(cv_fit)
optimal_lambda <- cv_fit$lambda.min # Optimal value of λ
# cat("Optimal Lambda:", optimal_lambda, "\n")

predicted_risk <- predict(fit, newx = x, s = optimal_lambda, type = "link") # Risk prediction
# Grouping of patients by risk (median risk value)
risk_groups <- ifelse(predicted_risk > median(predicted_risk), "High Risk", "Low Risk")

# Creating and visualizing survival curves
surv_fit <- survfit(y ~ risk_groups)
pdf("SA.pdf", width = 10, height = 6)
ggsurvplot(surv_fit, data = data.frame(time, status, risk_groups), pval = TRUE)
dev.off()


# Chunk 20
#_______________________________________________
de_genes <- rownames(res[res$padj < 0.05, ])  
x <- t(assay(vsd)[de_genes, ])
surv <- data_patient[, c('X.Patient.Identifier', 'Overall.Survival.Status', 'Overall.Survival..Months.')]
surv$X.Patient.Identifier <- substr(gsub("-", "\\.", surv$X.Patient.Identifier), 1, 12)
rownames(x) <- substr(rownames(x), 1, 12) 
surv <- surv[surv$X.Patient.Identifier %in% rownames(x), ]

# Оставляем только идентификаторы в metadata, которые есть в rownames(x)
metadata <- metadata[metadata$Patient_ID %in% rownames(x), ]
# Упорядочиваем metadata так, чтобы идентификаторы совпадали с rownames(x)
metadata <- metadata[match(rownames(x), metadata$Patient_ID), ]


time <- as.numeric(surv$Overall.Survival..Months.)  
status <- as.numeric(ifelse(surv$Overall.Survival.Status == "1:DECEASED", 1, 0))  
set.seed(123) # For reproducibility
valid_indices <- time > 0  
time <- time[valid_indices]
status <- status[valid_indices]
x <- x[valid_indices, ] 

metadata <- metadata[valid_indices, ] # Filter metadata to only include rows that correspond to patients with a valid survival time 

y <- Surv(time, status) 
set.seed(123) # For reproducibility
fit <- glmnet(x, y, family = "cox")
set.seed(123) # For reproducibility
cv_fit <- cv.glmnet(x, y, family = "cox", type.measure = "deviance")
# plot(cv_fit)
optimal_lambda <- cv_fit$lambda.min
# cat("Optimal Lambda:", optimal_lambda, "\n")
predicted_risk <- predict(fit, newx = x, s = optimal_lambda, type = "link")

# Creating risk groups based on patients' HER2 status
risk_groups <- ifelse(metadata$HER2_Status == 1, "HER2 Amplified", "HER2 Not Amplified")

surv_fit <- survfit(y ~ risk_groups)
ggsurvplot(surv_fit, data = data.frame(time, status, risk_groups), pval = TRUE)