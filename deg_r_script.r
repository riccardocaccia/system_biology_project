# -------------------------
# Load packages
# -------------------------
library(dplyr)
library(edgeR)
library(data.table)

# -------------------------
# Load raw count matrix
# -------------------------
# Assuming first column = gene_id, other columns = samples
raw <- fread("GSE133206_raw_counts_GRCh38.p13_NCBI.tsv")

# Set gene IDs as rownames
genes <- raw[[1]]
raw_mat <- as.data.frame(raw[,-1])
rownames(raw_mat) <- genes

# -------------------------
# remove unwanted genes
# -----------------------
# If you do NOT have annotation, skip removal or provide annotation separately.

# mito   <- grepl("^MT-", genes)
# pseudo <- grepl("pseudogene", raw$gene_biotype, ignore.case=TRUE)
# rRNA   <- grepl("rRNA", raw$gene_type, ignore.case=TRUE)
# short  <- raw$length < 200
# keep   <- !(mito | pseudo | rRNA | short)
# raw_mat <- raw_mat[keep,]

# -------------------------
# Create group factor
# -------------------------
group <- factor(c(rep("GSM3902322", N1),  #normal
                  rep("GSM3902321",  N2)))  #cancer

# -------------------------
# Build DGEList
# -------------------------
y <- DGEList(counts = raw_mat, group = group)

# -------------------------
# Filter low-expressed genes
# -------------------------
keep <- filterByExpr(y, group = group)
y <- y[keep,, keep.lib.sizes = FALSE]

# -------------------------
# Normalize TMM
# -------------------------
y <- calcNormFactors(y)

# -------------------------
# Design matrix
# -------------------------
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# -------------------------
# Estimate dispersions
# -------------------------
y <- estimateDisp(y, design)

# -------------------------
# Fit GLM
# -------------------------
fit <- glmQLFit(y, design)

# cancer vs normal
contrast <- makeContrasts(Cancer_vs_normaly = GSM3902321 - GSM3902322,
                          levels = design)

# -------------------------
# Test
# -------------------------
qlf <- glmQLFTest(fit, contrast = contrast)

# -------------------------
# extract results
# -------------------------
results <- topTags(qlf, n = Inf, adjust.method = "BH")

# save the table
write.csv(results$table, "DEG_results.csv", row.names = TRUE)

# -------------------------
# output summary
# -------------------------
summary(decideTests(qlf))
