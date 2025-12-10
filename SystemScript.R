# ----------------
# import packages
# ----------------

library(dplyr)
library(recount3)
library(recount)
library(edgeR)
library(patchwork)
library(data.table)
library(Matrix)

# ---------------
# load data and TPM
# ---------------

data <- read.csv("GSE133206_raw_counts_GRCh38.p13_NCBI.tsv", sep="\t")

assays(data)$TPM <- recount::getTPM(data)


# ---------------
# remove pseudo, short,  mito, rRNA
# ---------------

rr <- rowRanges(data)
mito <- seqnames(rr) == "chrM"
pseudo <- grepl("pseudogene", rr$gene_biotype, ignore.case = TRUE)
rRNA <- grepl("rRNA", rr$gbkey, ignore.case = TRUE)
short <- rr$bp_length < 200

to_remove <- mito | pseudo | rRNA | short
data_filtered <- data[!to_remove,]
counts_data_selected <- assays(data_selected)$counts

# ---------------
# remove 0 and do logcpm
# ---------------

y <- DGEList(counts = data)

keep.exprs <- filterByExpr(y, group = group)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]

logcpm_before <- cpm(y, log = TRUE)
y <- calcNormFactors(y, method = "TMM")
y

# ---------------
# extract results
# ---------------

group <- as.factor(c("Healthy", "Cancer"))

design <- model.matrix(~0+group, data = y$samples)

qlf <- glmQLFTest(fit, contrast = c(-1,1))

results <- topTags(qlf, n = 10000000, adjust.method = "BH", 
                     sort.by = "PValue", p.value = 0.01)

write.csv(results, file = "results.csv")

