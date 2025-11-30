#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: mismatch_single.R <trait_name> <path/to/gwas.QC.gz>")
}

trait <- args[1]
gwas_path <- args[2]

# Paths (relative for GitHub)
bim_path <- "data/1000G_QC/target_1000G_finalQC.bim"
qc_path  <- "data/1000G_QC/target_1000G.QC.snplist"
out_dir  <- "data/mismatch/"

# Load target BIM
bim <- read.table(bim_path, stringsAsFactors = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")

# Load QC SNP list
qc <- read.table(qc_path, stringsAsFactors = FALSE)

# Load GWAS summary
gwas <- read.table(gzfile(gwas_path), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Uppercase alleles
gwas$A1 <- toupper(gwas$ALT)
gwas$A2 <- toupper(gwas$REF)
bim$B.A1 <- toupper(bim$B.A1)
bim$B.A2 <- toupper(bim$B.A2)

# Merge base + target
info <- merge(bim, gwas, by = c("SNP", "CHR", "BP"))

# Keep QC-passed SNPs
info <- info[info$SNP %in% qc$V1, ]

# Complement function
complement <- function(x) {
  switch(x, "A"="T", "T"="A", "C"="G", "G"="C", NA)
}

# Exact matches
info.match <- subset(info, A1 == B.A1 & A2 == B.A2)

# Complement matches
info$C.A1 <- sapply(info$B.A1, complement)
info$C.A2 <- sapply(info$B.A2, complement)
info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)

comp.idx <- bim$SNP %in% info.complement$SNP
bim[comp.idx, "B.A1"] <- sapply(bim[comp.idx, "B.A1"], complement)
bim[comp.idx, "B.A2"] <- sapply(bim[comp.idx, "B.A2"], complement)

# Recode
info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
rec.idx <- bim$SNP %in% info.recode$SNP
tmp <- bim[rec.idx, "B.A1"]
bim[rec.idx, "B.A1"] <- bim[rec.idx, "B.A2"]
bim[rec.idx, "B.A2"] <- tmp

# Complement + recode
info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
crec.idx <- bim$SNP %in% info.crecode$SNP
tmp <- bim[crec.idx, "B.A1"]
bim[crec.idx, "B.A1"] <- sapply(bim[crec.idx, "B.A2"], complement)
bim[crec.idx, "B.A2"] <- sapply(tmp, complement)

# Output updated A1 allele file
write.table(
  bim[, c("SNP", "B.A1")],
  file = paste0(out_dir, trait, ".a1"),
  col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t"
)

# Output mismatched SNPs
mismatch <- bim$SNP[
  !(bim$SNP %in% info.match$SNP |
      bim$SNP %in% info.complement$SNP |
      bim$SNP %in% info.recode$SNP |
      bim$SNP %in% info.crecode$SNP)
]

write.table(
  mismatch,
  file = paste0(out_dir, trait, ".mismatch"),
  col.names = FALSE, row.names = FALSE, quote = FALSE
)
