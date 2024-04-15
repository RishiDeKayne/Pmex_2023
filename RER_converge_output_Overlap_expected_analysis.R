#analysis of overlap between RER shifted genes and known candidat loci and divergently expressed genes

##analyse RERconverge output
setwd("/Users/rishide-kayne/Dropbox/RishiMac2/Kelley_Lab/Kerry_project/rer_converge/")

#got list of genes from kara
KARA_DF <- read.csv("Pmex_target_genes.tsv", sep = "\t")

KS <- subset(KARA_DF, KARA_DF$Gene.Set == "Sulfide Detox/Response")
KN <- subset(KARA_DF, KARA_DF$Gene.Set == "NuOX")

KARA_DF <- subset(KARA_DF, nchar(KARA_DF$Human_ID)>0)

KARA_DF_Sulfide_Detox <- subset(KARA_DF, KARA_DF$Gene.Set == "Sulfide Detox/Response")
KARA_DF_oxphos <- subset(KARA_DF, KARA_DF$Gene.Set == "NuOX")

nrow(KARA_DF_Sulfide_Detox)
nrow(KARA_DF_oxphos)

Kara_test <- rbind(KARA_DF_Sulfide_Detox, KARA_DF_oxphos)
kara_test_unique <- unique(Kara_test$Human_ID)

#load gene data TABLE 1
RERgenes <- read.csv("../RER_Manuscript_GBE_2024_02_20/TableS2.txt", sep = "")

RERgenes_with_anno <- subset(RERgenes, RERgenes$Annotation1 != "no_orthology")

RERgenes_with_anno$overlap <- "no"

for(i in 1:nrow(RERgenes_with_anno)){
  gene_name <- as.character(RERgenes_with_anno$Annotation1)[i]
  if(gene_name %in% KARA_DF_Sulfide_Detox$Human_ID){
    RERgenes_with_anno$overlap[i] <- "yes"
  }
  if(gene_name %in% KARA_DF_oxphos$Human_ID){
    RERgenes_with_anno$overlap[i] <- "yes"
  }
}

RERgenes_with_anno_overlap_yes <- subset(RERgenes_with_anno, RERgenes_with_anno$overlap == "yes")
RERgenes_with_anno_overlap_no <- subset(RERgenes_with_anno, RERgenes_with_anno$overlap == "no")

nrow(RERgenes_with_anno)
nrow(RERgenes_with_anno_overlap_yes)

prop_overlap_all <- nrow(RERgenes_with_anno_overlap_yes)/nrow(RERgenes_with_anno)
prop_overlap_all

sig_RER <- subset(RERgenes_with_anno, RERgenes_with_anno$Significant == "yes")
nrow(sig_RER)

sigRERgenes_with_anno_overlap_yes <- subset(sig_RER, sig_RER$overlap == "yes")

nrow(sigRERgenes_with_anno_overlap_yes)

prop_overlap_RERsig <- nrow(sigRERgenes_with_anno_overlap_yes)/nrow(sig_RER)
prop_overlap_RERsig

#fisher input
sig_no_overlap <- nrow(sig_RER)-nrow(sigRERgenes_with_anno_overlap_yes)
non_sig_overlap <- nrow(RERgenes_with_anno_overlap_yes)-nrow(sigRERgenes_with_anno_overlap_yes)
non_sig_no_overlap <- nrow(RERgenes_with_anno)-non_sig_overlap-nrow(sig_RER)
sig_overlap <- nrow(sigRERgenes_with_anno_overlap_yes)

sum(non_sig_no_overlap, sig_no_overlap, non_sig_overlap, sig_overlap)

# dat <- data.frame(
#   "No_overlap_with_candidates" = c(non_sig_no_overlap, sig_no_overlap),
#   "Overlap_with_candidates" = c(non_sig_overlap, sig_overlap),
#   row.names = c("nonsig_RER", "sig_RER"),
#   stringsAsFactors = FALSE
# )
# colnames(dat) <- c("No Overlap", "Overlap")

dat <- data.frame(
  "Candidate" = c(sig_overlap, sig_no_overlap),
  "Not Candidate" = c(non_sig_overlap, non_sig_no_overlap),
  row.names = c("Candidate", "Not Candidate"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("rateShift", "NoRateShift")

dat

mosaicplot(dat,
           main = "Mosaic plot",
           color = TRUE
)

fisher.test(dat)

#load gene data TABLE 2
RERgenes <- read.csv("../RER_Manuscript_GBE_2024_02_20/TableS3.tsv", sep = "")

RERgenes_with_anno <- subset(RERgenes, RERgenes$Annotation1 != "no_orthology")

RERgenes_with_anno$overlap <- "no"

for(i in 1:nrow(RERgenes_with_anno)){
  gene_name <- as.character(RERgenes_with_anno$Annotation1)[i]
  if(gene_name %in% KARA_DF_Sulfide_Detox$Human_ID){
    RERgenes_with_anno$overlap[i] <- "yes"
  }
  if(gene_name %in% KARA_DF_oxphos$Human_ID){
    RERgenes_with_anno$overlap[i] <- "yes"
  }
}

RERgenes_with_anno_overlap_yes <- subset(RERgenes_with_anno, RERgenes_with_anno$overlap == "yes")
nrow(RERgenes_with_anno)
nrow(RERgenes_with_anno_overlap_yes)

prop_overlap_all <- nrow(RERgenes_with_anno_overlap_yes)/nrow(RERgenes_with_anno)
prop_overlap_all

sig_RER <- subset(RERgenes_with_anno, RERgenes_with_anno$Significant == "yes")
nrow(sig_RER)

sigRERgenes_with_anno_overlap_yes <- subset(sig_RER, sig_RER$overlap == "yes")

nrow(sigRERgenes_with_anno_overlap_yes)

prop_overlap_RERsig <- nrow(sigRERgenes_with_anno_overlap_yes)/nrow(sig_RER)
prop_overlap_RERsig

#fisher input
sig_no_overlap <- nrow(sig_RER)-nrow(sigRERgenes_with_anno_overlap_yes)
non_sig_overlap <- nrow(RERgenes_with_anno_overlap_yes)-nrow(sigRERgenes_with_anno_overlap_yes)
non_sig_no_overlap <- nrow(RERgenes_with_anno)-non_sig_overlap-nrow(sig_RER)
sig_overlap <- nrow(sigRERgenes_with_anno_overlap_yes)

sum(non_sig_no_overlap, sig_no_overlap, non_sig_overlap, sig_overlap)

# dat <- data.frame(
#   "No_overlap_with_candidates" = c(non_sig_no_overlap, sig_no_overlap),
#   "Overlap_with_candidates" = c(non_sig_overlap, sig_overlap),
#   row.names = c("nonsig_RER", "sig_RER"),
#   stringsAsFactors = FALSE
# )
# colnames(dat) <- c("No Overlap", "Overlap")

dat <- data.frame(
  "Candidate" = c(sig_overlap, sig_no_overlap),
  "Not Candidate" = c(non_sig_overlap, non_sig_no_overlap),
  row.names = c("Candidate", "Not Candidate"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("rateShift", "NoRateShift")

dat

mosaicplot(dat,
           main = "Mosaic plot",
           color = TRUE
)

fisher.test(dat)


