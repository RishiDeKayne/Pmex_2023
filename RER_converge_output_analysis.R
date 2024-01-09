#analyse RERconverge output
setwd("/Users/rishide-kayne/Dropbox/RishiMac2/Kelley_Lab/Kerry_project/rer_converge/")

#for genes I have a gene number and a corresponding rna-XM number
#for promoters I have multiple GENEIDs which correspond to a specific rna-XM number

#############ANNOTATION LOAD
##make lookup table with all genes in gff
gff_raw <- read.delim("gff_annotation/GCF_003331165.1_Xiphophorus_hellerii-4.1_genomic.gff", header=F, comment.char="#")

#Get genes
gff.genes <- gff_raw[gff_raw[,3]=="gene",]
#row9_genes <- as.data.frame(gff.genes$V9)
row9_genes_split <- data.frame(do.call('rbind', strsplit(as.character(gff.genes$V9), ":", fixed = TRUE)))
row9_genes_splitsplit <- as.data.frame(do.call('rbind', strsplit(as.character(row9_genes_split$X2), "Name", fixed = T)))
row9_genes_splitsplit_ID <- gsub(";", "", row9_genes_splitsplit$V1)
gff.genes$ID_number <- row9_genes_splitsplit_ID

gff.genes$promoter_window_value <- "empty"
for(i in 1:nrow(gff.genes)){
  if(as.character(gff.genes$V7[i]) == "+"){
    gff.genes$promoter_window_value[i] <- paste(gff.genes$V1[i], (gff.genes$V4[i]-500), gff.genes$V4[i], sep = "_")
  }
  if(as.character(gff.genes$V7[i]) == "-"){
    gff.genes$promoter_window_value[i] <- paste(gff.genes$V1[i], gff.genes$V5[i], (gff.genes$V5[i]+500), sep = "_")
  }
}
gff.genes$midpoint <- as.integer(((gff.genes$V5+gff.genes$V4)/2))

promoter_input <- read.csv("promoters_500bp_xmac.merged.bed", header = F, sep = "\t")


gff.exon <- gff_raw[gff_raw[,3]=="exon",]
#row9_genes <- as.data.frame(gff.genes$V9)
row9_genes_split <- data.frame(do.call('rbind', strsplit(as.character(gff.exon$V9), "GeneID:", fixed = TRUE)))
row9_genes_splitsplit <- as.data.frame(do.call('rbind', strsplit(as.character(row9_genes_split$X2), ",Genbank", fixed = T)))
gff.exon$ID_number <- row9_genes_splitsplit$V1
row9_genes_split <- data.frame(do.call('rbind', strsplit(as.character(gff.exon$V9), "transcript_id=", fixed = TRUE)))
gff.exon$rna_name <- row9_genes_split$X2
gff.exon$rna_name_with_beg <- gsub("X", "rna-X", gff.exon$rna_name)

#then load orthology
blastp_output <- read.csv("blast_output/xhel_human_blastp_results.sorted.top1.tab", sep = "\t", header = F)
blastp_output$locusNUMB <- gsub("id-LOC", "", blastp_output$V1)
#need to split column 2 based on | 
blastp_output_split <- data.frame(do.call('rbind', strsplit(as.character(blastp_output$V2), "|", fixed = TRUE)))
blastp_output$human_gene1 <- blastp_output_split$X2
blastp_output$human_gene2 <- blastp_output_split$X3

##############################################################
#                     ALL
##############################################################
par(mfrow=c(2,2))
#with new analysis 
output_stats <- read.csv("allvalues_RER.csv")

#get rid of gene bit
rownames(output_stats) <- gsub("gene_", "", rownames(output_stats))

#plot histogram
hist(output_stats$Rho, breaks = 18, xlim = c(-0.6, 0.6), ylim = c(0, 3000), main = '', xlab = "relative evolutionary rate (Rho)")

#work out upper and lower 5% bounds and plot ablines
a <- quantile(x = output_stats$Rho, probs = c(0.05,0.95))
abline(v=as.numeric(a[1]), col = "red", lty = 2, lwd = 2)
abline(v=as.numeric(a[2]), col = "red", lty = 2, lwd = 2)

#plot p value histogram and add 0.05 line
hist(output_stats$P, breaks = 18, main = '', xlab = "p-value", ylim = c(0, 2000))
abline(v=as.numeric(0.05), col = "red", lty = 2, lwd = 2)

#make subset of significant observations
sig <- subset(output_stats, output_stats$P < 0.05)

#plot relationship of these with rho
#plot(sqrt(sig$Rho^2), sig$P, ylab = "p-value", xlab = "|Rho|")

#add gene numbers to sig outliers
sig$gene_numb <- rownames(sig)

#load gene name info with ncbi gene names
names <- read.csv("gene_names_and_numbers.txt", sep = " ", head = FALSE)

sig$gene_name <- c()

#go in loop and add ncbi gene name to each row
for(i in 1:length(sig$Rho)){
  numby_line <- sig[i,]
  numby <- numby_line$gene_numb
  line <- subset(names, as.character(names$V1) == as.character(numby))
  sig$gene_name[i] <- line$V2
}

#do the same for full dataset too
output_stats$gene_numb <- rownames(output_stats)
output_stats$gene_name <- c()
for(i in 1:length(output_stats$Rho)){
  numby_line <- output_stats[i,]
  numby <- numby_line$gene_numb
  line <- subset(names, as.character(names$V1) == as.character(numby))
  output_stats$gene_name[i] <- line$V2
}

#get rid of rna bit of name
output_stats$gene_name <- gsub("rna-","", output_stats$gene_name)

#get gene annotation info
product <- read.csv("gene_name_and_product_from_gff_longest.txt", sep = "\t", head = FALSE)

#append this to the df too
output_stats$product <- c()
for(i in 1:length(output_stats$Rho)){
  output_line <- output_stats[i,]
  output_gene_name <- output_line$gene_name
  product_name <- subset(product, as.character(product$V1) == as.character(output_gene_name))
  output_stats$product[i] <- product_name$V2
}

#now load permutation info
perm_output_stats_rho <- read.csv("permuted_all_RHO_RER.csv")
perm_output_stats_rho$V1 <- gsub("gene_", "", perm_output_stats_rho$V1)

perm_output_stats_p <- read.csv("permuted_all_P_RER.csv")
perm_output_stats_p$V1 <- gsub("gene_", "", perm_output_stats_p$V1)

#count the number of test-statistics as or more extreme than our initial test statistic
#divide that number by the total number of test-statistics we calculated

#take the first column of gene names
gene_names_list <- as.data.frame(perm_output_stats_p[,1])
#then get all the permuted combos
p_numbs_df <- perm_output_stats_p[,2:128]
#then get just our last proper sulfidic test values which will be used as threshold
threshold <- perm_output_stats_p[,129]

#do test
permute_P_list <- c()
for(j in 1:nrow(p_numbs_df)){
  threshold_value <- threshold[j]
  vec <- as.numeric(p_numbs_df[j,])
  more_sig_count <- sum(vec <= threshold_value)
  permute_P <- (more_sig_count+1)/128
  permute_P_list[j] <- permute_P
}

#count how many are empirically significant = 1689
sum(permute_P_list <= 0.05)

#make new df with rho, p and empirical p
new_all_df <- cbind(gene_names_list, threshold, permute_P_list, perm_output_stats_rho[,129])
colnames(new_all_df) <- c("gene", "P", "permuted_P", "Rho")

#and now find outliers that had a significant P and are empirically significant
new_all_df$symbol <- 1
for(i in 1:nrow(new_all_df)){
  if(new_all_df$P[i] < 0.05 && new_all_df$permuted_P[i] < 0.05){
    new_all_df$symbol[i] <- 16
  }
}

#and make fancy plot
outliers_red <- subset(new_all_df, new_all_df$symbol == 16)
outliers_red <- subset(outliers_red, outliers_red$Rho > 0)
outliers_blue <- subset(new_all_df, new_all_df$symbol == 16)
outliers_blue <- subset(outliers_blue, outliers_blue$Rho < 0)

plot(new_all_df$Rho, -log10(new_all_df$P), yaxt = 'n', ylab = 'Original P-value', xlab = "Rho", xlim = c(-0.8, 0.8), ylim = c(0,4), main = "all genes", col = "grey")
points(outliers_red$Rho, -log10(outliers_red$P), col = "darkred", pch = 16)
points(outliers_blue$Rho, -log10(outliers_blue$P), col = "darkblue", pch = 16)
axis(side = 2, at = c(0.0, 1, 2, 3), labels = c("1", "0.1", "0.01", "0.001"), tick = TRUE)

plot(new_all_df$Rho, -log10(new_all_df$permuted_P), yaxt = 'n', ylab = 'Empirical P-value', xlab = "Rho", xlim = c(-0.8, 0.8), ylim = c(0,3), main = "all genes", col = "grey")
points(outliers_red$Rho, -log10(outliers_red$permuted_P), col = "darkred", pch = 16)
points(outliers_blue$Rho, -log10(outliers_blue$permuted_P), col = "darkblue", pch = 16)
axis(side = 2, at = c(0.0, 1, 2, 3), labels = c("1", "0.1", "0.01", "0.001"), tick = TRUE)

print(paste("total number of genes analysed: ", as.character(nrow(output_stats))))
print(paste("in all genes there are ", as.character(nrow(outliers_blue)), " genes evolving slower than expected"))
print(paste("in all genes there are ", as.character(nrow(outliers_red)), " genes evolving faster than expected"))

#standard gene set enrichment test:
#prepare GO file


ALL_GENES <- new_all_df

##############################################################
#                     Promoters
##############################################################
#with new analysis 
output_stats <- read.csv("allvalues_promoters_RER.csv")

#get rid of gene bit
rownames(output_stats) <- gsub("gene_", "", rownames(output_stats))

#plot histogram
hist(output_stats$Rho, breaks = 18, xlim = c(-0.6, 0.6), ylim = c(0, 2500), main = '', xlab = "relative evolutionary rate (Rho)")

#work out upper and lower 5% bounds and plot ablines
a <- quantile(x = output_stats$Rho, probs = c(0.05,0.95))
abline(v=as.numeric(a[1]), col = "red", lty = 2, lwd = 2)
abline(v=as.numeric(a[2]), col = "red", lty = 2, lwd = 2)

#plot p value histogram and add 0.05 line
hist(output_stats$P, breaks = 18, main = '', xlab = "p-value", ylim = c(0, 2000))
abline(v=as.numeric(0.05), col = "red", lty = 2, lwd = 2)

#make subset of significant observations
sig <- subset(output_stats, output_stats$P < 0.05)

#plot relationship of these with rho
#plot(sqrt(sig$Rho^2), sig$P, ylab = "p-value", xlab = "|Rho|")

#add gene numbers to sig outliers
sig$gene_numb <- rownames(sig)

#now load permutation info
perm_output_stats_rho <- read.csv("permuted_promoters_RHO_RER.csv")

perm_output_stats_p <- read.csv("permuted_promoters_P_RER.csv")

#count the number of test-statistics as or more extreme than our initial test statistic
#divide that number by the total number of test-statistics we calculated

#take the first column of gene names
gene_names_list <- as.data.frame(perm_output_stats_p[,1])
#then get all the permuted combos
p_numbs_df <- perm_output_stats_p[,2:128]
#then get just our last proper sulfidic test values which will be used as threshold
threshold <- perm_output_stats_p[,129]

#do test
permute_P_list <- c()
for(j in 1:nrow(p_numbs_df)){
  threshold_value <- threshold[j]
  vec <- as.numeric(p_numbs_df[j,])
  more_sig_count <- sum(vec <= threshold_value)
  permute_P <- (more_sig_count+1)/128
  permute_P_list[j] <- permute_P
}

#count how many are empirically significant = 1810
sum(permute_P_list <= 0.05)

#make new df with rho, p and empirical p
new_all_df <- cbind(gene_names_list, threshold, permute_P_list, perm_output_stats_rho[,129])
colnames(new_all_df) <- c("gene", "P", "permuted_P", "Rho")

#and now find outliers that had a significant P and are empirically significant
new_all_df$symbol <- 1
for(i in 1:nrow(new_all_df)){
  if(new_all_df$P[i] < 0.05 && new_all_df$permuted_P[i] < 0.05){
    new_all_df$symbol[i] <- 16
  }
}

#and make fancy plot
outliers_red <- subset(new_all_df, new_all_df$symbol == 16)
outliers_red <- subset(outliers_red, outliers_red$Rho > 0)
outliers_blue <- subset(new_all_df, new_all_df$symbol == 16)
outliers_blue <- subset(outliers_blue, outliers_blue$Rho < 0)

plot(new_all_df$Rho, -log10(new_all_df$P), yaxt = 'n', ylab = 'Original P-value', xlab = "Rho", xlim = c(-0.8, 0.8), ylim = c(0,4), main = "promoters", col = "grey")
points(outliers_red$Rho, -log10(outliers_red$P), col = "darkred", pch = 16)
points(outliers_blue$Rho, -log10(outliers_blue$P), col = "darkblue", pch = 16)
axis(side = 2, at = c(0.0, 1, 2, 3), labels = c("1", "0.1", "0.01", "0.001"), tick = TRUE)

plot(new_all_df$Rho, -log10(new_all_df$permuted_P), yaxt = 'n', ylab = 'Empirical P-value', xlab = "Rho", xlim = c(-0.8, 0.8), ylim = c(0,3), main = "promoters", col = "grey")
points(outliers_red$Rho, -log10(outliers_red$permuted_P), col = "darkred", pch = 16)
points(outliers_blue$Rho, -log10(outliers_blue$permuted_P), col = "darkblue", pch = 16)
axis(side = 2, at = c(0.0, 1, 2, 3), labels = c("1", "0.1", "0.01", "0.001"), tick = TRUE)

print(paste("total number of promoter regions analysed: ", as.character(nrow(output_stats))))
print(paste("in promoters there are ", as.character(nrow(outliers_blue)), " genes evolving slower than expected"))
print(paste("in promoters there are ", as.character(nrow(outliers_red)), " genes evolving faster than expected"))


ALL_PROMOTERS <- new_all_df

prom_window_split <- data.frame(do.call('rbind', strsplit(as.character(ALL_PROMOTERS$gene), "_", fixed = TRUE)))
ALL_PROMOTERS$prom_value_summary <- prom_window_split$X4
ALL_PROMOTERS$prom_value_chr <- paste(prom_window_split$X1, prom_window_split$X2, sep = "_")

ALL_GENES$rna_name <- "empty"
ALL_GENES$GENE_ID_numb <- "empty"
ALL_GENES$promoter_window <- "empty"
ALL_GENES$promoter_rho <- "empty"
ALL_GENES$promoter_perm_p <- "empty"
ALL_GENES$promoter_symbol <- "empty"
ALL_GENES$human_gene_name1 <- "empty"
ALL_GENES$human_gene_name2 <- "empty"
ALL_GENES$E_val <- "empty"

for(i in 1:(nrow(ALL_GENES))){
  gene_chr <- as.character(ALL_GENES$gene[i])
  name_row <- subset(names, as.character(names$V1) == gene_chr)
  if(nrow(name_row) > 0){
    ALL_GENES$rna_name[i] <- name_row$V2
  }
}

for(i in 1:(nrow(ALL_GENES))){
  rna_name_chr <- as.character(ALL_GENES$rna_name[i])
  exons_row <- subset(gff.exon, as.character(gff.exon$rna_name_with_beg) == rna_name_chr)
  if(nrow(exons_row) > 0){
    ALL_GENES$GENE_ID_numb[i] <- exons_row$ID_number[1]
    gene_gff_sub <- subset(gff.genes, gff.genes$ID_number == exons_row$ID_number[1])
    ALL_GENES$promoter_window[i] <- gene_gff_sub$promoter_window_value
  }
}

for(i in 1:(nrow(ALL_GENES))){
  prom_window_ID <- as.character(ALL_GENES$promoter_window[i])
  rna_id <- as.character(ALL_GENES$rna_name[i])
  promoter_df <- subset(ALL_PROMOTERS, ALL_PROMOTERS$gene == prom_window_ID)
  if(nrow(promoter_df) > 0){
    ALL_GENES$promoter_rho[i] <- promoter_df$Rho
    ALL_GENES$promoter_perm_p[i] <- promoter_df$permuted_P
    ALL_GENES$promoter_symbol[i] <- promoter_df$symbol
  }
  blastp_row <- subset(blastp_output, as.character(blastp_output$locusNUMB) == rna_id)
  if(nrow(blastp_row) > 0){
    ALL_GENES$human_gene_name1[i] <- blastp_row$human_gene1
    ALL_GENES$human_gene_name2[i] <- blastp_row$human_gene2
    ALL_GENES$E_val[i] <- blastp_row$V11
  }
}

ALL_GENES$promoter_rho <- as.numeric(ALL_GENES$promoter_rho)
length(which(ALL_GENES$promoter_symbol == 16))

sig_gene_sig_prom <- subset(ALL_GENES, ALL_GENES$symbol == 16 & ALL_GENES$promoter_symbol == 16)
non_gene_sig_prom <- subset(ALL_GENES, ALL_GENES$symbol == 1 & ALL_GENES$promoter_symbol == 16)
sig_gene_non_prom <- subset(ALL_GENES, ALL_GENES$symbol == 16 & ALL_GENES$promoter_symbol == 1)
non_gene_non_prom <- subset(ALL_GENES, ALL_GENES$symbol == 1 & ALL_GENES$promoter_symbol == 1)

nrow(sig_gene_sig_prom)
nrow(non_gene_sig_prom)
nrow(sig_gene_non_prom)
nrow(non_gene_non_prom)

length(which(non_gene_sig_prom$promoter_rho > 0))
length(which(non_gene_sig_prom$promoter_rho < 0))

length(which(sig_gene_non_prom$Rho > 0))
length(which(sig_gene_non_prom$Rho < 0))

total <- sum(nrow(sig_gene_sig_prom), nrow(non_gene_sig_prom), nrow(sig_gene_non_prom), nrow(non_gene_non_prom))
total

ALL_GENES$E_val <- as.numeric(ALL_GENES$E_val)
ALL_GENES_CONFIDENT_E <- subset(ALL_GENES, ALL_GENES$E_val < 0.05)

webgestalt_2023_11_26 <- as.data.frame(cbind(ALL_GENES_CONFIDENT_E$human_gene_name1, ALL_GENES_CONFIDENT_E$Rho))
webgestalt_2023_11_26 <- subset(webgestalt_2023_11_26, webgestalt_2023_11_26$V1 != "empty")
webgestalt_2023_11_26$V2 <- as.numeric(webgestalt_2023_11_26$V2)
webgestalt_2023_11_26 <- webgestalt_2023_11_26[order(webgestalt_2023_11_26$V2),]
write.table(webgestalt_2023_11_26,file="webgestalt_2023_11_26.txt",sep="\t",col.names = F,row.names = F, quote = F)

OUTLIER_WEBGESTALT <- ALL_GENES_CONFIDENT_E
OUTLIER_WEBGESTALT_BACKGROUND <- as.data.frame(OUTLIER_WEBGESTALT$human_gene_name1)
write.table(OUTLIER_WEBGESTALT_BACKGROUND,file="OUTLIER_WEBGESTALT_BACKGROUND.txt",sep="\t",col.names = F,row.names = F, quote = F)

OUTLIER_WEBGESTALT_FAST <- as.data.frame(subset(OUTLIER_WEBGESTALT, OUTLIER_WEBGESTALT$permuted_P < 0.05 & OUTLIER_WEBGESTALT$Rho > 0 & OUTLIER_WEBGESTALT$P < 0.05))
OUTLIER_WEBGESTALT_FAST_LIST <- as.data.frame(OUTLIER_WEBGESTALT_FAST$human_gene_name1)
write.table(OUTLIER_WEBGESTALT_FAST_LIST,file="OUTLIER_WEBGESTALT_FAST.txt",sep="\t",col.names = F,row.names = F, quote = F)

OUTLIER_WEBGESTALT_SLOW <- as.data.frame(subset(OUTLIER_WEBGESTALT, OUTLIER_WEBGESTALT$permuted_P < 0.05 & OUTLIER_WEBGESTALT$Rho < 0 & OUTLIER_WEBGESTALT$P < 0.05))
OUTLIER_WEBGESTALT_SLOW_LIST <- as.data.frame(OUTLIER_WEBGESTALT_SLOW$human_gene_name1)
write.table(OUTLIER_WEBGESTALT_SLOW_LIST,file="OUTLIER_WEBGESTALT_SLOW.txt",sep="\t",col.names = F,row.names = F, quote = F)


ALL_PROMOTERS$human_gene_name1 <- "empty"
ALL_PROMOTERS$human_gene_name2 <- "empty"
ALL_PROMOTERS$E_val <- "empty"

for(i in 1:nrow(ALL_PROMOTERS)){
  gff_row <- subset(gff.genes, as.character(gff.genes$promoter_window_value) == as.character(ALL_PROMOTERS$gene[i]))
  if(nrow(gff_row) > 0){
    gene_number <- gff_row$ID_number
    gff_exon_row <- subset(gff.exon, as.character(gff.exon$ID_number) == as.character(gene_number))
    if(nrow(gff_exon_row) > 0){
      RNA_name <- gff_exon_row[1,]$rna_name_with_beg
      blast_result <- subset(blastp_output, blastp_output$locusNUMB == RNA_name)
      if(nrow(blast_result) > 0){
        ALL_PROMOTERS$human_gene_name1[i] <- blast_result$human_gene1
        ALL_PROMOTERS$human_gene_name2[i] <- blast_result$human_gene2
        ALL_PROMOTERS$E_val[i] <- blast_result$V11
      }
    }
  }
}

ALL_PROMOTERS$E_val <- as.numeric(ALL_PROMOTERS$E_val)
ALL_PROMOTERS_CONFIDENT_E <- subset(ALL_PROMOTERS, ALL_PROMOTERS$E_val < 0.05)

webgestalt_promoter_2023_11_26 <- as.data.frame(cbind(ALL_PROMOTERS_CONFIDENT_E$human_gene_name1, ALL_PROMOTERS_CONFIDENT_E$Rho))
webgestalt_promoter_2023_11_26 <- subset(webgestalt_promoter_2023_11_26, webgestalt_promoter_2023_11_26$V1 != "empty")
webgestalt_promoter_2023_11_26$V2 <- as.numeric(webgestalt_promoter_2023_11_26$V2)
webgestalt_promoter_2023_11_26 <- webgestalt_promoter_2023_11_26[order(webgestalt_promoter_2023_11_26$V2, decreasing = TRUE),]
write.table(webgestalt_promoter_2023_11_26,file="webgestalt_promoter_2023_11_26.txt",sep="\t",col.names = F,row.names = F, quote = F)

OUTLIER_WEBGESTALT_P <- ALL_PROMOTERS_CONFIDENT_E
OUTLIER_WEBGESTALT_P_BACKGROUND <- as.data.frame(OUTLIER_WEBGESTALT_P$human_gene_name1)
write.table(OUTLIER_WEBGESTALT_P_BACKGROUND,file="OUTLIER_WEBGESTALT_P_BACKGROUND.txt",sep="\t",col.names = F,row.names = F, quote = F)

OUTLIER_WEBGESTALT_P_FAST <- as.data.frame(subset(OUTLIER_WEBGESTALT_P, OUTLIER_WEBGESTALT_P$permuted_P < 0.05 & OUTLIER_WEBGESTALT_P$Rho > 0 & OUTLIER_WEBGESTALT_P$P < 0.05))
OUTLIER_WEBGESTALT_P_FAST_LIST <- as.data.frame(OUTLIER_WEBGESTALT_P_FAST$human_gene_name1)
write.table(OUTLIER_WEBGESTALT_P_FAST_LIST,file="OUTLIER_WEBGESTALT_P_FAST_LIST.txt",sep="\t",col.names = F,row.names = F, quote = F)

OUTLIER_WEBGESTALT_P_SLOW <- as.data.frame(subset(OUTLIER_WEBGESTALT_P, OUTLIER_WEBGESTALT_P$permuted_P < 0.05 & OUTLIER_WEBGESTALT_P$Rho < 0 & OUTLIER_WEBGESTALT_P$P < 0.05))
OUTLIER_WEBGESTALT_P_SLOW_LIST <- as.data.frame(OUTLIER_WEBGESTALT_P_SLOW$human_gene_name1)
write.table(OUTLIER_WEBGESTALT_P_SLOW_LIST,file="OUTLIER_WEBGESTALT_P_SLOW_LIST.txt",sep="\t",col.names = F,row.names = F, quote = F)

##############################################################
#                     Enhancer
##############################################################
#with new analysis 
output_stats <- read.csv("allvalues_enhancers_RER.csv")

#get rid of gene bit
rownames(output_stats) <- gsub("gene_", "", rownames(output_stats))

#plot histogram
hist(output_stats$Rho, breaks = 18, xlim = c(-0.6, 0.6), ylim = c(0, 1500), main = '', xlab = "relative evolutionary rate (Rho)")

#work out upper and lower 5% bounds and plot ablines
a <- quantile(x = output_stats$Rho, probs = c(0.05,0.95))
abline(v=as.numeric(a[1]), col = "red", lty = 2, lwd = 2)
abline(v=as.numeric(a[2]), col = "red", lty = 2, lwd = 2)

#plot p value histogram and add 0.05 line
hist(output_stats$P, breaks = 18, main = '', xlab = "p-value", ylim = c(0, 1000))
abline(v=as.numeric(0.05), col = "red", lty = 2, lwd = 2)

#make subset of significant observations
sig <- subset(output_stats, output_stats$P < 0.05)

#plot relationship of these with rho
#plot(sqrt(sig$Rho^2), sig$P, ylab = "p-value", xlab = "|Rho|")

#add gene numbers to sig outliers
sig$gene_numb <- rownames(sig)

#now load permutation info
perm_output_stats_rho <- read.csv("permuted_enhancers_RHO_RER.csv")

perm_output_stats_p <- read.csv("permuted_enhancers_P_RER.csv")

#count the number of test-statistics as or more extreme than our initial test statistic
#divide that number by the total number of test-statistics we calculated

#take the first column of gene names
gene_names_list <- as.data.frame(perm_output_stats_p[,1])
#then get all the permuted combos
p_numbs_df <- perm_output_stats_p[,2:128]
#then get just our last proper sulfidic test values which will be used as threshold
threshold <- perm_output_stats_p[,129]

#do test
permute_P_list <- c()
for(j in 1:nrow(p_numbs_df)){
  threshold_value <- threshold[j]
  vec <- as.numeric(p_numbs_df[j,])
  more_sig_count <- sum(vec <= threshold_value)
  permute_P <- (more_sig_count+1)/128
  permute_P_list[j] <- permute_P
}

#count how many are empirically significant = 686
sum(permute_P_list <= 0.05)

#make new df with rho, p and empirical p
new_all_df <- cbind(gene_names_list, threshold, permute_P_list, perm_output_stats_rho[,129])
colnames(new_all_df) <- c("gene", "P", "permuted_P", "Rho")

#and now find outliers that had a significant P and are empirically significant
new_all_df$symbol <- 1
for(i in 1:nrow(new_all_df)){
  if(new_all_df$P[i] < 0.05 && new_all_df$permuted_P[i] < 0.05){
    new_all_df$symbol[i] <- 16
  }
}

#and make fancy plot
outliers_red <- subset(new_all_df, new_all_df$symbol == 16)
outliers_red <- subset(outliers_red, outliers_red$Rho > 0)
outliers_blue <- subset(new_all_df, new_all_df$symbol == 16)
outliers_blue <- subset(outliers_blue, outliers_blue$Rho < 0)

plot(new_all_df$Rho, -log10(new_all_df$P), yaxt = 'n', ylab = 'Original P-value', xlab = "Rho", xlim = c(-0.8, 0.8), ylim = c(0,4), main = "enhancers", col = "grey")
points(outliers_red$Rho, -log10(outliers_red$P), col = "darkred", pch = 16)
points(outliers_blue$Rho, -log10(outliers_blue$P), col = "darkblue", pch = 16)
axis(side = 2, at = c(0.0, 1, 2, 3), labels = c("1", "0.1", "0.01", "0.001"), tick = TRUE)

plot(new_all_df$Rho, -log10(new_all_df$permuted_P), yaxt = 'n', ylab = 'Empirical P-value', xlab = "Rho", xlim = c(-0.8, 0.8), ylim = c(0,3), main = "enhancers", col = "grey")
points(outliers_red$Rho, -log10(outliers_red$permuted_P), col = "darkred", pch = 16)
points(outliers_blue$Rho, -log10(outliers_blue$permuted_P), col = "darkblue", pch = 16)
axis(side = 2, at = c(0.0, 1, 2, 3), labels = c("1", "0.1", "0.01", "0.001"), tick = TRUE)

print(paste("total number of enhancer regions analysed: ", as.character(nrow(output_stats))))
print(paste("in enhancers there are ", as.character(nrow(outliers_blue)), " genes evolving slower than expected"))
print(paste("in enhancers there are ", as.character(nrow(outliers_red)), " genes evolving faster than expected"))

ALL_ENHANCERS <- new_all_df

enhancer_input <- read.csv("../regulatory/putativeEnhancerPeaks_xmac_10.27.23.merged.bed", head=FALSE, sep = "\t")
enhancer_metadata <- read.csv("./PutativeEnhancerGeneAssociations_11.06.23.tsv", sep = "\t")

library(stringr)
enhancer_input$count <- str_count(enhancer_input$V4,"\\|")

enhancer_input$window <- paste(enhancer_input$V1, enhancer_input$V2, enhancer_input$V3, sep = "_")
  
row9_genes_split_NAME <- data.frame(do.call('rbind', strsplit(as.character(gff.genes$V9), "Dbxref", fixed = TRUE)))
row9_genes_split_NAME_COL <- gsub("ID=gene-", "", row9_genes_split_NAME$X1)
row9_genes_split_NAME_COL <- gsub(";", "", row9_genes_split_NAME_COL)

gff.genes$NAME <- row9_genes_split_NAME_COL

enhancer_input_one <- subset(enhancer_input, enhancer_input$count < 1)

ALL_ENHANCERS$matches <- "empty"
ALL_ENHANCERS$locus_numb <- "empty"
ALL_ENHANCERS$rna_log2FC <- "empty"
ALL_ENHANCERS$rna_corrPval <- "empty"
ALL_ENHANCERS$csrna_log2FC <- "empty"
ALL_ENHANCERS$srna_corrPval <- "empty"
ALL_ENHANCERS$rna_direction <- "empty"
ALL_ENHANCERS$csrna_direction <- "empty"
ALL_ENHANCERS$relationship <- "empty"
ALL_ENHANCERS$human_gene_name1 <- "empty"
ALL_ENHANCERS$human_gene_name2 <- "empty"
ALL_ENHANCERS$E_val <- "empty"


# for(i in 1:nrow(ALL_ENHANCERS)){
#   window <- ALL_ENHANCERS$gene[1]
#   enhancer_input_sub <- subset(enhancer_input, enhancer_input$window == window)
#   if(nrow(enhancer_input_sub)>0){
#     county <- enhancer_input_sub$count
#     ALL_ENHANCERS$matches[i] <- county
#     if(county == 0){
#       region <- enhancer_input_sub$V4
#       enhancer_metadata_sub <- subset(enhancer_metadata, enhancer_metadata$id_peak == region)
#       if(nrow(enhancer_metadata_sub)>0){
#         ALL_ENHANCERS$locus_numb[i] <- enhancer_metadata_sub$symbol
#         ALL_ENHANCERS$rna_log2FC[i] <- enhancer_metadata_sub$rna_log2FC
#         ALL_ENHANCERS$rna_corrPval[i] <- enhancer_metadata_sub$rna_corrPval
#         ALL_ENHANCERS$csrna_log2FC[i] <- enhancer_metadata_sub$csrna_log2FC
#         ALL_ENHANCERS$srna_corrPval[i] <- enhancer_metadata_sub$srna_corrPval
#         ALL_ENHANCERS$rna_direction[i] <- enhancer_metadata_sub$rna_direction
#         ALL_ENHANCERS$csrna_direction[i] <- enhancer_metadata_sub$csrna_direction
#         ALL_ENHANCERS$relationship[i] <- enhancer_metadata_sub$relationship
#         
#         gff_sub <- subset(gff.genes, gff.genes$NAME)
#       }
#     }
#   }
# }

#got list of genes from kara
KARA_DF <- read.csv("Pmex_target_genes.tsv", sep = "\t")

KARA_DF$gene_count_sig <- "empty"
KARA_DF$gene_count_nonsig <- "empty"

for(i in 1:nrow(KARA_DF)){
  gene_name <- KARA_DF$Human_ID[i]
  gene_df <- subset(ALL_GENES_CONFIDENT_E, as.character(ALL_GENES_CONFIDENT_E$human_gene_name1) == as.character(gene_name))
  gene_df_sig <- subset(gene_df, gene_df$symbol == 16)
  if(nrow(gene_df_sig)>0){
    newstring <- as.character(gene_df_sig$Rho)
    test <- paste(newstring, collapse = ";")
    KARA_DF$gene_count_sig[i] <- test
  }
  gene_df_non_sig <- subset(gene_df, gene_df$symbol == 1)
  if(nrow(gene_df_non_sig)>0){
    newstring <- as.character(gene_df_non_sig$Rho)
    test <- paste(newstring, collapse = ";")
    KARA_DF$gene_count_nonsig[i] <- test
  }
}

KARA_DF$prom_count_sig <- "empty"
KARA_DF$prom_count_nonsig <- "empty"


for(i in 1:nrow(KARA_DF)){
  gene_name <- KARA_DF$Human_ID[i]
  gene_df <- subset(ALL_PROMOTERS_CONFIDENT_E, as.character(ALL_PROMOTERS_CONFIDENT_E$human_gene_name1) == as.character(gene_name))
  gene_df_sig <- subset(gene_df, gene_df$symbol == 16)
  if(nrow(gene_df_sig)>0){
    newstring <- as.character(gene_df_sig$Rho)
    test <- paste(newstring, collapse = ";")
    KARA_DF$prom_count_sig[i] <- test
  }
  gene_df_nonsig <- subset(gene_df, gene_df$symbol == 1)
  if(nrow(gene_df_nonsig)>0){
    newstring <- as.character(gene_df_nonsig$Rho)
    test <- paste(newstring, collapse = ";")
    KARA_DF$prom_count_nonsig[i] <- test
  }
}

KARA_DF_Sulfide_Detox_Response_sig <- subset(KARA_DF, KARA_DF$Gene.Set == "Sulfide Detox/Response")
KARA_DF_Sulfide_Detox_Response_sig_genes <- subset(KARA_DF_Sulfide_Detox_Response_sig, KARA_DF_Sulfide_Detox_Response_sig$gene_count_sig != "empty")
KARA_DF_Sulfide_Detox_Response_sig_promoters <- subset(KARA_DF_Sulfide_Detox_Response_sig, KARA_DF_Sulfide_Detox_Response_sig$prom_count_sig != "empty")

unique(KARA_DF_Sulfide_Detox_Response_sig_genes$Human_ID)
unique(KARA_DF_Sulfide_Detox_Response_sig_genes$Human_gene_name)

unique(KARA_DF_Sulfide_Detox_Response_sig_promoters$Human_ID)
unique(KARA_DF_Sulfide_Detox_Response_sig_promoters$Human_gene_name)


KARA_DF_Sulfide_Detox <- subset(KARA_DF, KARA_DF$Gene.Set == "Sulfide Detox/Response")
KARA_DF_oxphos <- subset(KARA_DF, KARA_DF$Gene.Set == "NuOX")

nrow(KARA_DF_Sulfide_Detox)
nrow(KARA_DF_oxphos)

##############################################################
#                     SUPP FIGURES
##############################################################

#FIGURE S1
par(mfrow=c(1,3))
par(mar=c(4.5,4,2,2))

set <- ALL_GENES
hist(set$Rho, breaks = 18, xlim = c(-0.6, 0.6), ylim = c(0, 3000), main = 'Genes', xlab = "Relative evolutionary rate (Rho)")

#work out upper and lower 5% bounds and plot ablines
a <- quantile(x = set$Rho, probs = c(0.05,0.95))
abline(v=as.numeric(a[1]), col = "red", lty = 2, lwd = 2)
abline(v=as.numeric(a[2]), col = "red", lty = 2, lwd = 2)

set <- ALL_PROMOTERS
hist(set$Rho, breaks = 18, xlim = c(-0.6, 0.6), ylim = c(0, 3000), main = 'Promoters', xlab = "Relative evolutionary rate (Rho)")

#work out upper and lower 5% bounds and plot ablines
a <- quantile(x = set$Rho, probs = c(0.05,0.95))
abline(v=as.numeric(a[1]), col = "red", lty = 2, lwd = 2)
abline(v=as.numeric(a[2]), col = "red", lty = 2, lwd = 2)

set <- ALL_ENHANCERS
hist(set$Rho, breaks = 18, xlim = c(-0.6, 0.6), ylim = c(0, 3000), main = 'Enhancers', xlab = "Relative evolutionary rate (Rho)")

#work out upper and lower 5% bounds and plot ablines
a <- quantile(x = set$Rho, probs = c(0.05,0.95))
abline(v=as.numeric(a[1]), col = "red", lty = 2, lwd = 2)
abline(v=as.numeric(a[2]), col = "red", lty = 2, lwd = 2)


#FIGURE S2
par(mfrow=c(1,3))
par(mar=c(4.5,4,2,2))

set <- ALL_GENES
outliers_red <- subset(set, set$symbol == 16)
outliers_red <- subset(outliers_red, outliers_red$Rho > 0)
outliers_blue <- subset(set, set$symbol == 16)
outliers_blue <- subset(outliers_blue, outliers_blue$Rho < 0)

plot(set$Rho, -log10(set$permuted_P), yaxt = 'n', ylab = 'Permuted P-value', xlab = "Rho", xlim = c(-0.8, 0.8), ylim = c(0,3), col = "grey", main = "Genes")
points(outliers_red$Rho, -log10(outliers_red$permuted_P), col = "darkred", pch = 16)
points(outliers_blue$Rho, -log10(outliers_blue$permuted_P), col = "darkblue", pch = 16)
axis(side = 2, at = c(0.0, 1, 2, 3), labels = c("1", "0.1", "0.01", "0.001"), tick = TRUE)

set <- ALL_PROMOTERS
outliers_red <- subset(set, set$symbol == 16)
outliers_red <- subset(outliers_red, outliers_red$Rho > 0)
outliers_blue <- subset(set, set$symbol == 16)
outliers_blue <- subset(outliers_blue, outliers_blue$Rho < 0)

plot(set$Rho, -log10(set$permuted_P), yaxt = 'n', ylab = 'Permuted P-value', xlab = "Rho", xlim = c(-0.8, 0.8), ylim = c(0,3), col = "grey", main = "Promoters")
points(outliers_red$Rho, -log10(outliers_red$permuted_P), col = "darkred", pch = 16)
points(outliers_blue$Rho, -log10(outliers_blue$permuted_P), col = "darkblue", pch = 16)
axis(side = 2, at = c(0.0, 1, 2, 3), labels = c("1", "0.1", "0.01", "0.001"), tick = TRUE)

set <- ALL_ENHANCERS
outliers_red <- subset(set, set$symbol == 16)
outliers_red <- subset(outliers_red, outliers_red$Rho > 0)
outliers_blue <- subset(set, set$symbol == 16)
outliers_blue <- subset(outliers_blue, outliers_blue$Rho < 0)

plot(set$Rho, -log10(set$permuted_P), yaxt = 'n', ylab = 'Permuted P-value', xlab = "Rho", xlim = c(-0.8, 0.8), ylim = c(0,3), col = "grey", main = "Enhancers")
points(outliers_red$Rho, -log10(outliers_red$permuted_P), col = "darkred", pch = 16)
points(outliers_blue$Rho, -log10(outliers_blue$permuted_P), col = "darkblue", pch = 16)
axis(side = 2, at = c(0.0, 1, 2, 3), labels = c("1", "0.1", "0.01", "0.001"), tick = TRUE)

#Figure S3
par(mfrow=c(1,1))

ALL_GENES_WITHPROMOTERS <- subset(ALL_GENES, ALL_GENES$promoter_perm_p != "empty")
plot(ALL_GENES_WITHPROMOTERS$Rho, ALL_GENES_WITHPROMOTERS$promoter_rho, ylab = "Promoter relative evolutionary rate (Rho)", xlab = "Gene relative evolutionary rate (Rho)")

corr <- cor.test(x=ALL_GENES_WITHPROMOTERS$Rho, y=ALL_GENES_WITHPROMOTERS$promoter_rho, method = 'spearman')
corr

##############################################################
#                     SUPP FILES
##############################################################

#Table S2 - Gene RER results
#gene, rna name, gene ID number, Rho, p, permuted p, sig, annotation1, annotation2, Eval
S2_DF <- as.data.frame(cbind(ALL_GENES$gene, ALL_GENES$rna_name, ALL_GENES$GENE_ID_numb, ALL_GENES$Rho, ALL_GENES$P, ALL_GENES$permuted_P, ALL_GENES$symbol, ALL_GENES$human_gene_name1, ALL_GENES$human_gene_name2, ALL_GENES$E_val))
colnames(S2_DF) <- c("Gene_Number", "RNA_Name", "GENE_ID", "Rho", "P", "Permuted_P", "Significant", "Annotation1", "Annotation2", "BLASTp_Eval")
S2_DF$Significant <- gsub("16", "yes", S2_DF$Significant)
S2_DF$Significant <- gsub("1", "no", S2_DF$Significant)
S2_DF$Annotation1 <- gsub("empty", "no_orthology", S2_DF$Annotation1)
S2_DF$Annotation2 <- gsub("empty", "no_orthology", S2_DF$Annotation2)

write.table(S2_DF,file="Table_S2.tsv",sep="\t",col.names = T,row.names = F, quote = F)

#Table S3 - Promoter RER results
ALL_PROMOTERS$GENE_ID <- "no_gene_id"
for(i in 1:nrow(ALL_PROMOTERS)){
  gff_row <- subset(gff.genes, as.character(gff.genes$promoter_window_value) == as.character(ALL_PROMOTERS$gene[i]))
  if(nrow(gff_row) > 0){
    ALL_PROMOTERS$GENE_ID[i] <- gff_row$ID_number
  }
}
S3_DF <- as.data.frame(cbind(ALL_PROMOTERS$gene, ALL_PROMOTERS$Rho, ALL_PROMOTERS$P, ALL_PROMOTERS$permuted_P, ALL_PROMOTERS$symbol, ALL_PROMOTERS$GENE_ID, ALL_PROMOTERS$human_gene_name1, ALL_PROMOTERS$human_gene_name2, ALL_PROMOTERS$E_val))
colnames(S3_DF) <- c("Window", "Rho", "P", "Permuted_P", "Significant", "Gene_ID", "Annotation1", "Annotation2", "BLASTp_Eval")
S3_DF$Significant <- gsub("16", "yes", S3_DF$Significant)
S3_DF$Significant <- gsub("1", "no", S3_DF$Significant)
S3_DF$Annotation1 <- gsub("empty", "no_orthology", S3_DF$Annotation1)
S3_DF$Annotation2 <- gsub("empty", "no_orthology", S3_DF$Annotation2)

write.table(S3_DF,file="Table_S3.tsv",sep="\t",col.names = T,row.names = F, quote = F)

#Table S4 - Enhancer RER results
S4_DF <- as.data.frame(cbind(ALL_ENHANCERS$gene, ALL_ENHANCERS$Rho, ALL_ENHANCERS$P, ALL_ENHANCERS$permuted_P, ALL_ENHANCERS$symbol))
colnames(S4_DF) <- c("Window", "Rho", "P", "Permuted_P", "Significant")
S4_DF$Significant <- gsub("16", "yes", S4_DF$Significant)
S4_DF$Significant <- gsub("1", "no", S4_DF$Significant)

write.table(S4_DF,file="Table_S4.tsv",sep="\t",col.names = T,row.names = F, quote = F)

write.table(blastp_output,file="blastp_output.tsv",sep="\t",col.names = T,row.names = F, quote = F)

#Table S5 - genes and promoters
S5_DF <- as.data.frame(cbind(ALL_GENES_WITHPROMOTERS$gene, ALL_GENES_WITHPROMOTERS$rna_name, ALL_GENES_WITHPROMOTERS$GENE_ID_numb, ALL_GENES_WITHPROMOTERS$Rho, ALL_GENES_WITHPROMOTERS$P, ALL_GENES_WITHPROMOTERS$permuted_P, ALL_GENES_WITHPROMOTERS$promoter_window, ALL_GENES_WITHPROMOTERS$Rho, ALL_GENES_WITHPROMOTERS$promoter_perm_p, ALL_GENES_WITHPROMOTERS$human_gene_name1, ALL_GENES_WITHPROMOTERS$human_gene_name2, ALL_GENES_WITHPROMOTERS$E_val))
colnames(S5_DF) <- c("Gene_Number", "RNA_Name", "GENE_ID", "Rho", "P", "Permuted_P", "Promoter_window", "Promoter_Rho", "Promoter_Permuted_P", "Annotation1", "Annotation2", "BLASTp_Eval")
write.table(S5_DF,file="Table_S5.tsv",sep="\t",col.names = T,row.names = F, quote = F)

#Table S6 - overlap with previous work (done in excel)
#Table S7 - GSEA/GO input gene lists (done in excel)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#                 OLD CODE SNIPPETS
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

library(ape)
test_tree <- read.tree(text = "(MX15-627,(MX15-734,(((MX15-706,MX15-698)1:3.317559746479158,((FL15-049,FL15-025)1:2.185863612463782,(MX15-650,MX15-739)1:0.42883480509971506)1:2.5950159412430422)1:0.2717517145091955,((CR15-031,CR15-039)1:3.2696745939783423,((FL15-043,FL15-032)1:2.208420613733077,(vs,esperanza)1:1.7730915877528486)1:1.9351398460512297)1:2.900536294688129)1:3.0523799368871):0.0);")
plot(test_tree)

#############################################################################################
library(tidyverse)
# library(clusterProfiler)
# BiocManager::install("org.Dr.eg.db")
library(org.Dr.eg.db)

# Read in and reformat orthogroup table -----------------------------------

orthos <- read_tsv('~/Downloads/RER.protein__v__Danio_rerio.GRCz11.pep.all.tsv') 

orthos.long <- orthos %>% 
  separate_longer_delim(cols = RER.protein,delim = ', ') %>%  
  mutate(first_danio_id = str_split_fixed(Danio_rerio.GRCz11.pep.all,', ',2)[,1]) %>% 
  dplyr::select(xhel_id=RER.protein,first_danio_id,Orthogroup) %>% 
  mutate(first_danio_id = str_split_fixed(first_danio_id,'[.]',2)[,1])

orthos.long


# Read in accelerated genes, get orthologs --------------------------------
genes.acc <- read_csv('~/Downloads/GO_test_sig_genes_faster.csv') %>% 
  dplyr::select(xhel_id=1) %>% 
  left_join(orthos.long) %>% 
  mutate(xhel_id = str_split_fixed(xhel_id,'[.]',2)[,1])

genes.acc  

acc.entrez <- AnnotationDbi::select(org.Dr.eg.db, 
                                    keys = genes.acc$first_danio_id,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "ENSEMBL")

genes.acc.fullInfo <- genes.acc %>% left_join(acc.entrez,by=c('first_danio_id'='ENSEMBL')) %>% 
  dplyr::select(xhel_id,danio_ensembl=2,danio_entrez=4,danio_symbol=5,Orthogroup)


################################################################################
################################################################################

output_stats <- read.csv("allvalues_enhancers_RER.csv")
output_stats <- na.omit(output_stats)
hist(output_stats$Rho, breaks = 18, xlim = c(-0.6, 0.6), ylim = c(0, 1500), main = '', xlab = "relative evolutionary rate (Rho)")
a <- quantile(x = output_stats$Rho, probs = c(0.05,0.95))
abline(v=as.numeric(a[1]), col = "red", lty = 2, lwd = 2)
abline(v=as.numeric(a[2]), col = "red", lty = 2, lwd = 2)

hist(output_stats$P, breaks = 18, main = '', xlab = "p-value")
abline(v=as.numeric(0.05), col = "red", lty = 2, lwd = 2)

sig <- subset(output_stats, output_stats$P < 0.05)

plot(sqrt(sig$Rho^2), sig$P, ylab = "p-value", xlab = "|Rho|")


####
for(i in 1:(nrow(ALL_GENES))){
  gene_id_chr <- as.character(ALL_GENES$GENE_ID_numb[i])
  rna_id <- as.character(ALL_GENES$rna_name[i])
  #rna_id <- "rna-XM_032549199.1"
  genes_row <- subset(gff.genes, as.character(gff.genes$ID_number) == gene_id_chr)
  if(nrow(genes_row) > 0){
    prom_end <- genes_row$promoter_value
    prom_chrom <- genes_row$V1
    promoter_df <- subset(ALL_PROMOTERS, ALL_PROMOTERS$prom_value_summary == prom_end)
    if(nrow(promoter_df) > 0 && nrow(promoter_df) < 2){
      ALL_GENES$promoter_window[i] <- promoter_df$gene
      ALL_GENES$promoter_rho[i] <- promoter_df$Rho
      ALL_GENES$promoter_perm_p[i] <- promoter_df$permuted_P
      ALL_GENES$promoter_symbol[i] <- promoter_df$symbol
    }
    if(nrow(promoter_df) > 1){
      promoter_df2 <- subset(promoter_df, promoter_df$prom_value_chr == prom_chrom)
      ALL_GENES$promoter_window[i] <- promoter_df2$gene
      ALL_GENES$promoter_rho[i] <- promoter_df2$Rho
      ALL_GENES$promoter_perm_p[i] <- promoter_df2$permuted_P
      ALL_GENES$promoter_symbol[i] <- promoter_df2$symbol
    }
  }
  blastp_row <- subset(blastp_output, as.character(blastp_output$locusNUMB) == rna_id)
  if(nrow(blastp_row) > 0){
    ALL_GENES$human_gene_name1[i] <- blastp_row$human_gene1
    ALL_GENES$human_gene_name2[i] <- blastp_row$human_gene2
    ALL_GENES$E_val[i] <- blastp_row$V11
  }
}

#####ALL OLD GO approach

GO_terms <- read.csv("GO.oneline.out", sep = " ", colClasses = c(goid = "character"))
unique_genes <- as.data.frame(unique(GO_terms$qpid))
GO_terms$goid <- gsub("^", "GO:", GO_terms$goid )

newgene_col <- rep("name", nrow(unique_genes))
geneGO_col <- rep("unknown", nrow(unique_genes))

for(i in 1:nrow(unique_genes)){
  a <- as.character(unique_genes[i,])
  a_df <- subset(GO_terms, GO_terms$qpid == a)
  if(nrow(a_df)>0){
    newstring <- as.character(a_df$goid)
    test <- paste(newstring, collapse = ";")
  }
  newgene_col[i] <- a
  geneGO_col[i] <- test
}

GO_input <- as.data.frame(cbind(newgene_col, geneGO_col))
colnames(GO_input) <- c("V1", "V2")

write.table(GO_input,file="GO_enrich_allgenes_input.csv",sep="\t",col.names = F,row.names = F, quote = F)

#now filter for only genes in our sampling
output_stats_gene_list <- as.vector(gsub("^", "rna-", output_stats$gene_name))
GO_input_filt <- subset(GO_input, GO_input$V1 %in% output_stats_gene_list)

write.table(GO_input_filt,file="GO_enrich_only_genes_that_ran_input.csv",sep="\t",col.names = F,row.names = F, quote = F)

new_gene_list <- as.vector(gsub("rna-", "", GO_input_filt$V1))

#now make value input
output_stats_filt <- subset(output_stats, output_stats$gene_name %in% new_gene_list)
output_stats_filt$gene_name <- gsub("^", "rna-", output_stats_filt$gene_name)
output_stats_filt$permuted_p <- c()

for(i in 1:nrow(output_stats_filt)){
  gene_numb_value <- output_stats_filt$gene_numb[i]
  new_df_df <- subset(new_all_df, new_all_df$gene == gene_numb_value)
  output_stats_filt$permuted_p[i] <- new_df_df$permuted_P
}

outliers_all <- subset(output_stats_filt, output_stats_filt$P < 0.05)
outliers_all <- subset(outliers_all, outliers_all$permuted_p < 0.05)
GO_test_genes_faster <- subset(outliers_all, outliers_all$Rho > 0)
write.table(GO_test_genes_faster$gene_name,file="GO_test_sig_genes_faster.csv",sep=",",col.names = T,row.names = F, quote = F)

GO_test_genes_slower <- subset(outliers_all, outliers_all$Rho < 0)
write.table(GO_test_genes_slower$gene_name,file="GO_test_sig_genes_slower.csv",sep=",",col.names = T,row.names = F, quote = F)


GO_test_RHO <- as.data.frame(cbind(output_stats_filt$gene_name, output_stats_filt$Rho))
write.table(GO_test_RHO,file="GO_test_RHO.csv",sep=",",col.names = T,row.names = F, quote = F)

GO_test_P <- as.data.frame(cbind(output_stats_filt$gene_name, output_stats_filt$P))
write.table(GO_test_P,file="GO_test_P.csv",sep=",",col.names = T,row.names = F, quote = F)

GO_test_permP <- as.data.frame(cbind(output_stats_filt$gene_name, output_stats_filt$permuted_p))
write.table(GO_test_permP,file="GO_test_permP.csv",sep=",",col.names = T,row.names = F, quote = F)

GO_test_bin <- GO_test_permP
GO_test_bin$bin <- 0

for(i in 1:nrow(GO_test_bin)){
  if(GO_test_bin$V2[i]<0.05)
    GO_test_bin$bin[i] <- 1
}

GO_test_bin_done <- as.data.frame(cbind(GO_test_bin$V1, GO_test_bin$bin))
write.table(GO_test_bin_done,file="GO_test_bin_done.csv",sep=",",col.names = T,row.names = F, quote = F)

testy <- read.csv(file = "../regulatory/DE.out", sep = "\t")

uniq_testy <- as.data.frame(unique(testy$qpid))

geany <- rep("name", nrow(uniq_testy))
geanyname <- rep("unknown", nrow(uniq_testy))

for(i in 1:nrow(uniq_testy)){
  a <- as.character(uniq_testy[i,])
  a_df <- subset(testy, testy$qpid == a)
  if(nrow(a_df)>0){
    newstring <- as.character(a_df$genename[1])
  }
  geany[i] <- a
  geanyname[i] <- newstring
}

newtesty <- as.data.frame(cbind(geany, geanyname))
newtesty$count <- nchar(newtesty$geanyname)
newtesty_filt <- subset(newtesty, newtesty$count>0)


#now for GSEA
#read in files from blair and want to append the correc rho value
GSEA_background <- read.csv("GO_enrich_only_genes_that_ran_input_danioInfo (1).tsv", sep = "\t")
GSEA_background$rho <- "NA"

for(i in 1:nrow(GSEA_background)){
  g_o_i <- as.character(GSEA_background$xhel_id[i])
  g_o_i_gene_df <- subset(output_stats_filt, as.character(output_stats_filt$gene_name) == g_o_i)
  GSEA_background$rho[i] <- as.character(g_o_i_gene_df$Rho)
}

GSEA_background$rho <- as.numeric(GSEA_background$rho)

GSEA_background$first_danio_id_new <- gsub("\\..*", "", GSEA_background$first_danio_id)

GSEA_background_for_webgestalt <- as.data.frame(cbind(GSEA_background$first_danio_id_new, GSEA_background$rho))
GSEA_background_for_webgestalt$V2 <- as.numeric(GSEA_background_for_webgestalt$V2)
GSEA_background_for_webgestalt <- GSEA_background_for_webgestalt[order(GSEA_background_for_webgestalt$V2),]

write.table(GSEA_background_for_webgestalt,file="GSEA_background_for_webgestalt.txt",sep="\t",col.names = T,row.names = F, quote = F)

a <- as.data.frame(unique(GSEA_background$xhel_id))
a <- a[order(a$`unique(GSEA_background$xhel_id)`),]
b <- as.data.frame(unique(output_stats_filt$gene_name))
b <- b[order(b$`unique(output_stats_filt$gene_name)`),]

ab <- cbind(a, b)

gestalt_background_genes <- newtesty_filt$geanyname
write.table(gestalt_background_genes,file="gestalt_background_genes.csv",sep=",",col.names = T,row.names = F, quote = F)

gestalt_sig_high_rho <- 
  gestalt_sig_low_rho <- 
  
  
  ortho_output <- read.csv("RER.protein__v__Danio_rerio.GRCz11.pep.all.tsv", sep = "\t")
library(stringr)

# Count the number of 'a's in each element of string
ortho_output$homot_term_count <- str_count(ortho_output$Danio_rerio.GRCz11.pep.all, ",")
one_term_ortho <- subset(ortho_output, ortho_output$homot_term_count == 0)
