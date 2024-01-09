#make input for regulatory RERconverge run
setwd("/Users/rishide-kayne/Dropbox/RishiMac2/Kelley_Lab/Kerry_project/regulatory/")
#want to have 'regions' input and then a filename for the region

part1 <- read.csv("putativeEnhancerPeaks_xmac_10.27.23.merged.bed", head=FALSE, sep = "\t")
part2 <- read.csv("promoters_500bp_xmac.merged.bed", head=FALSE, sep = "\t")

part2$No_of_elements <- nchar(gsub("[^|]","",part2$V4))

part2 <- subset(part2, part2$No_of_elements < 1)

newpart1 <- part1[,1:3]
newpart2 <- part2[,1:3]

full_parts <- rbind(newpart1, newpart2)

ord1full_parts <- full_parts[order(full_parts$V2),]
ord2full_parts <- ord1full_parts[order(ord1full_parts$V1),]


newpart1$region <- paste(paste(newpart1$V1, newpart1$V2, sep = ":"), newpart1$V3, sep = "-")
newpart1$filename <- paste(paste(newpart1$V1, newpart1$V2, sep = "_"), newpart1$V3, sep = "_")

ord1newpart1 <- newpart1[order(newpart1$V2),]
ord2newpart1 <- ord1newpart1[order(ord1newpart1$V1),]

enhancer_file <- ord2newpart1[,4:5]

write.table(enhancer_file, file = "enhancers_rer_input.csv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


newpart2$region <- paste(paste(newpart2$V1, newpart2$V2, sep = ":"), newpart2$V3, sep = "-")
newpart2$filename <- paste(paste(newpart2$V1, newpart2$V2, sep = "_"), newpart2$V3, sep = "_")

ord1newpart2 <- newpart2[order(newpart2$V2),]
ord2newpart2 <- ord1newpart2[order(ord1newpart2$V1),]

promoter_file <- ord2newpart2[,4:5]

write.table(promoter_file, file = "promoters_rer_input.csv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


