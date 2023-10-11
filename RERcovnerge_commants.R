library(RERconverge)
toyTrees=readTrees("/hb/home/rdekayne/11_RERconverge/test_run/NucGeneTrees.trees")

mamRERw = getAllResiduals(toyTrees, transform = "sqrt", weighted = T, scale = T)

saveRDS(mamRERw, file="mamRERw.rds") 
newmamRERw = readRDS("mamRERw.rds")

#make average and gene tree plots
pdf(file="test_plot1.pdf", width=8, height=8)
par(mfrow=c(1,2))
avgtree=plotTreeHighlightBranches(toyTrees$masterTree, hlspecies=c("pso", "lluvia", "esperanza", "banos", "gloria", "FL15-032", "CR15-031", "MX15-627", "MX15-706", "FL15-025", "MX15-650", "MX15-747"), hlcols=c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue"), main="Average tree") #plot average tree
gene_12686=plotTreeHighlightBranches(toyTrees$trees$gene_12686, hlspecies=c("pso", "lluvia", "esperanza", "banos", "gloria", "FL15-032", "CR15-031", "MX15-627", "MX15-706", "FL15-025", "MX15-650", "MX15-747"), hlcols=c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue"), main="BEND3 tree") #plot individual gene tree
dev.off()

#plot RERs
pdf(file="test_plot2.pdf", width=8, height=8)
par(mfrow=c(1,1))
phenvExample <- foreground2Paths(c("pso", "lluvia", "esperanza", "banos", "gloria", "FL15-032", "CR15-031", "MX15-627", "MX15-706", "FL15-025", "MX15-650", "MX15-747"),toyTrees,clade="terminal")
plotRers(mamRERw,"gene_12686",phenv=phenvExample) #plot RERs
dev.off()

#plot RERs as tree
pdf(file="test_plot3.pdf", width=8, height=8)
par(mfrow=c(1,1))
gene_12686rers = returnRersAsTree(toyTrees, mamRERw, "gene_12686", plot = TRUE, phenv=phenvExample) #plot RERs
dev.off()

#strwrap(gsub(":",write.tree(gene_12686rers),replacement=": "))

#write.tree(gene_12686rers, file='gene_12686rersRER.nwk')

multirers = returnRersAsTreesAll(toyTrees,mamRERw)
write.tree(multirers, file='toyRERs.nwk', tree.names=TRUE)

#visualize RERs along branches as a heatmap
pdf(file="test_plot4.pdf", width=20, height=20)
newgene_12686rers = treePlotRers(treesObj=toyTrees, rermat=mamRERw, index="gene_12686", type="c", nlevels=9, figwid=2)
dev.off()


######################
BINARY TRAIT ANALYSIS
######################

marineextantforeground = c("pso", "lluvia", "esperanza", "banos", "gloria", "FL15-032", "CR15-031", "MX15-627", "MX15-706", "FL15-025", "MX15-650", "MX15-747")
marineb2a = foreground2Tree(marineextantforeground, toyTrees, clade="ancestral")
marineb2b = foreground2Tree(marineextantforeground, toyTrees, clade="terminal")
marineb2c = foreground2Tree(marineextantforeground, toyTrees, clade="all")
marineb2d = foreground2Tree(marineextantforeground, toyTrees, clade="all", weighted = TRUE)

pdf(file="test_plot5.pdf", width=20, height=20)
par(mfrow=c(2,2))
plot(marineb2a,cex=0.6,main="ancestral")
plot(marineb2b,cex=0.6,main="terminal")
plot(marineb2c,cex=0.6,main="all unweighted", x.lim=c(0,2.5))
labs2c = round(marineb2c$edge.length,3)
labs2c[labs2c==0] = NA
edgelabels(labs2c, col = 'black', bg = 'transparent', adj = c(0.5,-0.5),cex = 0.4,frame='n')
plot(marineb2d,cex=0.6,main="all weighted", x.lim=c(0,2.5))
labs2d = round(marineb2d$edge.length,3)
labs2d[labs2d==0] = NA
edgelabels(labs2d, col = 'black', bg = 'transparent', adj = c(0.5,-0.5),cex = 0.4,frame='n')
dev.off()

phenvMarine2b=tree2Paths(marineb2b, toyTrees)
corMarine=correlateWithBinaryPhenotype(mamRERw, phenvMarine2b, min.sp=10, min.pos=2,weighted="auto")
head(corMarine[order(corMarine$P),])
allvalues_test_DF <- as.data.frame(corMarine[order(corMarine$P),])
write.table(allvalues_test_DF, "allvalues_test.csv", sep = ',', quote = FALSE)

pdf(file="test_plot6.pdf", width=8, height=8)
plotRers(mamRERw,"gene_4539",phenv=phenvMarine2b) 
dev.off()

pdf(file="test_plot7.pdf", width=8, height=8)
hist(corMarine$P, breaks=15, xlab="Kendall P-value", main="P-values for correlation between all genes and peno")
dev.off()
