 ############################################################################################################
#  This script clusters/projects a set of L1000 samples using NMF
#  in order to define "components" representing transcriptional responses.
#  It also creates a 3D biplot of the H matrix and makes macthes of each perturbation vs. the rest.
#  Oct 24, 2013
############################################################################################################

  # Pre-processing: rank normalization for CC datasets
  
   source("/xchip/cogs/hogstrom/analysis/pablos_NMF_analysis/TA/CNMF.4.R")
   # parse args
   # path1 <- "/xchip/cogs/projects/NMF/NMF_parameter_evaluation2/A549_c9_lm_epsilon"
   # prefix1 <- "clique_compound_classes_n612x978"
   args <- commandArgs(trailingOnly = TRUE)
   path1 <- args[1]
   prefix1 <- args[2]
   
   MSIG.Preprocess.Dataset(
      input.ds            = paste(path1,"/",prefix1,".gct",sep=""),
      output.ds           = paste(path1,"/",prefix1,".NORM.gct",sep=""),
      normalization       = 6)   # replace values with rank/total genes.

   ##  CMAP Compound Classes ----------------------------------------------------------------------------------------------------------  
      # Input Files
      L1000.file   <- paste(path1,"/",prefix1,".NORM.gct",sep="")
      # annot.file   <- paste(path1,"/",cell1,"_top_intra_connecting_compound_classes.v2.txt",sep="")
      # annot.file   <- paste(path1,"/clique_compound_classes.v2.txt",sep="")
      # annot.file   <- paste(path1,"/shRNA_drug_target_genes.v2.txt",sep="")
      # Parameters
      k.comp       <- 20    # Optimal number of components: 9, 20
      name.column  <- 1     # Column # in annot.file containing the perturbation name
      class.column <- 7     # Column # in annot.file containing the class or category name
      use.prefix   <- F     # Use only prefix before "_" to find association between perturbation names in Input File vs. annot.file
      n.top        <- 35         # Number of top/highest IC associations to display in heatmap
      n.bottom     <- 10      # Number of bottom/lowest IC associations to display in heatmap
      # Output Files
      pdf.file     <- paste(path1,"/",prefix1, ".k", k.comp, ".pdf", sep="")
      W.file       <- paste(path1,"/",prefix1, ".W.k", k.comp, ".gct", sep="")
      H.file       <- paste(path1,"/",prefix1, ".H.k", k.comp, ".gct", sep="")
      MI.in.file       <- paste(path1,"/",substr(prefix1,1,23), ".MI.input_space.gct", sep="")
      MI.k.file       <- paste(path1,"/",substr(prefix1,1,23), ".MI.k", k.comp, ".gct", sep="")
      movie.file   <- paste(path1,"/",prefix1, ".Biplot.Movie.k", k.comp, ".gct", sep="")
   ##  --------------------------------------------------------------------------------------------------------------------------------------------

   # Libraries

   # pdf(file=pdf.file, height=8.5, width=11)

   source("/xchip/cogs/hogstrom/analysis/pablos_NMF_analysis/TA/CNMF.4_lh.R")
   source("/xchip/cogs/hogstrom/analysis/pablos_NMF_analysis/TA/OPAM.library.v7.R")   
   source("/xchip/cogs/hogstrom/analysis/pablos_NMF_analysis/TA/FS.library.v8.6.R")     
   library(RColorBrewer)
   library(MASS)
   library(smacof)
   library(quadprog)

   # Define color map

   mycol.pinko <- vector(length=512, mode = "numeric")
   for (k in 1:256) mycol.pinko[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
   for (k in 257:512) mycol.pinko[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
   mycol.pinko <- rev(mycol.pinko)
   ncolors <- length(mycol.pinko)

   mycol <- vector(length=512, mode = "numeric")
   for (k in 1:512) mycol[k] <- rgb(470, 470 - 470*(k - 1)/512, 511, maxColorValue=511)
   ncolors <- length(mycol)

   col.classes <- brewer.pal(7, "Set1")

   # Read L1000 expression dataset

   dataset.1 <- MSIG.Gct2Frame(filename = L1000.file)
   m.2 <- data.matrix(dataset.1$ds)
   dim(m.2)

   # Read annotation file

   # annot.table <- read.table(annot.file, header=F, sep="\t", skip=0, colClasses = "character")
   # annot.table <- read.table(annot.file, header=T, sep="\t", skip=0, colClasses = "character")
   # gene.table <- annot.table[, name.column]
   # pathway.table <- annot.table[, class.column]
   # gene.set <- vector(length=ncol(m.2), mode="character")
   # if (use.prefix == T) {
   #    for (i in 1:ncol(m.2)) {
   #       gene.set[i] <- strsplit(colnames(m.2)[i], split="_")[[1]]
   #    }
   # } else {
   #    gene.set <- colnames(m.2)
   # }
   # locs <- match(gene.set, gene.table)
   # gene.class <- pathway.table[locs]
   # for (k in 1:length(gene.class)) gene.class[k] <- substr(gene.class[k], 1, 10)
     
   # table(gene.class)
   # all.classes <- unique(gene.class)

   # Visualize input matrix

   # heatmap(m.2, scale="row", col=mycol, margins=c(15, 15), cexRow=0.10, cexCol=0.5, main="A Matrix", xlab = " ", ylab= " ")

   # Perform NMF clustering

   NMF.out <- NMF.div(V = m.2, k = k.comp, maxniter = 2000, seed = 123, stopconv = 40, stopfreq = 10)
   W <- NMF.out$W
   H_orig <- NMF.out$H
   colnames(W) <- paste("c", seq(1, k.comp), sep="")
   row.names(W) <- row.names(m.2)

   # Obtain H via W: Project original and additional dataset using non-negative solver

   H <- matrix(0, nrow=k.comp, ncol= ncol(m.2), dimnames=list(colnames(W), colnames(m.2)))
   for (i in 1:ncol(H)) H[, i] <- nnls.fit(W, m.2[, i], wsqrt=1, eps=0, rank.tol=1e-07)

   # Save W and H matrices

   write.gct.2(gct.data.frame = W, descs = row.names(W), filename = W.file)
   write.gct.2(gct.data.frame = H, descs = row.names(H), filename = H.file)

   # Visualize W

   # hc <- hclust(dist(t(W)), "complete")
   # d1.W <- as.dendrogram(hc)
   # hc2 <- hclust(dist(W), "complete")
   # d2 <- as.dendrogram(hc2)
 
   # heatmap(W, Colv=d1.W, Rowv = d2,  scale="row", col=mycol, margins=c(15, 15), cexRow=0.10, cexCol=0.5, main="Sorted W Matrix", xlab = "Meta-Genes", ylab= " ")

  # Visualize H 

   # hc <- hclust(dist(t(H)), "complete")
   # d1 <- as.dendrogram(hc)
   # hc2 <- hclust(dist(H), "complete")
   # d2 <- as.dendrogram(hc2)
   # heatmap(t(H), Colv=d1.W, Rowv = d1,  scale="row", col=mycol, margins=c(15, 15), cexRow=0.35, cexCol=0.55, main="Sorted H Matrix", xlab = "Meta-Genes", ylab= " ")

   # Signatures association plot
       
   nf <- layout(matrix(c(1), 1, 1, byrow=T), 1, 1, FALSE)
         
   MI.matrix <- matrix(0, nrow=ncol(m.2), ncol=ncol(m.2), dimnames = list(colnames(m.2), colnames(m.2)))
   for (i in 1:ncol(m.2)) for (j in 1:ncol(m.2)) MI.matrix[i, j] <- IC.v1(m.2[,i], m.2[,j])
   dist.matrix <- as.dist(1 - MI.matrix)
   HC <- hclust(dist.matrix, method="complete")
   MI.matrix <- MI.matrix[HC$order, HC$order]
   write.gct.2(gct.data.frame = MI.matrix, descs = row.names(MI.matrix), filename = MI.in.file)
   # gc <- gene.class[HC$order]
   max.v <- max(max(MI.matrix), -min(MI.matrix))
   MI.matrix <- ceiling(ncolors * (MI.matrix - (- max.v))/(1.001*(max.v - (- max.v))))
   # V <- apply(MI.matrix, MARGIN=2, FUN=rev)
   # par(mar = c(2, 8, 14, 8))
   # image(1:dim(V)[2], 1:dim(V)[1], t(V), main = "IC Association Matrix using landmark genes", zlim = c(0, ncolors), col=mycol,
   #       axes=FALSE, sub = "", xlab= "", ylab="")
   # mtext(row.names(V), at=1:nrow(V), side = 2, cex=0.30, col=rev(col.classes[match(gc, all.classes)]), line=0, las=1, font=2, family="")
   # mtext(colnames(V), at=1:ncol(V), side = 3, cex=0.30, col=col.classes[match(gc, all.classes)], line=0, las=3, font=2, family="")

   MI.matrix <- matrix(0, nrow=ncol(H), ncol=ncol(H), dimnames = list(colnames(H), colnames(H)))
   for (i in 1:ncol(H)) for (j in 1:ncol(H)) MI.matrix[i, j] <- IC.v1(H[,i], H[,j])
   dist.matrix <- as.dist(1 - MI.matrix)
   HC <- hclust(dist.matrix, method="complete")
   MI.matrix <- MI.matrix[HC$order, HC$order]
   write.gct.2(gct.data.frame = MI.matrix, descs = row.names(MI.matrix), filename = MI.k.file)
   # gc <- gene.class[HC$order]
   max.v <- max(max(MI.matrix), -min(MI.matrix))
   MI.matrix <- ceiling(ncolors * (MI.matrix - (- max.v))/(1.001*(max.v - (- max.v))))
   # V <- apply(MI.matrix, MARGIN=2, FUN=rev)
   # par(mar = c(2, 8, 14, 8))
   # image(1:dim(V)[2], 1:dim(V)[1], t(V), main = "IC Association Matrix using components", zlim = c(0, ncolors), col=mycol,
   #       axes=FALSE, sub = "", xlab= "", ylab="")
   # mtext(row.names(V), at=1:nrow(V), side = 2, cex=0.30, col=rev(col.classes[match(gc, all.classes)]), line=0, las=1, font=2, family="")
   # mtext(colnames(V), at=1:ncol(V), side = 3, cex=0.30, col=col.classes[match(gc, all.classes)], line=0, las=3, font=2, family="")

   # Merge H matrices

   # Normalize H

   # Ht <- t(H)
   # sum.H <- apply(Ht, MARGIN=1, FUN=sum)
   # for (i in 1:nrow(Ht)) {
   #    Ht[i,] <- Ht[i,]/sum.H[i]
   # }

   # for (i in 1:nrow(Ht)) {
   #    nf <- layout(matrix(c(1, 2, 3, 4), 4, 1, byrow=T), 1, c(4, ceiling(n.top/2) + 2, ceiling(n.bottom/2) + 2, 3),  FALSE)
   #    ind <- order(Ht[i,], decreasing=T)
   #    Ht.temp <- Ht[, ind]
   #    mutinf.m <- NULL
   #    for (k in 1:nrow(Ht.temp)) mutinf.m <- c(mutinf.m, IC.v1(Ht.temp[i,], Ht.temp[k,]))
   #    ind <- order(mutinf.m, decreasing=T)
   #    mutinf.m <- signif(mutinf.m[ind], 3)
   #    Ht.temp <- Ht.temp[ind,]
   #    gene.class.temp <- gene.class[ind]
     
   #    # Target profile
    
   #    V <- Ht.temp[1,]
   #    gc <- gene.class.temp[1]
   #    V <- V*ncolors
   #    par(mar = c(1, 20, 4, 6))    
   #    image(1:length(V), 1:1, as.matrix(V), zlim = c(0, ncolors), col=mycol, axes=FALSE, main=row.names(Ht)[i], sub = "", cex.main = 2, xlab= "", ylab="")     
   #    mtext(paste(row.names(Ht)[i], " (", gc, ") ", sep=""), at=1, side = 2, cex=0.50, col=col.classes[match(gc, all.classes)], line=0, las=1, font=2, family="")
   #    axis(4, at=1+0.2, labels="        IC  ", adj= 0.5, tick=FALSE, las = 1, cex.axis=0.85, line=-1, font=2, family="")
   #    axis(4, at=1, labels="   (Information ", adj= 0.5, tick=FALSE, las = 1, cex.axis=0.75, line=-1, font=2, family="")
   #    axis(4, at=1-0.2, labels="    Coefficient) ", adj= 0.5, tick=FALSE, las = 1, cex.axis=0.75, line=-1, font=2, family="")     
   #    axis(3, at=1:length(V), labels=colnames(Ht.temp), adj= 0.5, tick=FALSE, las = 1, cex.axis=ifelse(k.comp >= 25, 0.6, 0.85),
   #         font.axis=1, line=-1, font=2, family="")
     
   #    # HIGHEST IC ASSOCIATIONS

   #    V <- Ht.temp[2:(n.top+1),]
   #    gc <- gene.class.temp[2:(n.top+1)]
   #    row.names(V) <- paste(row.names(V), " (", gc, ") ", sep="")
   #    muti <- mutinf.m[2:(n.top+1)]
   #    V <- V*ncolors
   #    V <- apply(V, MARGIN=2, FUN=rev)
   #    gc <- rev(gc)
   #    muti <- rev(muti)
   #    par(mar = c(1, 20, 3, 6))
   #    image(1:dim(V)[2], 1:dim(V)[1], t(V), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="HIGHEST IC ASSOCIATIONS", sub = "", xlab= "", ylab="")     
   #    mtext(row.names(V), at=1:nrow(V), side = 2, cex=0.50, col=col.classes[match(gc, all.classes)], line=0, las=1, font=2, family="")
   #    axis(4, at=1:nrow(V), labels=muti, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.85, font.axis=1, line=0, font=2, family="")
   #    axis(3, at=1:ncol(V), labels=colnames(V), adj= 0.5, tick=FALSE, las = 1, cex.axis=ifelse(k.comp >= 25, 0.6, 0.85),
   #         font.axis=1, line=-1, font=2, family="")

   #    # LOWEST IC ASSOCIATIONS
     
   #    V <- Ht.temp[seq(nrow(Ht.temp) - n.bottom + 1, nrow(Ht.temp)),]
   #    gc <- gene.class.temp[seq(nrow(Ht.temp) - n.bottom + 1, nrow(Ht.temp))]
   #    row.names(V) <- paste(row.names(V), " (", gc, ") ", sep="")
   #    muti <- mutinf.m[seq(nrow(Ht.temp) - n.bottom + 1, nrow(Ht.temp))]
   #    V <- V*ncolors
   #    V <- apply(V, MARGIN=2, FUN=rev)
   #    gc <- rev(gc)
   #    muti <- rev(muti)     
   #    par(mar = c(2, 20, 3, 6))     
   #    image(1:dim(V)[2], 1:dim(V)[1], t(V), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="LOWEST IC ASSOCIATIONS", sub = "", xlab= "", ylab="")     
   #    mtext(row.names(V), at=1:nrow(V), side = 2, cex=0.50, col=col.classes[match(gc, all.classes)], line=0, las=1, font=2, family="")    
   #    axis(4, at=1:nrow(V), labels=muti, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.85, font.axis=1, line=0, font=2, family="")
   #    axis(3, at=1:ncol(V), labels=colnames(V), adj= 0.5, tick=FALSE, las = 1, cex.axis=ifelse(k.comp >= 25, 0.6, 0.85),
   #         font.axis=1, line=-1, font=2, family="")

   #    par.mar <- par("mar")
   #    par(mar = c(4, 45, 1, 5))
   #    leg.set <- seq(0, 1, 0.01)
   #    image(1:length(leg.set), 1:1, as.matrix(leg.set), zlim=c(0, 1), col=mycol, axes=FALSE, main=paste("Legend"), sub = "", xlab= "", ylab="",font=2, family="")
   #    ticks <- seq(0, 1, 0.1)
   #    tick.cols <- rep("black", 5)
   #    tick.lwd <- 2
   #    locs <- NULL
   #    for (k in 1:length(ticks)) locs <- c(locs, which.min(abs(ticks[k] - leg.set)))
   #    axis(1, at=locs, labels=ticks, adj= 0.5, tick=T, cex=0.7, cex.axis=1, line=0, font=2, family="")
   #    mtext("Component Amplitude", cex=0.85, side = 1, line = 2.5, outer=F)     
   #    par(mar = par.mar)

   # }
   dev.off()

   # ### Biplot MDS projection

   # library(rgl)
   # library(bpca)

   # row.names <- row.names(H)
   # col.names <- colnames(H)
   # row.names(H) <- paste("    ", row.names, sep="")
   # colnames(H) <- paste("    ", col.names, sep="")

   # obj.col <- col.classes[match(gene.class, all.classes)]
   # colnames(H) <- paste(colnames(H), " (", gene.class, ")", sep="")

   # bpca.1 <- bpca(t(H), method='hj', lambda.end=3)

   # open3d(windowRect = c(10, 10, 1200, 1200), zoom=0.35)

   # plot(bpca.1, ref.lines=F, rgl.use=TRUE, var.col='black', var.factor=0.35, var.cex=1.2, 
   #       obj.names=T, obj.cex=0.70, obj.col=obj.col, simple.axes=FALSE, box=FALSE)

   # s <- strsplit(movie.file, split="/")
   # movie.name <- s[[1]][length(s[[1]])]
   # dir <- s[[1]][seq(1, length(s[[1]])-1)]

   # my.movie3d(spin3d(axis=c(1,1,1), rpm=2), duration=5, fps = 10, movie = movie.name, dir = paste(dir, collapse="/"),
   #         convert = TRUE, clean = TRUE, verbose=TRUE, top = TRUE, type = "gif", startTime = 0) 
