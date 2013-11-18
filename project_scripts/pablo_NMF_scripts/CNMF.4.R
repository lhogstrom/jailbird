
consensusNMF.2 <- function(input.ds, k.init, k.final, num.clusterings, maxniter, error.function, rseed=123456789, directory = "", stopconv = 40, stopfreq = 10, non.interactive.run = F, doc.string = "", ...) {

#
#  GenePattern Methodology for:
#
#  Metagenes and Molecular Pattern Discovery using Matrix Factorization
#  Jean-Philippe Brunet, Pablo Tamayo, Todd. R. Golub, and Jill P. Mesirov
# 
#  Author:  Pablo Tamayo (tamayo@genome.wi.mit.edu)
#
#  Based on the original matlab version written by Jean-Philippe Brunet (brunet@broad.mit.edu) and
#  with additional contributions from: Ted Liefeld (liefeld@broad.mit.edu)   
#  Date:  November 27, 2003
#
#  Last change March 3, 2005: modifications to make the output more readable.
#
#  Execute from an R console window with this command:
#  source("<this file>", echo = TRUE)
#  E.g. someoutput <- mynmf2(input.ds="c:\\nmf\\all_aml.res",k.init=2,k.final=5,num.clusterings=20,maxniter=500) 
#
#  For details on the method see:
#
#  Proc. Natl. Acad. Sci. USA 2004 101: 4164-4169
#  http://www.broad.mit.edu/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=89
#
#  Input parameters
#
#   input.ds
#                       input gene expression dataset in GCT or RES format
#   k.init
#                       initial value of k
#   k.final
#                       final value of k
#   num.clusterings
#                       number of NMF clusterings to build consensus matrix
#   maxniter
#                       maximum number of NMF iterations
#   error.function
#                       NMF error function: "divergence" of "euclidean"
#   rseed
#                       random number generator seed
#   directory
#                       file directory where to store the result files
#   stopconv
#                       how many no change checks are needed to stop NMF iterations (convergence)
#   stopfreq
#                       frequency (NMF iterations) of "no change" checks 
#   non.interactive.run 
#                       flag controlling if the plots are produced interatively (Rgui and saved) or only saved in files
#   doc.string
#                       prefix to be added to the output files
#
#  Output files are (prefix with "doc.string")
#
#   params.txt 
#                       run parameters and time of execution
#   membership.gct		
#			membership results for samples at all values of K
#   cophenetic.txt 
#			cophenetic values for each K
#   cophenetic.plot.jpeg
#			plot of cophenetic for each value of K		
#   consensus.k#.gct (for each value of K)
#			consensus matrix for k=#
#   consensus.plot.k#.jpeg (for each value of K)
#			plot of consensus matrix for k=#
#   graphs.k#.jpeg (for each value of K)

# save input parameters

filename <- paste(directory, doc.string, ".params.txt", sep="", collapse="")  

time.string <- as.character(as.POSIXlt(Sys.time(),"GMT"))
write(paste("Run of NMF on ", time.string), file=filename)

write(paste("input.ds =", input.ds, sep=" "), file=filename, append=T) 
write(paste("k.init = ", k.init, sep=" "), file=filename, append=T) 
write(paste("k.final =", k.final, sep=" "), file=filename, append=T) 
write(paste("num.clusterings =", num.clusterings, sep=" "), file=filename, append=T) 
write(paste("maxniter =", maxniter, sep=" "), file=filename, append=T) 
write(paste("error.function =", error.function, sep=" "), file=filename, append=T) 
write(paste("rseed =", rseed, sep=" "), file=filename, append=T) 
write(paste("directory =", directory, sep=" "), file=filename, append=T) 
write(paste("stopconv =", stopconv, sep=" "), file=filename, append=T) 
write(paste("stopfreq =", stopfreq, sep=" "), file=filename, append=T)
write(paste("non.interctive.run =", non.interactive.run, sep=" "), file=filename, append=T) 
write(paste("doc.string =", doc.string, sep=" "), file=filename, append=T) 


k.init<-as.integer(k.init)
k.final<-as.integer(k.final)
num.clusterings<-as.integer(num.clusterings)
n.iter<-as.integer(maxniter)
if (!is.na(rseed)){
     seed <- as.integer(rseed)
}


# library(mva)
# library(MASS)
# library(GenePattern)

D <- CNMF.read.dataset(input.ds)
A <- data.matrix(D)

# Threshold negative values to small quantity 

eps <- .Machine$double.eps
A[A < 0] <- eps



cols <- length(A[1,])
rows <- length(A[,1])

col.names <- names(D)

num.k <- k.final - k.init + 1

rho <- vector(mode = "numeric", length = num.k)
k.vector <- vector(mode = "numeric", length = num.k)

k.index <- 1

connect.matrix.ordered <- array(0, c(num.k, cols, cols))

filename <- paste(directory, doc.string, ".", "graphs.pdf", sep="", collapse="")
pdf(file=filename, width = 9, height = 11)

for (k in k.init:k.final) { 

   nf <- layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow=T), c(1, 1), c(1, 1, 1, 1), TRUE)
   assign <- matrix(0, nrow = num.clusterings, ncol = cols)

   for (i in 1:num.clusterings) {
	  
        print(paste("Computing clustering number=", i, " for k=", k, sep=""))

        if (error.function == "divergence"){
	    NMF.out <- NMF.div(V = A, k = k, maxniter = n.iter, seed = seed + i, stopconv = stopconv, stopfreq = stopfreq)
	} else if (error.function == "euclidean"){
	    NMF.out <- NMF(V = A, k = k, maxniter = n.iter, seed = seed + i, stopconv = stopconv, stopfreq = stopfreq)
	} else {
            stop(paste("Un-supported error function=", error.function, sep=""))
        }
        print(paste(NMF.out$t, " NMF iterations performed", sep=""))

        for (j in 1:cols) { # Find membership
            class <- order(NMF.out$H[,j], decreasing=T)
            assign[i, j] <- class[1]
        }

	if (i == 1) {  # Plot example for first clustering iteration
            H.saved <- NMF.out$H
            sub.string <- paste(doc.string, " k=", k, sep="")
            plot(1:NMF.out$t, NMF.out$error.v[1:NMF.out$t], pch = 20, cex = 1.5, col = 1, xlab="time", ylab="NMF error", sub=sub.string, main=paste("Convergence plot k=", k, " example", sep=""))


            if (rows < 1000) {
               W <- NMF.out$W
            } else {
               W <- NMF.out$W[sample(x = 1:rows, size = 1000),]
            }
            sub.string <- paste(doc.string, " k=", k, sep="")
            CNMF.matrix.abs.plot(W, sub = sub.string, log = F, main = "Example W matrix (orig. order)", ylab = "genes", xlab ="metasamples")
            CNMF.matrix.abs.plot(H.saved, sub = sub.string, log = F, main = "Example H matrix (orig. order)", ylab = "metagenes", xlab ="samples")
            CNMF.metagene.plot(H = H.saved, main = "Metagenes Example (orig. order)", sub = sub.string, xlab = "samples", ylab = "metagenes")

        }

        rm(NMF.out)

     }  ## end  for (i in 1:num.clusterings)

   
     # compute consensus matrix
     connect.matrix <- matrix(0, nrow = cols, ncol = cols)

     for (i in 1:num.clusterings) {
       for (j in 1:cols) {
          for (p in 1:cols) {
             if (j != p) {
                  if (assign[i, j] == assign[i, p]) {
                    connect.matrix[j, p] <- connect.matrix[j, p] + 1
                  } 
              } else {
                    connect.matrix[j, p] <- connect.matrix[j, p] + 1
              }
           }
       }
     }

     connect.matrix <- connect.matrix / num.clusterings

     dist.matrix <- 1 - connect.matrix
     dist.matrix <- as.dist(dist.matrix)
     HC <- hclust(dist.matrix, method="average")

     dist.coph <- cophenetic(HC)
     k.vector[k.index] <- k
     rho[k.index] <- cor(dist.matrix, dist.coph)
     rho[k.index] <- signif(rho[k.index], digits = 4)
   
#     connect.matrix.ordered <- matrix(0, nrow=cols, ncol = cols)

     for (i in 1:cols) {
        for (j in 1:cols) {
           connect.matrix.ordered[k.index, i, j] <- connect.matrix[HC$order[i], HC$order[j]]
         }
     }

     # compute consensus clustering membership

     membership <- cutree(HC, k = k)

     max.k <- max(membership)
     items.names.ordered <- col.names[HC$order]
     membership.ordered <- membership[HC$order]
     results <- data.frame(cbind(membership.ordered, items.names.ordered))

     if (k > k.init){
          all.membership <- cbind(all.membership, membership);
     } else {
          all.membership <- cbind(membership);
     }

     sub.string <- paste(doc.string, " k=", k, sep="")
     CNMF.matrix.abs.plot(connect.matrix.ordered[k.index,,], sub=sub.string, log = F, main = "Ordered Consensus Matrix", 
                          ylab = "samples", xlab ="samples")
     plot(HC, xlab="samples", cex = 0.75, labels = col.names, sub = sub.string, col = "blue", 
          main = paste("Ordered Linkage Tree. Coph=", rho[k.index]))

     matrixGct <- data.frame(connect.matrix.ordered[k.index,,])
     filename <- paste(directory, doc.string, ".", "matrix.k.",k, ".gct", sep="", collapse="")
     CNMF.write.gct.2(matrixGct, descs = "", filename) 

     resultsGct <- data.frame(membership.ordered)
     row.names(resultsGct) <- items.names.ordered
     filename <- paste(directory, doc.string, ".", "consensus.k.",k, ".gct", sep="", collapse="")
     CNMF.write.gct.2(resultsGct, descs = "", filename) 

     H.sorted <- H.saved[,HC$order]
     sub.string <- paste(doc.string, " k=", k, sep="")
     CNMF.matrix.abs.plot(H.sorted, sub = sub.string, log = F, main = "Example H matrix (ordered)", ylab = "metagenes", xlab ="samples")
     CNMF.metagene.plot(H = H.sorted, sub = sub.string, main = "Metagenes Example (ordered)", xlab = "samples", ylab = "metagenes")

     nf <- layout(matrix(c(1), 1, 1, byrow=T), c(1, 1), c(1, 1), TRUE)

     conlabel <- paste("Consensus k =", k, sep=" ", collapse="")

     sub.string <- paste("Consensus matrix k=", k, "; dataset= ", input.ds, sep="")
     CNMF.ConsPlot(connect.matrix.ordered[k.index,,], col.labels = membership.ordered, col.names = items.names.ordered, 
                   main = " ", sub=sub.string, xlab=" ", ylab=" ")

     k.index <- k.index + 1

} # end of loop over k


# Save consensus matrices in one file

  nf <- layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), 4, 4, byrow=T), c(1, 1, 1, 1), c(1, 1, 1, 1), TRUE)

  for (k in 1:num.k) { 
     CNMF.matrix.abs.plot(connect.matrix.ordered[k,,], log = F, main = paste("k=", k.vector[k]), 
                          sub = paste("Cophenetic coef.=", rho[k]), ylab = "samples", xlab ="samples")
  }
   
  y.range <- c(1 - 2*(1 - min(rho)), 1)
  plot(k.vector, rho, main ="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, 
                          xlab = "k", ylab="Cophenetic correlation", type = "n")
  lines(k.vector, rho, type = "l", col = "black")
  points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")

# Write the membership matrix

resultsmembership <- data.frame(all.membership)
row.names(resultsmembership) <- col.names
colnames(resultsmembership) <- paste("k=", seq(k.init, k.final))

print("Membership:")

print(resultsmembership)

filename <- paste(directory, doc.string, ".", "membership", ".txt", sep="", collapse="")

col.names <- paste(colnames(resultsmembership), collapse = "\t")
col.names <- paste("SAMPLE", col.names, sep= "\t")
write(noquote(col.names), file = filename, append = F, ncolumns = length(col.names))
write.table(resultsmembership, file=filename, quote=F, col.names = F, row.names = T, append = T, sep="\t")

  nf <- layout(matrix(c(1), 1, 1, byrow=T), c(1), c(1), TRUE)
y.range <- c(1 - 2*(1 - min(rho)), 1)
plot(k.vector, rho, main ="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, xlab = "k", ylab="Cophenetic correlation", type = "n")
lines(k.vector, rho, type = "l", col = "black")
points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")

filename <- paste(directory, doc.string, ".", "cophenetic.txt", sep="")

xx <- cbind(k.vector, rho)
write(noquote(c("k", "\t", "Cophenetic Coefficient")), file = filename, append = F, ncolumns = 1000)
write.table(xx, file = filename, append = T, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = F)

dev.off()
}


CNMF.read.dataset <- function(file) {
	result <- regexpr(paste(".gct","$",sep=""), tolower(file))
	if(result[[1]] != -1)
		return(CNMF.read.gct(file))
	result <- regexpr(paste(".res","$",sep=""), tolower(file))
	if(result[[1]] != -1)
		return(CNMF.read.res(file))
	stop("Input is not a res or gct file.")	
}

CNMF.matrix.abs.plot <- function(V, axes = F, log = F, norm = T, transpose = T, matrix.order = T, max.v = 1, min.v = 0, main = " ", sub = " ", xlab = " ", ylab = "  ") {
      rows <- length(V[,1])
      cols <- length(V[1,])
      if (log == T) {
         V <- log(V)
      }
      B <- matrix(0, nrow=rows, ncol=cols)
	for (i in 1:rows) {
           for (j in 1:cols) {
                if (matrix.order == T) {
                   k <- rows - i + 1
                } else {
                   k <- i
                }
                if (norm == T) {
                  if ((max.v == 1) && (min.v == 0)) {
                     max.val <- max(V)
                     min.val <- min(V)
                  } else {
		     	   max.val = max.v
                     min.val = min.v
                  }
               }
	     B[k, j] <-  max.val - V[i, j] + min.val
           }
      }
	if (transpose == T) {
	  B <- t(B)
        }
	if (norm == T) {
#            image(z = B, zlim = c(min.val, max.val), axes = axes, col = rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), main = main, sub = sub, xlab = xlab, ylab = ylab)
            image(z = B, zlim = c(min.val, max.val), axes = axes, col = rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75), main = main, sub = sub, xlab = xlab, ylab = ylab)           
      } else {
#            image(z = B, axes = axes, col = rainbow(100, s = 1, v = 0.6, start = 0.1, end = 0.9, gamma = 1), main = main, sub = sub, xlab = xlab, ylab = ylab)
            image(z = B, axes = axes, col = rainbow(100, s = 1, v = 0.6, start = 0.1, end = 0.9), main = main, sub = sub, xlab = xlab, ylab = ylab)         
      }
      return(list(B, max.val, min.val))
    }

CNMF.metagene.plot <- function(H, main = " ", sub = " ", xlab = "samples ", ylab = "amplitude") {
	k <- length(H[,1])
	S <- length(H[1,])
	index <- 1:S
	maxval <- max(H)
        minval <- min(H)
	plot(index, H[1,], xlim=c(1, S), ylim=c(minval, maxval), main = main, sub = sub, ylab = ylab, xlab = xlab, type="n")
	for (i in 1:k) {
	    lines(index, H[i,], type="l", col = i, lwd=2)
        }
}


CNMF.ConsPlot <- function(V, col.labels, col.names, main = " ", sub = " ", xlab=" ", ylab=" ") {

# Plots a heatmap plot of a consensus matrix

     cols <- length(V[1,])
     B <- matrix(0, nrow=cols, ncol=cols)
     max.val <- max(V)
     min.val <- min(V)
     for (i in 1:cols) {
         for (j in 1:cols) {
             k <- cols - i + 1
	     B[k, j] <-  max.val - V[i, j] + min.val
          }
     }

     col.names2 <- rev(col.names)
     col.labels2 <- rev(col.labels)
     D <- matrix(0, nrow=(cols + 1), ncol=(cols + 1))

     col.tag <- vector(length=cols, mode="numeric")
     current.tag <- 0
     col.tag[1] <- current.tag
     for (i in 2:cols) {
        if (col.labels[i] != col.labels[i - 1]) {
             current.tag <- 1 - current.tag
        }
        col.tag[i] <- current.tag
     }
     col.tag2 <- rev(col.tag)
     D[(cols + 1), 2:(cols + 1)] <- ifelse(col.tag %% 2 == 0, 1.02, 1.01)
     D[1:cols, 1] <- ifelse(col.tag2 %% 2 == 0, 1.02, 1.01)
     D[(cols + 1), 1] <- 1.03
     D[1:cols, 2:(cols + 1)] <- B[1:cols, 1:cols]

#     col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), "#BBBBBB", "#333333", "#FFFFFF")
     col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75), "#BBBBBB", "#333333", "#FFFFFF")     
     image(1:(cols + 1), 1:(cols + 1), t(D), col = col.map, axes=FALSE, main=main, sub=sub, xlab= xlab, ylab=ylab)
     for (i in 1:cols) {
         col.names[i]  <- paste("      ", substr(col.names[i], 1, 12), sep="")
         col.names2[i] <- paste(substr(col.names2[i], 1, 12), "     ", sep="")
     }

     axis(2, at=1:cols, labels=col.names2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.50, font.axis=1, line=-1)
     axis(2, at=1:cols, labels=col.labels2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.65, font.axis=1, line=-1)

     axis(3, at=2:(cols + 1), labels=col.names, adj= 1, tick=FALSE, las = 3, cex.axis=0.50, font.axis=1, line=-1)
     axis(3, at=2:(cols + 1), labels=as.character(col.labels), adj = 1, tick=FALSE, las = 1, cex.axis=0.65, font.axis=1, line=-1)

     return()
   }

CNMF.read.res <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in RES format and converts it into an R data frame
#
   header.cont <- readLines(filename, n = 1)
   temp <- unlist(strsplit(header.cont, "\t"))
   colst <- length(temp)
   header.labels <- temp[seq(3, colst, 2)]
   ds <- read.delim(filename, header=F, row.names = 2, sep="\t", skip=3, blank.lines.skip=T, comment.char="", as.is=T)
   colst <- length(ds[1,])
   cols <- (colst - 1)/2
   rows <- length(ds[,1])
   A <- matrix(nrow=rows - 1, ncol=cols)
   A <- ds[1:rows, seq(2, colst, 2)]
   table1 <- data.frame(A)
   names(table1) <- header.labels
   return(table1)
}

CNMF.read.gct <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in GCT format and converts it into an R data frame
#
   ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T)
   ds <- ds[-1]
   return(ds)
}

CNMF.write.gct.2 <- function(gct.data.frame, descs = "", filename) 
{
    f <- file(filename, "w")
    cat("#1.2", "\n", file = f, append = TRUE, sep = "")
    cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
    cat("Name", "\t", file = f, append = TRUE, sep = "")
    cat("Description", file = f, append = TRUE, sep = "")

    colnames <- colnames(gct.data.frame)
    cat("\t", colnames[1], file = f, append = TRUE, sep = "")

    if (length(colnames) > 1) {
       for (j in 2:length(colnames)) {
           cat("\t", colnames[j], file = f, append = TRUE, sep = "")
       }
     }
    cat("\n", file = f, append = TRUE, sep = "\t")

    oldWarn <- options(warn = -1)
    m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] +  2)
    m[, 1] <- row.names(gct.data.frame)
    if (length(descs) > 1) {
        m[, 2] <- descs
    } else {
        m[, 2] <- row.names(gct.data.frame)
    }
    index <- 3
    for (i in 1:dim(gct.data.frame)[2]) {
        m[, index] <- gct.data.frame[, i]
        index <- index + 1
    }
    write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
    close(f)
    options(warn = 0)

}

NMF.div <- function(V, k, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {

        N <- length(V[,1])
        M <- length(V[1,])
        set.seed(seed)
        W <- matrix(runif(N*k), nrow = N, ncol = k)  # Initialize W and H with random numbers
        H <- matrix(runif(k*M), nrow = k, ncol = M)
        VP <- matrix(nrow = N, ncol = M)
        error.v <- vector(mode = "numeric", length = maxniter)
        new.membership <- vector(mode = "numeric", length = M)
        old.membership <- vector(mode = "numeric", length = M)
        no.change.count <- 0
        eps <- .Machine$double.eps
        for (t in 1:maxniter) {
                VP = W %*% H
                W.t <- t(W)
                H <- H * (W.t %*% (V/VP)) + eps
                norm <- apply(W, MARGIN=2, FUN=sum)
                for (i in 1:k) {
                    H[i,] <- H[i,]/norm[i]
                }
                VP = W %*% H
                H.t <- t(H)
                W <- W * ((V/VP) %*% H.t) + eps
                norm <- apply(H, MARGIN=1, FUN=sum)
                for (i in 1:k) {
                    W[,i] <- W[,i]/norm[i]
                }
               error.v[t] <- sum(V * log((V + eps)/(VP + eps)) - V + VP)/(M * N)
               if (t %% stopfreq == 0) {

                    for (j in 1:M) {
                        class <- order(H[,j], decreasing=T)
                        new.membership[j] <- class[1]
                     }
                     if (sum(new.membership == old.membership) == M) {
                        no.change.count <- no.change.count + 1
                     } else {
                        no.change.count <- 0
                     }
                     if (no.change.count == stopconv) break
                     old.membership <- new.membership
               }
        }
        return(list(W = W, H = H, t = t, error.v = error.v))
}

MSIG.Preprocess.Dataset <- function(
   input.ds, 
   output.ds,
   thres = NULL, 
   ceil = NULL, 
   shift = NULL,
   fold = NULL, 
   delta = NULL, 
   normalization = NULL,
   cntrl.genes = NULL) {

   print(c("Running MSIG.Preprocess.Dataset... on:", input.ds))
   print(c("output file:", output.ds))
   print(c("normalization =", normalization))
   
# Read dataset

   dataset <- MSIG.Gct2Frame(filename = input.ds)
   m <- data.matrix(dataset$ds)
   gs.names <- dataset$row.names
   gs.descs <- dataset$descs
   sample.names <- dataset$names

# threshold, ceiling and shift

   if (!is.null(thres)) {
     m[m < thres] <- thres
   }
   if (!is.null(ceil)) {
      m[m > ceil] <- ceil
   }
   if (!is.null(shift)) {
      m <- m + shift
   }

   # identify and save control genes

   if (!is.null(cntrl.genes)) {
      gene.names2 <- intersect(cntrl.genes, gs.names)
      locs <- match(gene.names2, gs.names, nomatch=0)
      msig.cntrl <- m[locs, ]
      msig.cntrl.genes <- gs.names[locs]
      msig.cntrl.descs <- gs.descs[locs]
      m <- m[-locs, ]
      gs.names <- gs.names[-locs]
      gs.descs <- gs.descs[-locs]
    }

   # variation filter

   if ((!is.null(fold)) && (!is.null(delta))) {
      temp <- MSIG.VarFilter(V = m, fold = fold, delta = delta, gene.names = gs.names, gene.descs = gs.descs) 
      m <- temp$V
      gs.names <- temp$new.gene.names
      gs.descs <- temp$new.gene.descs
      dim(m) 
   }

   # restore control genes

   if (!is.null(cntrl.genes)) {
      m <- rbind(m, msig.cntrl)
      gs.names <- c(gs.names, msig.cntrl.genes)
      gs.descs <- c(gs.descs, msig.cntrl.descs)
    }

# normalization

   if (!is.null(normalization)) {
      if (normalization == 1) {
         m <- MSIG.NormalizeCols.Rank(m)
      } else if (normalization == 2) {
         m <- MSIG.NormalizeCols.Rank(m)/length(m[,1])
      } else if (normalization == 3) {
         m <- GSEA.NormalizeCols(m) + 3
         m <- GSEA.Threshold(m, 0.001, 100000) 
      } else if (normalization == 4) {
         m <- MSIG.NormalizeCols.Rank(m)/length(m[,1])
      } else if (normalization == 5) {
         m <- MSIG.NormalizeCols.Rescale(m)
      } else if (normalization == 6) {
         cols <- length(m[1,])
         for (j in 1:cols) {  # column rank normalization from 0 to N - 1
            m[,j] <- rank(m[,j], ties.method = "average") - 1
         }
         m <- 10000*m/(length(m[,1]) - 1)
      } else if (normalization == 7) {
         m <- ((100*MSIG.NormalizeCols.Rank(m))%/%length(m[,1]) + 1)
      } else if (normalization == 8) { 
          row.mean <- apply(m, MARGIN=1, FUN=mean)
          for (i in 1:length(m[,1])) {
             m[i,] <- m[i,] / row.mean[i]
          }
      }
   }
   
   V <- data.frame(m)
   names(V) <- sample.names
   row.names(V) <- gs.names
   write.gct(gct.data.frame = V, descs = gs.descs, filename = output.ds)  

 }

GSEA.Threshold <- function(V, thres, ceil) { 
#
# Threshold and ceiling pre-processing for gene expression matrix
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

        V[V < thres] <- thres
        V[V > ceil] <- ceil
        return(V)
}

GSEA.VarFilter <- function(V, fold, delta, gene.names = "") { 
#
# Variation filter pre-processing for gene expression matrix
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

        cols <- length(V[1,])
        rows <- length(V[,1])
        row.max <- apply(V, MARGIN=1, FUN=max)
               row.min <- apply(V, MARGIN=1, FUN=min)
        flag <- array(dim=rows)
        flag <- (row.max /row.min >= fold) & (row.max - row.min >= delta)
        size <- sum(flag)
        B <- matrix(0, nrow = size, ncol = cols)
        j <- 1
        if (length(gene.names) == 1) {
           for (i in 1:rows) {
              if (flag[i]) {
                 B[j,] <- V[i,]
                 j <- j + 1
               }
           }
        return(B)
        } else {
            new.list <- vector(mode = "character", length = size)
            for (i in 1:rows) {
              if (flag[i]) {
                 B[j,] <- V[i,]
                 new.list[j] <- gene.names[i]
                 j <- j + 1
              }
            }
        return(list(V = B, new.list = new.list))
        }
}

MSIG.VarFilter <- function(V, fold, delta, gene.names = "", gene.descs = "") { 

# Variation filter pre-processing for gene expression matrix

        cols <- length(V[1,])
        rows <- length(V[,1])
        row.max <- apply(V, MARGIN=1, FUN=max)
        row.min <- apply(V, MARGIN=1, FUN=min)
        flag <- array(dim=rows)
        flag <- (row.max /row.min >= fold) & (row.max - row.min >= delta)
        size <- sum(flag)
        B <- matrix(0, nrow = size, ncol = cols)
        j <- 1
        if (length(gene.names) == 1) {
           for (i in 1:rows) {
              if (flag[i]) {
                 B[j,] <- V[i,]
                 j <- j + 1
               }
           }
        return(B)
        } else {
            new.gene.names <- vector(mode = "character", length = size)
            new.gene.descs <- vector(mode = "character", length = size)
            for (i in 1:rows) {
              if (flag[i]) {
                 B[j,] <- V[i,]
                 new.gene.names[j] <- gene.names[i]
                 new.gene.descs[j] <- gene.descs[i]
                 j <- j + 1
              }
            }
        return(list(V = B, new.gene.names = new.gene.names, new.gene.descs = new.gene.descs, locations = flag))
        }
}

GSEA.NormalizeRows <- function(V) { 
#
# Stardardize rows of a gene expression matrix
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

        row.mean <- apply(V, MARGIN=1, FUN=mean)
        row.sd <- apply(V, MARGIN=1, FUN=sd)

        row.n <- length(V[,1])
        for (i in 1:row.n) {
             if (row.sd[i] == 0) {
                  V[i,] <- 0
           } else {
              V[i,] <- (V[i,] - row.mean[i])/row.sd[i]
           }
        }
        return(V)
}

GSEA.NormalizeCols <- function(V) { 
#
# Stardardize columns of a gene expression matrix
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

        col.mean <- apply(V, MARGIN=2, FUN=mean)
               col.sd <- apply(V, MARGIN=2, FUN=sd)
        col.n <- length(V[1,])
        for (i in 1:col.n) {
             if (col.sd[i] == 0) {
                  V[i,] <- 0
           } else {
              V[,i] <- (V[,i] - col.mean[i])/col.sd[i]
           }
        }
        return(V)
}

GSEA.NormalizeCols.Rank <- function(V) { 
#
      cols <- length(V[1,])
      rows <- length(V[,1])
      for (j in 1:cols) {  # column rank normalization
         V[,j] <- rank(V[,j], ties.method = "average")
      }

      return(V)
}


MSIG.NormalizeCols.Rank <- function(V) { 

      cols <- length(V[1,])
      rows <- length(V[,1])
      for (j in 1:cols) {  # column rank normalization
         V[,j] <- rank(V[,j], ties.method = "average")
      }

      return(V)
}

MSIG.NormalizeCols.Rescale <- function(V) { 

      epsilon <- 0.00001
      cols <- length(V[1,])
      for (j in 1:cols) {  # column rank normalization
         max.v <- max(V[,j])
         min.v <- min(V[,j])
         V[,j] <- (V[,j] - min.v + epsilon)/(max.v - min.v)
      }

      return(V)
    }

MSIG.Gct2Frame <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in GCT format and converts it into an R data frame
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

   ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T, na.strings = "")
   descs <- ds[,1]
   ds <- ds[-1]
   row.names <- row.names(ds)
   names <- names(ds)
   return(list(ds = ds, row.names = row.names, descs = descs, names = names))
}


write.gct <- function(gct.data.frame, descs = "", filename) 
{
    f <- file(filename, "w")
    cat("#1.2", "\n", file = f, append = TRUE, sep = "")
    cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
    cat("Name", "\t", file = f, append = TRUE, sep = "")
    cat("Description", file = f, append = TRUE, sep = "")

    names <- names(gct.data.frame)
    cat("\t", names[1], file = f, append = TRUE, sep = "")

    if (length(names) > 1) {
       for (j in 2:length(names)) {
           cat("\t", names[j], file = f, append = TRUE, sep = "")
       }
     }
    cat("\n", file = f, append = TRUE, sep = "\t")

    oldWarn <- options(warn = -1)
    m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] +  2)
    m[, 1] <- row.names(gct.data.frame)
    if (length(descs) > 1) {
        m[, 2] <- descs
    } else {
        m[, 2] <- row.names(gct.data.frame)
    }
    index <- 3
    for (i in 1:dim(gct.data.frame)[2]) {
        m[, index] <- gct.data.frame[, i]
        index <- index + 1
    }
    write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
    close(f)
    options(warn = 0)

  }

cluster.and.resort <- function(V, rows, cols) {
        dist.matrix <- dist(t(V))
        HC <- hclust(dist.matrix, method="complete")
        V <- V[, HC$order]
        cols <- cols[HC$order]
        dist.matrix <- dist(V)
        HC <- hclust(dist.matrix, method="complete")
        V <- V[HC$order, ]
        rows <- rows[HC$order]
        return(list(V = V, rows = rows, cols = cols))
      }

HeatMapPlot <- function(
V, 
row.names = NULL,
col.names = NULL,
main = " ", 
sub = " ", 
xlab=" ", 
ylab=" ",
row.norm = TRUE,
cmap.type = 1)   # 1 = red/bluecologram, 2 = scale of violets
{
       n.rows <- length(V[,1])
       n.cols <- length(V[1,])

       if (row.norm == TRUE) {
          row.mean <- apply(V, MARGIN=1, FUN=mean)
          row.sd <- apply(V, MARGIN=1, FUN=sd)
          row.n <- length(V[,1])
          for (i in 1:n.rows) {
	     if (row.sd[i] == 0) {
    	         V[i,] <- 0
             } else {
	         V[i,] <- (V[i,] - row.mean[i])/(0.333 * row.sd[i])
             }
             V[i,] <- ifelse(V[i,] < -6, -6, V[i,])
             V[i,] <- ifelse(V[i,] > 6, 6, V[i,])
          }
        }

       if (cmap.type == 1) { 
           red.blue.palette <- colorRampPalette(c("red", "white", "blue"), space = "rgb")
           mycol <- rev(red.blue.palette(20))
        } else if (cmap.type == 2) {  # range of violet
           violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
           mycol <- rev(violet.palette(20))
        }
        ncolors <- length(mycol) - 2

        heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
        heatm[1:n.rows,] <- V[seq(n.rows, 1, -1),]
        maxv <- max(V)
        minv <- min(V)
        rangev <- maxv - minv
#        windows(width=14, height=9)
        par(mar = c(6, 12, 3, 3))
        image(1:n.cols, 1:n.rows, t(heatm), col=mycol, axes=FALSE, main=main, sub = sub, xlab= xlab, ylab=ylab)
        if (!is.null(row.names)) {
            numC <- nchar(row.names)
            
            size.row.char <- 35/(ifelse(n.rows > 50, 50, n.rows) + 12)
            for (i in 1:n.rows) {
               row.names[i] <- substr(row.names[i], 1, 35)
            }
            axis(2, at=1:n.rows, labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char, font.axis=2, line=-1)
        }
        if (!is.null(col.names)) {
            numC <- nchar(col.names)
            size.col.char <- 20/(n.cols + 15)
            for (i in 1:n.cols) {
               col.names[i] <- substr(col.names[i], 1, 35)
            }
           axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
        }

# version with sorted row and cols

#        dist.matrix <- dist(t(heatm))
#        HC <- hclust(dist.matrix, method="complete")
#        heatm <- heatm[, HC$order]
#        col.names <- col.names[HC$order]
#        dist.matrix <- dist(heatm)
#        HC <- hclust(dist.matrix, method="complete")
#        heatm <- heatm[HC$order, ]
#        row.names <- row.names[HC$order]
#        windows(width=14, height=9)
#        par(mar = c(6, 12, 3, 3))
#        image(1:n.cols, 1:n.rows, t(heatm), col=mycol, axes=FALSE, main=main, sub = sub, xlab= xlab, ylab=ylab)
#        if (!is.null(row.names)) {
#            axis(2, at=1:n.rows, labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char, font.axis=2, line=-1)
#        }
#        if (!is.null(col.names)) {
#            axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
#        }
       
	return()
}

NMF <- function(V, k, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {
        N <- length(V[,1])
        M <- length(V[1,])
        set.seed(seed)
        W <- matrix(runif(N*k), nrow = N, ncol = k)  # Initialize W and H with random numbers
        H <- matrix(runif(k*M), nrow = k, ncol = M)
        VP <- matrix(nrow = N, ncol = M)
        error.v <- vector(mode = "numeric", length = maxniter)
        new.membership <- vector(mode = "numeric", length = M)
        old.membership <- vector(mode = "numeric", length = M)
        eps <- .Machine$double.eps
        for (t in 1:maxniter) {
              VP = W %*% H
              H <- H * (crossprod(W, V)/crossprod(W, VP)) + eps
              VP = W %*% H
              H.t <- t(H)
              W <- W * (V %*% H.t)/(VP %*% H.t) + eps
              error.v[t] <- sqrt(sum((V - VP)^2))/(N * M)
               if (t %% stopfreq == 0) {
                    for (j in 1:M) {
                        class <- order(H[,j], decreasing=T)
                        new.membership[j] <- class[1]
                     }
                     if (sum(new.membership == old.membership) == M) {
                        no.change.count <- no.change.count + 1
                     } else {
                        no.change.count <- 0
                     }
                     if (no.change.count == stopconv) break
                     old.membership <- new.membership
               }
        }
        return(list(W = W, H = H, t = t, error.v = error.v))
}

nnls.fit <- function(x,y,wsqrt=1,eps=0,rank.tol=1e-07) {
  ## Purpose: Nonnegative Least Squares (similar to the S-Plus function
  ## with the same name) with the help of the R-library quadprog
  ## ------------------------------------------------------------------------
  ## Attention:
  ## - weights are square roots of usual weights
  ## - the constraint is coefficient>=eps
  ## ------------------------------------------------------------------------
  ## Author: Marcel Wolbers, July 99
  ##
  ##========================================================================
  require ("quadprog")
  m <- NCOL(x)
  if (length(eps)==1) eps <- rep(eps,m)
  x <- x * wsqrt
  y <- y * wsqrt
#  sometimes a rescaling of x and y helps (if solve.QP.compact fails otherwise)
  xscale <- apply(abs(x),2,mean)
  yscale <- mean(abs(y))
  x <- t(t(x)/xscale)
  y <- y/yscale
  Rinv <- backsolve(qr.R(qr(x)),diag(m))
  cf <- solve.QP.compact(Dmat=Rinv,dvec=t(x)%*%y,Amat=rbind(rep(1,m)),
                   Aind=rbind(rep(1,m),1:m),bvec=eps*xscale/yscale,
                         factorized=TRUE)$sol
  cf <- cf*yscale/xscale  #scale back
  cf
}
