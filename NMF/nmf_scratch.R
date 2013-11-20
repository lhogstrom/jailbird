
   # Perform NMF clustering

   NMF.out <- NMF.div(V = m.2, k = k.comp, maxniter = 2000, seed = 123, stopconv = 40, stopfreq = 10)
   W <- NMF.out$W
   H_orig <- NMF.out$H
   colnames(W) <- paste("c", seq(1, k.comp), sep="")
   row.names(W) <- row.names(m.2)

V = m.2
k = k.comp
maxniter = 2000
seed = 123
stopconv = 40
stopfreq = 10

N = 978
M = 78
k = 9

N = 100
M = 40
k = 9
#test matrix for failing point
g1 <- matrix(runif(N*M), nrow = N, ncol = M)
h1 <- matrix(runif(M*k), nrow = M, ncol = k)
# r = g1 %*% h1
r <- g1 %*% h1

S.train = g1
W.train = h1
#      Z <- S.train %*% W.train
# if the above matrix multiply runs out of memory try this: 



#NMF function
function(V, k, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {

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
                #W <- W * ((V/VP) %*% H.t) + eps
                # W <- W * (((V/VP)+ eps) %*% H.t) + eps #crasshing here
                Vd <- ((V/VP)+ eps)
                Dp <- matrix(0, nrow=dim(Vd)[1], ncol=dim(H.t)[2]) 
                for (i in 1:dim(Vd)[1]) { 
                    for (j in 1:dim(H.t)[2]) { 
                      Dp[i, j] <- Vd[i,] %*% H.t[,j] 
                    } 
                } 
                W <- W * Dp + eps
                norm <- apply(H, MARGIN=1, FUN=sum)
                for (i in 1:k) {
                    W[,i] <- W[,i]/norm[i]
                }
               error.v[t] <- sum(V * log((V + eps)/(VP + eps)) - V + VP)/(M * N) # cannot be Z for denom
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



## python code
# import cmap.io.gct as gct
# gt = gct.GCT()
# gctF = '/xchip/cogs/hogstrom/analysis/pablos_NMF_analysis/TA/CC_lh/MCF7_top_intra_connecting_compound_classes_n79x978.NORM.gct'        
# gt.read(gctF)


#      Z <- S.train %*% W.train
# if the above matrix multiply runs out of memory try this: 
     X <- matrix(0, nrow=dim(S.train)[1], ncol=dim(W.train)[2]) 
     for (i in 1:dim(S.train)[1]) { 
       for (j in 1:dim(W.train)[2]) { 
          X[i, j] <- S.train[i,] %*% W.train[,j] 
       } 
     } 

# loop - Pablo's scratch code
for (i in 1:N) {
    for (j in 1:M) {
            s = -
            for (k in 1:K.comp)
            s = s+A[i,k] * B[k,j]
    }
}


#matrix multiplication w/ scalapack
