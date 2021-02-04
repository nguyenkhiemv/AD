gradient <- function(B,W,y,index){
    ## B is M x K matrix of basis evaluations, at the M = \sum_{i=1}^N n_i timepoints.
    ## W is N x K current solution matrix
    ## y is M vector of responses (measurements)
    ## index is M vector, indicating which elements belong to which "subject" ie 1,1,1,2,2,3,3,3 ..., N,N,N
    resid <- y - drop( rowSums(B * W[index, ]))
    rB <- resid * B # elementwise
    lrB <- split(rB, index) # drops the dimensions
    K <- ncol(B)
    G <- sapply(lrB, function(x)colSums(matrix(x,ncol=K)))
    t(G)
    }


### Test
set.seed(1)
B <- matrix(rnorm(50),10,5)
y <- rnorm(10)
index=c(1,1,1,1,2,2,2,4,4,4)
W <- matrix(rnorm(15),4,5)
gradient(B,W,y,index)
gradient.old.lukasz(B,W,y,index)

## Do middle one by hand to check

mid <- c(5,6,7)
resid= drop(y[mid]-B[mid,]%*%W[2,])
t(resid)%*%B[mid,]
