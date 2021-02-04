source("new.fimpute.R")
library("fda")
library("fdapace")
library("fcomplete")

set.seed(0)

rep.res = c()
rep.time = c()
for (i in 1:1){

### DATA GENERATION
n = 200 # subjects
m = 7  # number of basis functions

# Create and plot fourier basis
basis = create.bspline.basis(nbasis=m)
plot(basis)

# Create random low-rank functional data in that basis
Wfull = matrix(rnorm(n*m),n,m)
pp = 3
ss = svd(Wfull, nu = pp, nv = pp)
truess = ss
Wtrue = ss$u %*% diag(ss$d[1:pp]) %*% t(ss$v)

# Sample observations
tpoints = 5

functions = fd(t(Wtrue), basis)

# Generate observations (points at some indexes) with noise
T = matrix(runif(n*tpoints),n,tpoints)
OBS = t(sapply(1:n, function(i) {eval.fd(T[i,], functions[i])})) + rnorm(n*tpoints)*0.2

data = data.frame(id = rep(1:n,tpoints), time = matrix(T,ncol = 1), measurement = matrix(OBS,ncol = 1))

id.var = "id"
time.var = "time"
value.var = "measurement"

#test.mask = 1:nobs %in% sample(1:nobs)[1:floor(nobs*0.1)]
lr = 1
lambdas = seq(0.5,5,0.5)
niter = 1e5
df = 7
tol = 1e-8

start_time <- Sys.time()
res = cv.fimpute(data,
                 value.var,
                 time.var,
                 id.var,
                 pp = pp,
                 lr = lr,
                 niter = niter,
                 df = df,
                 tol = tol,
                 lambdas = lambdas,
                 val.ratio = 0.1)

grid = 0:99/99
functions.on.grid = t(eval.fd(grid,functions))
plot(lambdas, res$loss)
bestLambda = lambdas[which.min(res$loss)]


res$model = fimpute.fit(data,
           value.var,
           time.var,
           id.var,
           pp = pp,
           lr = lr,
           niter = niter,
           df = df,
           tol = tol,
           lambda = bestLambda)
time.fimpute = Sys.time() - start_time

err.fimpute = mean((functions.on.grid - res$model$functions$measurement)**2) / mean((functions.on.grid - 0)**2)
cat(paste("MSE =", err.fimpute))

start_time <- Sys.time()
model = fcomplete::fregression(measurement ~ time | id, data = data,bins = 100,lambda = c(0.5,0.75,1,1.25,1.5),K=pp,d = 7)
err.fcomplete = mean((functions.on.grid - model$fit)**2) / mean((functions.on.grid - 0)**2)
cat(paste("MSE =",err.fcomplete))
time.fcomplete = Sys.time() - start_time

start_time <- Sys.time()
fpca.inputs = MakeFPCAInputs(IDs = data$id, tVec = data$time, yVec = data$measurement, sort = TRUE)
FPCAdense <- FPCA(fpca.inputs$Ly, fpca.inputs$Lt, list(nRegGrid=100,maxK=pp))
fpca.curves = (FPCAdense$xiEst %*% t(FPCAdense$phi)) + FPCAdense$bwMu
err.fpca = mean((fpca.curves - functions.on.grid)**2) / mean((functions.on.grid - 0)**2)
cat(paste("MSE =", err.fpca))
time.fpca = Sys.time() - start_time

rep.res = rbind(rep.res, c(err.fimpute, err.fcomplete, err.fpca))
rep.time = rbind(rep.time, c(time.fimpute, time.fcomplete, time.fpca))
}

pdf(paste0("images/simulation-performance.pdf"),width = 5, height = 5)
colnames(rep.res) = c("grid-free", "grid-based", "fPCA")
boxplot(rep.res)
dev.off()

pdf(paste0("images/simulation-time.pdf"),width = 5, height = 5)
colnames(rep.time) = c("grid-free", "grid-based", "fPCA")
boxplot(rep.time)
dev.off()

pdf("images/simulated-curves.pdf",width = 4, height = 4)
plot(functions[1:10],ylim=c(-1,1)) # Note that these functions are rank 2
title("10 sampled 'true' curves")
dev.off()

# Plot one curve and observed points
pdf("images/simulated-1curve.pdf",width = 4, height = 4)
plot(functions[1],ylim=c(-1,1))
points(T[1,],OBS[1,], col="grey50")
title("Subject 1 with observed points")
dev.off()


par(mfrow=c(1,1))
plot.func = function(id){
  pdf(paste0("images/subj-",id,".pdf"),width = 5, height = 5)
  plot(functions[id],ylim=c(-0.75,0.75),lwd=2)
  lines(grid, model$fit[id,],col=2,lwd=2)
  lines(grid, res$model$functions$measurement[id,],col=3,lwd=2)
  lines(grid, fpca.curves[id,],col=4,lwd=2)
  points(T[id,],OBS[id,],lwd=2)
  title(paste("Predictions subj.",id))
  legend(0.5, -0.1,
         c("true curve", "grid-based", "grid-free", "fpca", "observations"),
         lwd = c(1, rep(2,4)),
         pch = c(rep(NA,4),1),
         lty = c(rep(1,4), NA),
         col = c(1:4,1))
  dev.off()
}
plot.func(1)
plot.func(2)
plot.func(3)
plot.func(4)
plot.func(5)
plot.func(6)
plot.func(7)
plot.func(8)
plot.func(10)
plot.func(15)
plot.func(20)

plot(FPCAdense)

plot(model$mean_curves$measurement)
#----------------------


preds.fimpute = predict.fimpute(res$model$functions[[1]], grid, data[!test.mask,], time.var, id.var)
preds.fcomplete = predict.fimpute(model$fit, grid, data[!test.mask,], time.var, id.var)
preds.fpca = predict.fimpute(fpca.curves, grid, data[!test.mask,], time.var, id.var)

plot(data[!test.mask,value.var],preds.fimpute)
plot(data[!test.mask,value.var],preds.fcomplete)
plot(data[!test.mask,value.var],preds.fpca)
