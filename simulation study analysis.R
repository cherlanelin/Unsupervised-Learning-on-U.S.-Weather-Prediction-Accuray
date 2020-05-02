## This is the R file including all the content related to simulation study.
## Please run all of it from the beginning to the end by order
## The main content of this file includes:
##     0. Installing required packages 
##     1. Functions for Data simulation and Analysis
##     2. Smooth FPCA on B-spline Non-parametric Regression
##     3. Clustering Number Selection Validation
##     4. Clustering Validation Study
######################################################################################



######### Installing required packages #####################################
# ipak, a function for checking the the installation of the packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# add new packages here
packages <- c("ggplot2", "dplyr", "tidyverse", "tidyr", "geosphere",
              "ggridges","scales","geofacet","lubridate","fda","RColorBrewer",
              "weathermetrics","ggrepel","DT","forecast","tseries","maps",
              "ggthemes","gridExtra","factoextra", "randomForest", 
              "NbClust","funFEM","fda.usc","tclust","dtwclust")
ipak(packages)



######## Functions for Data simulation and Analysis ##########################
### Function for calculating the CV error for a given lambda 
get.mcv.smoothfpca = function(lambda, data, ngroup, n_harm){
  range = c(as.numeric(min(data$t)),as.numeric(max(data$t)))
  bks.1 = as.Date(quantile(unclass(data$t), seq(0,1,length = 15)), origin = "1970-01-01")
  err1.basis = create.bspline.basis(rangeval = c(min(data$t),max(data$t)), breaks = bks.1, norder = 4)
  mindiff.state.mat.1 = matrix(unlist(data[,-which(names(data) %in% c("t"))]),ncol = 4*ngroup)
  err1.fd<-smooth.basis(argvals = data$t, y = mindiff.state.mat.1, err1.basis)$fd
  bspFpcPen = fdPar(fdobj=err1.basis, Lfdobj=2, lambda = lambda)
  smooth.fpc.cv = c()
  for (i in 1:50){
    smooth.fpc.raw = err1.fd
    smooth.fpc.raw$coefs = err1.fd$coefs[,-i]
    smooth.fpc = pca.fd(smooth.fpc.raw,nharm=n_harm,harmfdPar=bspFpcPen)
    smooth.fpc.coefs = as.matrix(smooth.fpc$harmonics$coefs) #17*5
    smooth.fpc.score = as.matrix(smooth.fpc$scores) #49*5
    
    smooth.fpc.xi.fd = err1.fd
    smooth.fpc.xi.fd$coefs = err1.fd$coefs[,i]
    smooth.fpc.xi = eval.fd(data$t,smooth.fpc.xi.fd) - eval.fd(data$t,smooth.fpc$meanfd) #1033*1
    smooth.fpc.basismat = eval.fd(data$t,smooth.fpc$harmonics) 
    
    # estimate the SSE = y'y - y'X(X'X)X'y
    smooth.fpc.ri = t(smooth.fpc.xi)%*%smooth.fpc.xi - t(smooth.fpc.xi)%*%
      smooth.fpc.basismat%*%ginv(t(smooth.fpc.basismat)%*%smooth.fpc.basismat)%*%
      t(smooth.fpc.basismat)%*%smooth.fpc.xi
    smooth.fpc.cv = c(smooth.fpc.cv,smooth.fpc.ri)
  }
  return(sum(smooth.fpc.cv))
}

### Data Simulation Function with Common Noise
generate.curve<-function(t.min,t.max,t.by,n.group,err.sd){
  t = seq(t.min, t.max, by = t.by)
  
  # Assuming there are 4 clusters in the real situation
  # The following four equation is the mean functions respect to each cluster
  mean.clu1 = sin(2*t) 
  mean.clu2 = 2*sin(2*t) 
  mean.clu3 = sin(t)/2
  mean.clu4 = sin(4*t)
  
  simu.clu1 = data.frame(matrix(ncol = n.group, nrow = length(t)))
  name.clu1 = paste0("clu1.",1:n.group)
  colnames(simu.clu1) <- name.clu1
  
  for(column in name.clu1){
    simu.clu1[,column] = mean.clu1 + rnorm(length(t),mean = 0, sd = err.sd)
  }
  
  simu.clu2 = data.frame(matrix(ncol = n.group, nrow = length(t)))
  name.clu2 = paste0("clu2.",1:n.group)
  colnames(simu.clu2) <- name.clu2
  
  for(column in name.clu2){
    simu.clu2[,column] = mean.clu2 + rnorm(length(t),mean = 0, sd = err.sd)
  }
  
  simu.clu3 = data.frame(matrix(ncol = n.group, nrow = length(t)))
  name.clu3 = paste0("clu3.",1:n.group)
  colnames(simu.clu3) <- name.clu3
  
  for(column in name.clu3){
    simu.clu3[,column] = mean.clu3 + rnorm(length(t),mean = 0, sd = err.sd)
  }
  
  simu.clu4 = data.frame(matrix(ncol = n.group, nrow = length(t)))
  name.clu4 = paste0("clu4.",1:n.group)
  colnames(simu.clu4) <- name.clu4
  
  for(column in name.clu4){
    simu.clu4[,column] = mean.clu4 + rnorm(length(t),mean = 0, sd = err.sd)
  }
  
  simu.all = cbind(simu.clu1,simu.clu2,simu.clu3,simu.clu4)
  simu.all$t = t
  
  return (simu.all)
}


### Data Simulation Function with Absolute Value-Proportion Noise
generate.curve.ctvardiff<-function(t.min,t.max,t.by,n.group,err.sd){
  t = seq(t.min, t.max, by = t.by)
  
  
  mean.clu1 = sin(2*t) 
  mean.clu2 = 2*sin(2*t) 
  mean.clu3 = sin(t)/2
  mean.clu4 = sin(4*t)
  
  simu.clu1 = data.frame(matrix(ncol = n.group, nrow = length(t)))
  name.clu1 = paste0("clu1.",1:n.group)
  colnames(simu.clu1) <- name.clu1
  
  for(column in name.clu1){
    simu.clu1[,column] = mean.clu1 + rnorm(length(t),mean = 0, sd = err.sd)*abs(sin(2*t))
  }
  
  simu.clu2 = data.frame(matrix(ncol = n.group, nrow = length(t)))
  name.clu2 = paste0("clu2.",1:n.group)
  colnames(simu.clu2) <- name.clu2
  
  for(column in name.clu2){
    simu.clu2[,column] = mean.clu2 + rnorm(length(t),mean = 0, sd = err.sd)*abs(2*sin(2*t))
  }
  
  simu.clu3 = data.frame(matrix(ncol = n.group, nrow = length(t)))
  name.clu3 = paste0("clu3.",1:n.group)
  colnames(simu.clu3) <- name.clu3
  
  for(column in name.clu3){
    simu.clu3[,column] = mean.clu3 + rnorm(length(t),mean = 0, sd = err.sd)*abs(sin(t)/2)
  }
  
  simu.clu4 = data.frame(matrix(ncol = n.group, nrow = length(t)))
  name.clu4 = paste0("clu4.",1:n.group)
  colnames(simu.clu4) <- name.clu4
  
  for(column in name.clu4){
    simu.clu4[,column] = mean.clu4 + rnorm(length(t),mean = 0, sd = err.sd)*abs(sin(4*t))
  }
  
  simu.all = cbind(simu.clu1,simu.clu2,simu.clu3,simu.clu4)
  simu.all$t = t
  return (simu.all)
}


### Functions for Cluster Number Selection Study
cluster.select.num = function(simu.all,n.group,knots.num,nharm,plotting = FALSE){
  t = simu.all$t
  simu.all = simu.all[,-which(names(simu.all) %in% c("t"))]
  simu.all.mat = matrix(unlist(simu.all),ncol = 4*n.group)
  
  bks = quantile(t, seq(0,1,length = knots.num+1))
  simu.basis = create.bspline.basis(rangeval = c(min(t),max(t)), breaks = bks, norder = 4) 
  simu.fd<-smooth.basis(argvals = t, y = simu.all.mat, simu.basis)$fd
  
  if(plotting == TRUE){plot(simu.fd)}
  
  ## FunFEM
  fem.clu.bic = funFEM(simu.fd,K = 2:6,model = "all",crit = "bic", init = "kmeans",eps = 1e-5)
  fem.num.bic = fem.clu.bic$K
  
  fem.clu.icl = funFEM(simu.fd,K = 2:6,model = "all",crit = "icl", init = "kmeans",eps = 1e-5)
  fem.num.icl = fem.clu.icl$K
  
  
  ## kmeans on coef
  fd.coefs = t(simu.fd$coefs)
  # all index  
  coeff.kmeans.best <- NbClust(fd.coefs, distance = "euclidean", 
                               min.nc=2, max.nc=6, 
                               method = "kmeans", 
                               index = "all")
  bsp.clu.select =  max(coeff.kmeans.best$Best.partition)
  
  
  ## kmeans on smoothed FPC scores
  range = c(as.numeric(min(t)),as.numeric(max(t)))
  bspFpcPen = fdPar(fdobj=simu.basis, Lfdobj=2, lambda = 10^(-2))
  smooth.fpc = pca.fd(simu.fd,nharm=2,harmfdPar=bspFpcPen)
  fd.scores = smooth.fpc$scores
  # all index  
  fpc.kmeans.best <- NbClust(fd.scores, distance = "euclidean", 
                             min.nc=2, max.nc=6, 
                             method = "kmeans", 
                             index = "all")
  fpc.clu.select =  max(fpc.kmeans.best$Best.partition)
  
  
  simu.clunum.compare = rbind(fem.num.bic,fem.num.icl,
                              bsp.clu.select,fpc.clu.select)
  return(simu.clunum.compare)
}


### Functions for Validation
# Hungarian Algorithm 
minWeightBipartiteMatching <- function(clusteringA, clusteringB) {
  require(clue)
  idsA <- unique(clusteringA)  # distinct cluster ids in a
  idsB <- unique(clusteringB)  # distinct cluster ids in b
  nA <- length(clusteringA)  # number of instances in a
  nB <- length(clusteringB)  # number of instances in b
  if (length(idsA) != length(idsB) || nA != nB) {
    stop("number of cluster or number of instances do not match")
  }
  
  nC <- length(idsA)
  tupel <- c(1:nA)
  
  # computing the distance matrix
  assignmentMatrix <- matrix(rep(-1, nC * nC), nrow = nC)
  for (i in 1:nC) {
    tupelClusterI <- tupel[clusteringA == i]
    solRowI <- sapply(1:nC, function(i, clusterIDsB, tupelA_I) {
      nA_I <- length(tupelA_I)  # number of elements in cluster I
      tupelB_I <- tupel[clusterIDsB == i]
      nB_I <- length(tupelB_I)
      nTupelIntersect <- length(intersect(tupelA_I, tupelB_I))
      return((nA_I - nTupelIntersect) + (nB_I - nTupelIntersect))
    }, clusteringB, tupelClusterI)
    assignmentMatrix[i, ] <- solRowI
  }
  
  # optimization
  result <- solve_LSAP(assignmentMatrix, maximum = FALSE)
  attr(result, "assignmentMatrix") <- assignmentMatrix
  return(result)
}


# Function to get cluster result
get.df.cluster = function(simu.all,n.group,knots.num,nharm){
  t = simu.all$t
  simu.all = simu.all[,-which(names(simu.all) %in% c("t"))]
  simu.all.mat = matrix(unlist(simu.all),ncol = 4*n.group)
  
  bks = quantile(t, seq(0,1,length = knots.num+1))
  simu.basis = create.bspline.basis(rangeval = c(min(t),max(t)), breaks = bks, norder = 4) 
  simu.fd<-smooth.basis(argvals = t, y = simu.all.mat, simu.basis)$fd
  
  # fem on icl
  fem.icl <- funFEM(simu.fd,K = 4,model = "all",crit = "icl", init = "kmeans",eps = 1e-5) # 4 cluster solution
  df.fem.icl <- fem.icl$cls
  
  # fem on bic
  fem.bic <- funFEM(simu.fd,K = 4,model = "all",crit = "bic", init = "kmeans",eps = 1e-5) # 4 cluster solution
  df.fem.bic <- fem.bic$cls
  
  # kmeans on bsp coeff
  simu.fd.coef = t(simu.fd$coefs)
  coeff.k <- kmeans(simu.fd.coef, 4) # 4 cluster solution
  df.bsp.k <- coeff.k$cluster
  
  # Kmeans on smoothed pc scores
  range = c(as.numeric(min(t)),as.numeric(max(t)))
  bspFpcPen = fdPar(fdobj=simu.basis, Lfdobj=2, lambda = 10^(-2))
  ## We will explain this lambda selection in the next part
  smooth.fpc = pca.fd(simu.fd,nharm=nharm,harmfdPar=bspFpcPen)
  df.pcscore = smooth.fpc$scores
  fpc.k <- kmeans(df.pcscore, 4) # 4 cluster solution
  df.fpc.k <- fpc.k$cluster
  
  # Hungarian Alg
  real = rep(c(1,2,3,4),each = n.group)
  a = minWeightBipartiteMatching(df.bsp.k,real)
  df.bsp.k.1 = a[df.bsp.k]
  
  b = minWeightBipartiteMatching(df.fpc.k,real)
  df.fpc.k.1 = b[df.fpc.k]
  
  c = minWeightBipartiteMatching(df.fem.bic,real)
  df.fem.bic.1 = c[df.fem.bic]
  
  d = minWeightBipartiteMatching(df.fem.icl,real)
  df.fem.icl.1 = d[df.fem.icl]
  
  clu.result = data.frame(clu.real = real,
                          clu.bsp = df.bsp.k.1,
                          clu.fpc = df.fpc.k.1,
                          clu.fem.bic = df.fem.bic.1,
                          clu.fem.icl = df.fem.icl.1)
  return(clu.result)
}


# Accuracy test
get.clu.accuracy = function(df.clu, n.group){
  bsp.accu = sum(df.clu$clu.real == df.clu$clu.bsp)/(n.group*4) 
  fpc.accu = sum(df.clu$clu.real == df.clu$clu.fpc)/(n.group*4)
  fem.bic.accu = sum(df.clu$clu.real == df.clu$clu.fem.bic)/(n.group*4)
  fem.icl.accu = sum(df.clu$clu.real == df.clu$clu.fem.icl)/(n.group*4)
  return(c(bsp.accu,fpc.accu,fem.bic.accu,fem.icl.accu))
}


# Function on Average Within Cluster SD and Distance to Real Mean of Cluster
get.clu.value = function(df, df.clu, n.group, method){
  
  t = df$t
  bks = quantile(t, seq(0,1,length = 20+1))
  fd.basis = create.bspline.basis(rangeval = c(min(t),max(t)), breaks = bks, norder = 4)
  
  # real clusters
  real.clu1 = sin(2*t) 
  real.clu1.fd = smooth.basis(argvals = t, y = real.clu1, fd.basis)$fd
  
  real.clu2 = 2*sin(2*t) 
  real.clu2.fd = smooth.basis(argvals = t, y = real.clu2, fd.basis)$fd
  
  real.clu3 = sin(t)/2 
  real.clu3.fd = smooth.basis(argvals = t, y = real.clu3, fd.basis)$fd
  
  real.clu4 = sin(4*t) 
  real.clu4.fd = smooth.basis(argvals = t, y = real.clu4, fd.basis)$fd
  
  temp.clu = df.clu
  if (method == "fpc"){
    temp.clu$clu = df.clu$clu.fpc
  }else if (method == "bsp"){
    temp.clu$clu = df.clu$clu.bsp
  }else if (method == "fem.bic"){
    temp.clu$clu = df.clu$clu.fem.bic
  }else if (method == "fem.icl"){
    temp.clu$clu = df.clu$clu.fem.icl
  }
  
  ## est.clu1 
  temp.clu$curve = colnames(df[,-which(names(df) %in% c("t"))])
  clu1 = temp.clu[temp.clu$clu == 1,]$curve
  clu1.ts = matrix(unlist(df[,names(df)%in%clu1]),ncol = length(clu1))
  clu1.fd = smooth.basis(argvals = t, y = clu1.ts, fd.basis)$fd
  clu1.mu = mean.fd(clu1.fd)
  clu1.sd = int.simpson(sd.fd(clu1.fd),method="TRAPZ")
  clu1.mu.diff = mean((real.clu1 - eval.fd(t,clu1.mu))^2)
  
  ## est.clu2 
  temp.clu$curve = colnames(df[,-which(names(df) %in% c("t"))])
  clu2 = temp.clu[temp.clu$clu == 2,]$curve
  clu2.ts = matrix(unlist(df[,names(df)%in%clu2]),ncol = length(clu2))
  clu2.fd = smooth.basis(argvals = t, y = clu2.ts, fd.basis)$fd
  clu2.mu = mean.fd(clu2.fd)
  clu2.sd = int.simpson(sd.fd(clu2.fd),method="TRAPZ")
  clu2.mu.diff = mean((real.clu2 - eval.fd(t,clu2.mu))^2)
  
  ## est.clu3 
  temp.clu$curve = colnames(df[,-which(names(df) %in% c("t"))])
  clu3 = temp.clu[temp.clu$clu == 3,]$curve
  clu3.ts = matrix(unlist(df[,names(df)%in%clu3]),ncol = length(clu3))
  clu3.fd = smooth.basis(argvals = t, y = clu3.ts, fd.basis)$fd
  clu3.mu = mean.fd(clu3.fd)
  clu3.sd = int.simpson(sd.fd(clu3.fd),method="TRAPZ")
  clu3.mu.diff = mean((real.clu3 - eval.fd(t,clu3.mu))^2)
  
  ## est.clu4 
  temp.clu$curve = colnames(df[,-which(names(df) %in% c("t"))])
  clu4 = temp.clu[temp.clu$clu == 4,]$curve
  clu4.ts = matrix(unlist(df[,names(df)%in%clu4]),ncol = length(clu4))
  clu4.fd = smooth.basis(argvals = t, y = clu4.ts, fd.basis)$fd
  clu4.mu = mean.fd(clu4.fd)
  clu4.sd = int.simpson(sd.fd(clu4.fd),method="TRAPZ")
  clu4.mu.diff = mean((real.clu4 - eval.fd(t,clu4.mu))^2)
  
  return(c(mean(clu1.sd,clu2.sd,clu3.sd,clu4.sd),
           mean(clu1.mu.diff,clu2.mu.diff,clu3.mu.diff,clu4.mu.diff)))
}

######## B-spline Basis Selection and Smooth FPCA Lambda Selection ##########################
# We will discuss 2 scenarios: common noise V.S. absolute value-proportion noise
# Under each scenario, we consider:
# 1. The overall varianace of the noise: sigma= 1 V.S. 2
#  To conclude, we have 4 run of simulation 
# (noise pattern X noise overall variance X number of cluster)

## Determine number of PC functions and smoothing parameter lambda in each situation
# Round 1: n = 20, Sigma = 1, common noise
set.seed(888)
df = generate.curve(0,10,0.01,20,1)
t = df$t
df = df[,-which(names(df) %in% c("t"))]
df.mat = matrix(unlist(df),ncol = 4*20)
bks = quantile(t, seq(0,1,length = 20+1))
df.basis = create.bspline.basis(rangeval = c(min(t),max(t)), breaks = bks, norder = 4)
df.fd<-smooth.basis(argvals = t, y = df.mat, df.basis)$fd
raw.fpca = pca.fd(df.fd,nharm=2) ## Select 2: reach 94.83% at PC2 (76.811% at PC1)
sum(raw.fpca$varprop)

df$t = t
table.lambda.mcv<-c()
for (j in -40:0) { 
  smooth.fpc.cv  =  get.mcv.smoothfpca(10^(j/5),df,20,2)
  table.lambda.mcv <- rbind(table.lambda.mcv,c(j,smooth.fpc.cv))
}
plot(table.lambda.mcv[,1]/5, table.lambda.mcv[,2], type = "l", 
     main = "CV plot Lambda Selection Scenario 1 (n = 20, Sigma = 1)", 
     xlab = "log(Lambda) Value", ylab = "CV")

# Round 2: n = 20, Sigma = 2, common noise
set.seed(888)
df = generate.curve(0,10,0.01,20,2)
t = df$t
df = df[,-which(names(df) %in% c("t"))]
df.mat = matrix(unlist(df),ncol = 4*20)
bks = quantile(t, seq(0,1,length = 20+1))
df.basis = create.bspline.basis(rangeval = c(min(t),max(t)), breaks = bks, norder = 4)
df.fd<-smooth.basis(argvals = t, y = df.mat, df.basis)$fd
raw.fpca = pca.fd(df.fd,nharm=6) ## Select 6: reach 90.14% at PC6 (88.98% at PC5)
sum(raw.fpca$varprop) 

df$t = t
table.lambda.mcv<-c()
for (j in -40:0) { 
  smooth.fpc.cv  =  get.mcv.smoothfpca(10^(j/5),df,20,6)
  table.lambda.mcv <- rbind(table.lambda.mcv,c(j,smooth.fpc.cv))
}
plot(table.lambda.mcv[,1]/5, table.lambda.mcv[,2], type = "l", 
     main = "CV plot Lambda Selection Scenario 1 (n = 20, Sigma = 2)", 
     xlab = "log(Lambda) Value", ylab = "CV")

# Round 3: n = 20, Sigma = 1, absolute-value proportion noise
set.seed(888)
df = generate.curve.ctvardiff(0,10,0.01,20,1)
t = df$t
df = df[,-which(names(df) %in% c("t"))]
df.mat = matrix(unlist(df),ncol = 4*20)
bks = quantile(t, seq(0,1,length = 20+1))
df.basis = create.bspline.basis(rangeval = c(min(t),max(t)), breaks = bks, norder = 4)
df.fd<-smooth.basis(argvals = t, y = df.mat, df.basis)$fd
raw.fpca = pca.fd(df.fd,nharm=2) ## Select 2: reach 96.01% at PC2 (78% at PC1)
sum(raw.fpca$varprop)

df$t = t
table.lambda.mcv<-c()
for (j in -40:0){ 
  smooth.fpc.cv  =  get.mcv.smoothfpca(10^(j/5),df,20,2)
  table.lambda.mcv <- rbind(table.lambda.mcv,c(j,smooth.fpc.cv))
}
plot(table.lambda.mcv[,1]/5, table.lambda.mcv[,2], type = "l", 
     main = "CV plot Lambda Selection Scenario 2 (n = 20, Sigma = 1)", 
     xlab = "Lambda Value", ylab = "CV")

# Round 4: n = 20, Sigma = 2, absolute-value proportion noise
set.seed(888)
df = generate.curve.ctvardiff(0,10,0.01,20,2)
t = df$t
df = df[,-which(names(df) %in% c("t"))]
df.mat = matrix(unlist(df),ncol = 4*20)
bks = quantile(t, seq(0,1,length = 20+1))
df.basis = create.bspline.basis(rangeval = c(min(t),max(t)), breaks = bks, norder = 4)
df.fd<-smooth.basis(argvals = t, y = df.mat, df.basis)$fd
raw.fpca = pca.fd(df.fd,nharm=4) ## Select 4: reach 91.40% at PC4 (89.9% at PC3)
sum(raw.fpca$varprop) 

df$t = t
table.lambda.mcv<-c()
for (j in -40:0) { 
  smooth.fpc.cv  =  get.mcv.smoothfpca(10^(j/5),df,20,4)
  table.lambda.mcv <- rbind(table.lambda.mcv,c(j,smooth.fpc.cv))
}
plot(table.lambda.mcv[,1]/5, table.lambda.mcv[,2], type = "l", 
     main = "CV plot Lambda Selection Scenario 2 (n = 20, Sigma = 2)", 
     xlab = "Lambda Value", ylab = "CV")


#### From the plots, we find that the lambda is around 10^(-2)


######## Clustering Number Selection Validation ##########################
# Round 1: n = 20, Sigma = 1, common noise
simunum.round1.20= data.frame(matrix(ncol = 5, nrow = 1))
colnames(simunum.round1.20) <- c("method","index","clu.num","n.group","round")
seed = 888
for (i in  1:200){
  set.seed(seed)
  df = generate.curve(t.min=0,t.max=10,t.by=0.01,
                      n.group=20,err.sd=1)
  cluster.result = data.frame(matrix(ncol = 5, nrow = 4))
  colnames(cluster.result) <- c("method","index","clu.num","n.group","round")
  cluster.result$method = c("fem","fem","bsp","fpc")
  cluster.result$index = c("bic","icl","all","all")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 1, common noise",each = 4)
  cluster.result$clu.num = cluster.select.num(simu.all=df,n.group=20,knots.num=20,nharm = 2)
  
  simunum.round1.20 = rbind(simunum.round1.20,cluster.result)
  seed = seed+10000*i
}
simunum.round1.20 = na.omit(simunum.round1.20)
write.csv(simunum.round1.20, file = "simunum_round1_20.csv", row.names = FALSE)

##
# Round 2: n = 20, Sigma = 2, common noise
simunum.round2.20= data.frame(matrix(ncol = 5, nrow = 1))
colnames(simunum.round2.20) <- c("method","index","clu.num","n.group","round")
seed = 888
for (i in  1:200){
  set.seed(seed)
  df = generate.curve(t.min=0,t.max=10,t.by=0.01,
                      n.group=20,err.sd=2)
  cluster.result = data.frame(matrix(ncol = 5, nrow = 4))
  colnames(cluster.result) <- c("method","index","clu.num","n.group","round")
  cluster.result$method = c("fem","fem","bsp","fpc")
  cluster.result$index = c("bic","icl","all","all")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 2, common noise",each = 4)
  cluster.result$clu.num = cluster.select.num(simu.all=df,n.group=20,knots.num=20,nharm = 6)
  
  simunum.round2.20 = rbind(simunum.round2.20,cluster.result)
  seed = seed+10000*i
}
simunum.round2.20 = na.omit(simunum.round2.20)
write.csv(simunum.round2.20, file = "simunum_round2_20.csv", row.names = FALSE)


# Round 3: n = 20, Sigma = 1, absolute-value proportion noise
simunum.round3.20= data.frame(matrix(ncol = 5, nrow = 1))
colnames(simunum.round3.20) <- c("method","index","clu.num","n.group","round")
seed = 888
for (i in  1:200){
  set.seed(seed)
  df = generate.curve.ctvardiff(t.min=0,t.max=10,t.by=0.01,
                                n.group=20,err.sd=1)
  cluster.result = data.frame(matrix(ncol = 5, nrow = 4))
  colnames(cluster.result) <- c("method","index","clu.num","n.group","round")
  cluster.result$method = c("fem","fem","bsp","fpc")
  cluster.result$index = c("bic","icl","all","all")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 1, varies noise",each = 4)
  cluster.result$clu.num = cluster.select.num(simu.all=df,n.group=20,knots.num=20,nharm = 2)
  
  simunum.round3.20 = rbind(simunum.round3.20,cluster.result)
  seed = seed+10000*i
}
simunum.round3.20 = na.omit(simunum.round3.20)
write.csv(simunum.round3.20, file = "simunum_round3_20.csv", row.names = FALSE)

# Round 4: n = 20, Sigma = 2, absolute-value proportion noise
simunum.round4.20= data.frame(matrix(ncol = 5, nrow = 1))
colnames(simunum.round4.20) <- c("method","index","clu.num","n.group","round")
seed = 888
for (i in  1:200){
  set.seed(seed)
  df = generate.curve.ctvardiff(t.min=0,t.max=10,t.by=0.01,
                                n.group=20,err.sd=2)
  cluster.result = data.frame(matrix(ncol = 5, nrow = 4))
  colnames(cluster.result) <- c("method","index","clu.num","n.group","round")
  cluster.result$method = c("fem","fem","bsp","fpc")
  cluster.result$index = c("bic","icl","all","all")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 2, varies noise",each = 4)
  cluster.result$clu.num = cluster.select.num(simu.all=df,n.group=20,knots.num=20,nharm = 4)
  
  simunum.round4.20 = rbind(simunum.round4.20,cluster.result)
  seed = seed+10000*i
}
simunum.round4.20 = na.omit(simunum.round4.20)
write.csv(simunum.round4.20, file = "simunum_round4_20.csv", row.names = FALSE)

### Summary analysis of the cluster number detection.
r1.20.raw = read.csv("simunum_round1_20.csv",header = TRUE)

r2.20.raw = read.csv("simunum_round2_20.csv",header = TRUE)

r3.20.raw = read.csv("simunum_round3_20.csv",header = TRUE)

r4.20.raw = read.csv("simunum_round4_20.csv",header = TRUE)

#round 1
r1.20 = r1.20.raw%>%group_by(round, n.group, method,
                             index, clu.num)%>%summarize(count = n())
r1.20

#round 2
r2.20 = r2.20.raw%>%group_by(round, n.group, method,
                             index, clu.num)%>%summarize(count = n())
r2.20

#round 3
r3.20 = r3.20.raw%>%group_by(round, n.group, method,
                             index, clu.num)%>%summarize(count = n())
r3.20

#round 4
r4.20 = r4.20.raw%>%group_by(round, n.group, method,
                             index, clu.num)%>%summarize(count = n())
r4.20

##
########### Clustering Validation Study #####################################
# Round 1: n = 20, Sigma = 1, common noise
simuvid.round1.20= data.frame(matrix(ncol = 7, nrow = 1))
colnames(simuvid.round1.20) <- c("method","index","n.group","round",
                                 "accuracy","avg.sd","avg.mudis")
seed = 888
for (i in  1:200){
  set.seed(seed)
  df = generate.curve(t.min=0,t.max=10,t.by=0.01,
                      n.group=20,err.sd=1)
  cluster.result = data.frame(matrix(ncol = 7, nrow = 4))
  colnames(cluster.result) <- c("method","index","n.group","round",
                                "accuracy","avg.sd","avg.mudis")
  cluster.result$method = c("bsp","fpc","fem","fem")
  cluster.result$index = c("all","all","bic","icl")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 1, common noise",each = 4)
  df.clu = get.df.cluster(df,n.group = 20, knots.num = 20, nharm = 2)
  
  cluster.result$accuracy = get.clu.accuracy(df.clu, 20)
  
  df.bsp = get.clu.value(df, df.clu, 20, "bsp")
  df.fpc = get.clu.value(df, df.clu, 20, "fpc")
  df.bic = get.clu.value(df, df.clu, 20, "fem.bic")
  df.icl = get.clu.value(df, df.clu, 20, "fem.icl")
  cluster.result$avg.sd = c(df.bsp[1],df.fpc[1],df.bic[1],df.icl[1])
  cluster.result$avg.mudis = c(df.bsp[2],df.fpc[2],df.bic[2],df.icl[2])
  
  simuvid.round1.20 = rbind(simuvid.round1.20,cluster.result)
  seed = seed+10000*i
}
simuvid.round1.20 = na.omit(simuvid.round1.20)
write.csv(simuvid.round1.20, file = "simuvid_round1_20.csv", row.names = FALSE)

# Round 2: n = 20, Sigma = 2, common noise
simuvid.round2.20= data.frame(matrix(ncol = 7, nrow = 1))
colnames(simuvid.round2.20) <- c("method","index","n.group","round",
                                 "accuracy","avg.sd","avg.mudis")
seed = 888
for (i in  1:100){
  set.seed(seed)
  df = generate.curve(t.min=0,t.max=10,t.by=0.01,
                      n.group=20,err.sd=2)
  cluster.result = data.frame(matrix(ncol = 7, nrow = 4))
  colnames(cluster.result) <- c("method","index","n.group","round",
                                "accuracy","avg.sd","avg.mudis")
  cluster.result$method = c("bsp","fpc","fem","fem")
  cluster.result$index = c("all","all","bic","icl")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 2, common noise",each = 4)
  df.clu = get.df.cluster(df,n.group = 20, knots.num = 20, nharm = 6)
  
  cluster.result$accuracy = get.clu.accuracy(df.clu, 20)
  
  df.bsp = get.clu.value(df, df.clu, 20, "bsp")
  df.fpc = get.clu.value(df, df.clu, 20, "fpc")
  df.bic = get.clu.value(df, df.clu, 20, "fem.bic")
  df.icl = get.clu.value(df, df.clu, 20, "fem.icl")
  cluster.result$avg.sd = c(df.bsp[1],df.fpc[1],df.bic[1],df.icl[1])
  cluster.result$avg.mudis = c(df.bsp[2],df.fpc[2],df.bic[2],df.icl[2])
  
  simuvid.round2.20 = rbind(simuvid.round2.20,cluster.result)
  seed = seed+10000*i
}

simuvid.round2.20.1 = simuvid.round2.20
seed = 8888
for (i in  1:90){
  set.seed(seed)
  df = generate.curve(t.min=0,t.max=10,t.by=0.01,
                      n.group=20,err.sd=2)
  cluster.result = data.frame(matrix(ncol = 7, nrow = 4))
  colnames(cluster.result) <- c("method","index","n.group","round",
                                "accuracy","avg.sd","avg.mudis")
  cluster.result$method = c("bsp","fpc","fem","fem")
  cluster.result$index = c("all","all","bic","icl")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 2, common noise",each = 4)
  df.clu = get.df.cluster(df,n.group = 20, knots.num = 20, nharm = 6)
  
  cluster.result$accuracy = get.clu.accuracy(df.clu, 20)
  
  df.bsp = get.clu.value(df, df.clu, 20, "bsp")
  df.fpc = get.clu.value(df, df.clu, 20, "fpc")
  df.bic = get.clu.value(df, df.clu, 20, "fem.bic")
  df.icl = get.clu.value(df, df.clu, 20, "fem.icl")
  cluster.result$avg.sd = c(df.bsp[1],df.fpc[1],df.bic[1],df.icl[1])
  cluster.result$avg.mudis = c(df.bsp[2],df.fpc[2],df.bic[2],df.icl[2])
  
  simuvid.round2.20.1 = rbind(simuvid.round2.20.1,cluster.result)
  seed = seed+100*i
}

simuvid.round2.20.2 = simuvid.round2.20.1
seed = 88888
for (i in  1:10){
  set.seed(seed)
  df = generate.curve(t.min=0,t.max=10,t.by=0.01,
                      n.group=20,err.sd=2)
  cluster.result = data.frame(matrix(ncol = 7, nrow = 4))
  colnames(cluster.result) <- c("method","index","n.group","round",
                                "accuracy","avg.sd","avg.mudis")
  cluster.result$method = c("bsp","fpc","fem","fem")
  cluster.result$index = c("all","all","bic","icl")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 2, common noise",each = 4)
  df.clu = get.df.cluster(df,n.group = 20, knots.num = 20, nharm = 6)
  
  cluster.result$accuracy = get.clu.accuracy(df.clu, 20)
  
  df.bsp = get.clu.value(df, df.clu, 20, "bsp")
  df.fpc = get.clu.value(df, df.clu, 20, "fpc")
  df.bic = get.clu.value(df, df.clu, 20, "fem.bic")
  df.icl = get.clu.value(df, df.clu, 20, "fem.icl")
  cluster.result$avg.sd = c(df.bsp[1],df.fpc[1],df.bic[1],df.icl[1])
  cluster.result$avg.mudis = c(df.bsp[2],df.fpc[2],df.bic[2],df.icl[2])
  
  simuvid.round2.20.2 = rbind(simuvid.round2.20.2,cluster.result)
  seed = seed+100*i
}
simuvid.round2.20.2 = na.omit(simuvid.round2.20.2)
write.csv(simuvid.round2.20.2[1:800,], file = "simuvid_round2_20.csv", row.names = FALSE)

# Round 3: n = 20, Sigma = 1, absolute-value proportion noise
simuvid.round3.20= data.frame(matrix(ncol = 7, nrow = 1))
colnames(simuvid.round3.20) <- c("method","index","n.group","round",
                                 "accuracy","avg.sd","avg.mudis")
seed = 8888
for (i in  1:50){ 
  set.seed(seed)
  df = generate.curve.ctvardiff(t.min=0,t.max=10,t.by=0.01,
                                n.group=20,err.sd=1)
  cluster.result = data.frame(matrix(ncol = 7, nrow = 4))
  colnames(cluster.result) <- c("method","index","n.group","round",
                                "accuracy","avg.sd","avg.mudis")
  cluster.result$method = c("bsp","fpc","fem","fem")
  cluster.result$index = c("all","all","bic","icl")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 1, varies noise",each = 4)
  df.clu = get.df.cluster(df,n.group = 20, knots.num = 20, nharm = 2)
  
  cluster.result$accuracy = get.clu.accuracy(df.clu, 20)
  
  df.bsp = get.clu.value(df, df.clu, 20, "bsp")
  df.fpc = get.clu.value(df, df.clu, 20, "fpc")
  df.bic = get.clu.value(df, df.clu, 20, "fem.bic")
  df.icl = get.clu.value(df, df.clu, 20, "fem.icl")
  cluster.result$avg.sd = c(df.bsp[1],df.fpc[1],df.bic[1],df.icl[1])
  cluster.result$avg.mudis = c(df.bsp[2],df.fpc[2],df.bic[2],df.icl[2])
  
  simuvid.round3.20 = rbind(simuvid.round3.20,cluster.result)
  seed = seed+10000
}

simuvid.round3.20.1 = simuvid.round3.20
seed = 888
for (i in  1:50){ 
  set.seed(seed)
  df = generate.curve.ctvardiff(t.min=0,t.max=10,t.by=0.01,
                                n.group=20,err.sd=1)
  cluster.result = data.frame(matrix(ncol = 7, nrow = 4))
  colnames(cluster.result) <- c("method","index","n.group","round",
                                "accuracy","avg.sd","avg.mudis")
  cluster.result$method = c("bsp","fpc","fem","fem")
  cluster.result$index = c("all","all","bic","icl")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 1, varies noise",each = 4)
  df.clu = get.df.cluster(df,n.group = 20, knots.num = 20, nharm = 2)
  
  cluster.result$accuracy = get.clu.accuracy(df.clu, 20)
  
  df.bsp = get.clu.value(df, df.clu, 20, "bsp")
  df.fpc = get.clu.value(df, df.clu, 20, "fpc")
  df.bic = get.clu.value(df, df.clu, 20, "fem.bic")
  df.icl = get.clu.value(df, df.clu, 20, "fem.icl")
  cluster.result$avg.sd = c(df.bsp[1],df.fpc[1],df.bic[1],df.icl[1])
  cluster.result$avg.mudis = c(df.bsp[2],df.fpc[2],df.bic[2],df.icl[2])
  
  simuvid.round3.20.1 = rbind(simuvid.round3.20.1,cluster.result)
  seed = seed+1000
}

simuvid.round3.20.2 = simuvid.round3.20.1
seed = 88888
for (i in  1:50){ 
  set.seed(seed)
  df = generate.curve.ctvardiff(t.min=0,t.max=10,t.by=0.01,
                                n.group=20,err.sd=1)
  cluster.result = data.frame(matrix(ncol = 7, nrow = 4))
  colnames(cluster.result) <- c("method","index","n.group","round",
                                "accuracy","avg.sd","avg.mudis")
  cluster.result$method = c("bsp","fpc","fem","fem")
  cluster.result$index = c("all","all","bic","icl")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 1, varies noise",each = 4)
  df.clu = get.df.cluster(df,n.group = 20, knots.num = 20, nharm = 2)
  
  cluster.result$accuracy = get.clu.accuracy(df.clu, 20)
  
  df.bsp = get.clu.value(df, df.clu, 20, "bsp")
  df.fpc = get.clu.value(df, df.clu, 20, "fpc")
  df.bic = get.clu.value(df, df.clu, 20, "fem.bic")
  df.icl = get.clu.value(df, df.clu, 20, "fem.icl")
  cluster.result$avg.sd = c(df.bsp[1],df.fpc[1],df.bic[1],df.icl[1])
  cluster.result$avg.mudis = c(df.bsp[2],df.fpc[2],df.bic[2],df.icl[2])
  
  simuvid.round3.20.2 = rbind(simuvid.round3.20.2,cluster.result)
  seed = seed+1000
}

simuvid.round3.20.3 = simuvid.round3.20.2
seed = 888888
for (i in  1:43){ ## stop at 44, finish 1:43
  set.seed(seed)
  df = generate.curve.ctvardiff(t.min=0,t.max=10,t.by=0.01,
                                n.group=20,err.sd=1)
  cluster.result = data.frame(matrix(ncol = 7, nrow = 4))
  colnames(cluster.result) <- c("method","index","n.group","round",
                                "accuracy","avg.sd","avg.mudis")
  cluster.result$method = c("bsp","fpc","fem","fem")
  cluster.result$index = c("all","all","bic","icl")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 1, varies noise",each = 4)
  df.clu = get.df.cluster(df,n.group = 20, knots.num = 20, nharm = 2)
  
  cluster.result$accuracy = get.clu.accuracy(df.clu, 20)
  
  df.bsp = get.clu.value(df, df.clu, 20, "bsp")
  df.fpc = get.clu.value(df, df.clu, 20, "fpc")
  df.bic = get.clu.value(df, df.clu, 20, "fem.bic")
  df.icl = get.clu.value(df, df.clu, 20, "fem.icl")
  cluster.result$avg.sd = c(df.bsp[1],df.fpc[1],df.bic[1],df.icl[1])
  cluster.result$avg.mudis = c(df.bsp[2],df.fpc[2],df.bic[2],df.icl[2])
  
  simuvid.round3.20.3 = rbind(simuvid.round3.20.3,cluster.result)
  seed = seed+1000
}

simuvid.round3.20.4 = simuvid.round3.20.3
seed = 66
for (i in  1:7){ 
  set.seed(seed)
  df = generate.curve.ctvardiff(t.min=0,t.max=10,t.by=0.01,
                                n.group=20,err.sd=1)
  cluster.result = data.frame(matrix(ncol = 7, nrow = 4))
  colnames(cluster.result) <- c("method","index","n.group","round",
                                "accuracy","avg.sd","avg.mudis")
  cluster.result$method = c("bsp","fpc","fem","fem")
  cluster.result$index = c("all","all","bic","icl")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 1, varies noise",each = 4)
  df.clu = get.df.cluster(df,n.group = 20, knots.num = 20, nharm = 2)
  
  cluster.result$accuracy = get.clu.accuracy(df.clu, 20)
  
  df.bsp = get.clu.value(df, df.clu, 20, "bsp")
  df.fpc = get.clu.value(df, df.clu, 20, "fpc")
  df.bic = get.clu.value(df, df.clu, 20, "fem.bic")
  df.icl = get.clu.value(df, df.clu, 20, "fem.icl")
  cluster.result$avg.sd = c(df.bsp[1],df.fpc[1],df.bic[1],df.icl[1])
  cluster.result$avg.mudis = c(df.bsp[2],df.fpc[2],df.bic[2],df.icl[2])
  
  simuvid.round3.20.4 = rbind(simuvid.round3.20.4,cluster.result)
  seed = seed+1000
}

simuvid.round3.20.4 = na.omit(simuvid.round3.20.4)
write.csv(simuvid.round3.20.4, file = "simuvid_round3_20.csv", row.names = FALSE)

# Round 4: n = 20, Sigma = 2, absolute-value proportion noise
simuvid.round4.20= data.frame(matrix(ncol = 7, nrow = 1))
colnames(simuvid.round4.20) <- c("method","index","n.group","round",
                                 "accuracy","avg.sd","avg.mudis")
seed = 888
for (i in  1:200){
  set.seed(seed)
  df = generate.curve.ctvardiff(t.min=0,t.max=10,t.by=0.01,
                                n.group=20,err.sd=2)
  cluster.result = data.frame(matrix(ncol = 7, nrow = 4))
  colnames(cluster.result) <- c("method","index","n.group","round",
                                "accuracy","avg.sd","avg.mudis")
  cluster.result$method = c("bsp","fpc","fem","fem")
  cluster.result$index = c("all","all","bic","icl")
  cluster.result$n.group = rep(20,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 2, varies noise",each = 4)
  df.clu = get.df.cluster(df,n.group = 20, knots.num = 20, nharm = 4)
  
  cluster.result$accuracy = get.clu.accuracy(df.clu, 20)
  
  df.bsp = get.clu.value(df, df.clu, 20, "bsp")
  df.fpc = get.clu.value(df, df.clu, 20, "fpc")
  df.bic = get.clu.value(df, df.clu, 20, "fem.bic")
  df.icl = get.clu.value(df, df.clu, 20, "fem.icl")
  cluster.result$avg.sd = c(df.bsp[1],df.fpc[1],df.bic[1],df.icl[1])
  cluster.result$avg.mudis = c(df.bsp[2],df.fpc[2],df.bic[2],df.icl[2])
  
  simuvid.round4.20 = rbind(simuvid.round4.20,cluster.result)
  seed = seed+10000*i
}
simuvid.round4.20 = na.omit(simuvid.round4.20)
write.csv(simuvid.round4.20, file = "simuvid_round4_20.csv", row.names = FALSE)


### Summary analysis of the cluster validation
r1.20.raw = read.csv("simuvid_round1_20.csv",header = TRUE)

r2.20.raw = read.csv("simuvid_round2_20.csv",header = TRUE)

r3.20.raw = read.csv("simuvid_round3_20.csv",header = TRUE)

r4.20.raw = read.csv("simuvid_round4_20.csv",header = TRUE)

#round 1
r1.20 = r1.20.raw%>%group_by(round, n.group, method,
                             index)%>%summarize(mean.accuracy = mean(accuracy),
                                                sd.accuracy = sd(accuracy),
                                                mean.sd = mean(avg.sd),
                                                mean.mudis = mean(avg.mudis))
r1.20

#round 2
r2.20 = r2.20.raw%>%group_by(round, n.group, method,
                             index)%>%summarize(mean.accuracy = mean(accuracy),
                                                sd.accuracy = sd(accuracy),
                                                mean.sd = mean(avg.sd),
                                                mean.mudis = mean(avg.mudis))
r2.20

#round 3
r3.20 = r3.20.raw%>%group_by(round, n.group, method,
                             index)%>%summarize(mean.accuracy = mean(accuracy),
                                                sd.accuracy = sd(accuracy),
                                                mean.sd = mean(avg.sd),
                                                mean.mudis = mean(avg.mudis))
r3.20

#round 4
r4.20 = r4.20.raw%>%group_by(round, n.group, method,
                             index)%>%summarize(mean.accuracy = mean(accuracy),
                                                sd.accuracy = sd(accuracy),
                                                mean.sd = mean(avg.sd),
                                                mean.mudis = mean(avg.mudis))
r4.20

