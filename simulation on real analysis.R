## This is the R file including all the content related to the simulated data based on real data analysis.
## Please run all of it from the beginning to the end by order
## The main content of this file includes:
##     0. Installing required packages 
##     1. Functions for Data simulation and Analysis
##     2. Clustering Number Selection Validation
##     3. Clustering Validation Study
######################################################################################


################################## Installing required packages ######################
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
              "ggthemes","gridExtra","factoextra", "fiftystater","randomForest", 
              "NbClust","funFEM","fda.usc","tclust","dtwclust","fpc")
ipak(packages)
## The package fiftystater is for the U.S map visualization, and it may be unavailable on CRAN
## So its achieved file is attached in the Github repository with the path 



################################ Loading the unsmoothed non-missing data and cluster information####################
mindiff.date.1 = read.csv("mindiff_day1_nona.csv")
mindiff.date.1$date = as.Date(mindiff.date.1$x,format = "%Y-%m-%d")
mindiff.state.1 = read.csv("mindiff_state1_nona.csv")
cluster.result = read.csv("real_analysis_result.csv")

## Generate curves based on real data
generate.curve.real<-function(n.group, plotting = F){
  ## Simulating 4 clusters based on real application result
  # IN cluster 1
  clu.1 = cluster.result[cluster.result$fem.4 == 1,]$state
  clu.1.raw = mindiff.state.1[,names(mindiff.state.1)%in%clu.1]
  clu.1.rawmat.1 = matrix(unlist(clu.1.raw),ncol = ncol(clu.1.raw))
  bks.1 = as.Date(quantile(unclass(as.Date(mindiff.date.1$date)), seq(0,1,length = 15)), origin = "1970-01-01")
  err1.basis = create.bspline.basis(rangeval = c(min(mindiff.date.1$date),max(mindiff.date.1$date)), breaks = bks.1, norder = 4) 
  clu.1.fd<-smooth.basis(argvals = mindiff.date.1$date, y = clu.1.rawmat.1, err1.basis)$fd
  clu.1.fpca = pca.fd(clu.1.fd,nharm=3) 
  clu.1.fpca.loading = clu.1.fpca$harmonics$coefs
  clu.1.fpca.lambda = clu.1.fpca$values[1:3]
  clu.1.mean.coef = clu.1.fpca$meanfd$coefs
  
  ## Simulate n curves:
  clu.1.score1 = rnorm(n.group,0,sqrt(clu.1.fpca.lambda[1]))
  clu.1.score2 = rnorm(n.group,0,sqrt(clu.1.fpca.lambda[2]))
  clu.1.score3 = rnorm(n.group,0,sqrt(clu.1.fpca.lambda[3]))
  clu.1.scores = rbind(clu.1.score1,clu.1.score2,clu.1.score3)
  clu.1.simcoeff = t(clu.1.scores)%*%t(clu.1.fpca.loading)+matrix(clu.1.mean.coef,nrow=n.group,ncol=length(clu.1.mean.coef),byrow=TRUE)
  clu.1.sim.df = clu.1.fd
  clu.1.sim.df$coefs = t(clu.1.simcoeff)
  clu.1.sim.data = predict(clu.1.sim.df,mindiff.date.1$date)
 
  name.clu1 = paste0("clu1.",1:n.group)
  colnames(clu.1.sim.data) <- name.clu1
  clu.1.sim.data = data.frame(clu.1.sim.data)
  
  # IN cluster 2
  clu.2 = cluster.result[cluster.result$fem.4 == 2,]$state
  clu.2.raw = mindiff.state.1[,names(mindiff.state.1)%in%clu.2]
  clu.2.rawmat.1 = matrix(unlist(clu.2.raw),ncol = ncol(clu.2.raw))
  clu.2.fd<-smooth.basis(argvals = mindiff.date.1$date, y = clu.2.rawmat.1, err1.basis)$fd
  clu.2.fpca = pca.fd(clu.2.fd,nharm=3) 
  clu.2.fpca.loading = clu.2.fpca$harmonics$coefs
  clu.2.fpca.lambda = clu.2.fpca$values[1:3]
  clu.2.mean.coef = clu.2.fpca$meanfd$coefs
  
  ## Simulate n curves:
  clu.2.score1 = rnorm(n.group,0,sqrt(clu.2.fpca.lambda[1]))
  clu.2.score2 = rnorm(n.group,0,sqrt(clu.2.fpca.lambda[2]))
  clu.2.score3 = rnorm(n.group,0,sqrt(clu.2.fpca.lambda[3]))
  clu.2.scores = rbind(clu.2.score1,clu.2.score2,clu.2.score3)
  clu.2.simcoeff = t(clu.2.scores)%*%t(clu.2.fpca.loading)+matrix(clu.2.mean.coef,nrow=n.group,ncol=length(clu.2.mean.coef),byrow=TRUE)
  clu.2.sim.df = clu.2.fd
  clu.2.sim.df$coefs = t(clu.2.simcoeff)
  clu.2.sim.data = predict(clu.2.sim.df,mindiff.date.1$date)
  
  name.clu2 = paste0("clu2.",1:n.group)
  colnames(clu.2.sim.data) <- name.clu2
  clu.2.sim.data = data.frame(clu.2.sim.data)
  
  # IN cluster 3
  clu.3 = cluster.result[cluster.result$fem.4 == 3,]$state
  clu.3.raw = mindiff.state.1[,names(mindiff.state.1)%in%clu.3]
  clu.3.rawmat.1 = matrix(unlist(clu.3.raw),ncol = ncol(clu.3.raw))
  clu.3.fd<-smooth.basis(argvals = mindiff.date.1$date, y = clu.3.rawmat.1, err1.basis)$fd
  clu.3.fpca = pca.fd(clu.3.fd,nharm=3) 
  clu.3.fpca.loading = clu.3.fpca$harmonics$coefs
  clu.3.fpca.lambda = clu.3.fpca$values[1:3]
  clu.3.mean.coef = clu.3.fpca$meanfd$coefs
  
  ## Simulate n curves:
  clu.3.score1 = rnorm(n.group,0,sqrt(clu.3.fpca.lambda[1]))
  clu.3.score2 = rnorm(n.group,0,sqrt(clu.3.fpca.lambda[2]))
  clu.3.score3 = rnorm(n.group,0,sqrt(clu.3.fpca.lambda[3]))
  clu.3.scores = rbind(clu.3.score1,clu.3.score2,clu.3.score3)
  clu.3.simcoeff = t(clu.3.scores)%*%t(clu.3.fpca.loading)+matrix(clu.3.mean.coef,nrow=n.group,ncol=length(clu.3.mean.coef),byrow=TRUE)
  clu.3.sim.df = clu.3.fd
  clu.3.sim.df$coefs = t(clu.3.simcoeff)
  clu.3.sim.data = predict(clu.3.sim.df,mindiff.date.1$date)
  
  name.clu3 = paste0("clu3.",1:n.group)
  colnames(clu.3.sim.data) <- name.clu3
  clu.3.sim.data = data.frame(clu.3.sim.data)
  
  # IN cluster 4
  clu.4 = cluster.result[cluster.result$fem.4 == 4,]$state
  clu.4.raw = mindiff.state.1[,names(mindiff.state.1)%in%clu.4]
  clu.4.rawmat.1 = matrix(unlist(clu.4.raw),ncol = ncol(clu.4.raw))
  clu.4.fd<-smooth.basis(argvals = mindiff.date.1$date, y = clu.4.rawmat.1, err1.basis)$fd
  clu.4.fpca = pca.fd(clu.4.fd,nharm=3) 
  clu.4.fpca.loading = clu.4.fpca$harmonics$coefs
  clu.4.fpca.lambda = clu.4.fpca$values[1:3]
  clu.4.mean.coef = clu.4.fpca$meanfd$coefs
  
  ## Simulate n curves:
  clu.4.score1 = rnorm(n.group,0,sqrt(clu.4.fpca.lambda[1]))
  clu.4.score2 = rnorm(n.group,0,sqrt(clu.4.fpca.lambda[2]))
  clu.4.score3 = rnorm(n.group,0,sqrt(clu.4.fpca.lambda[3]))
  clu.4.scores = rbind(clu.4.score1,clu.4.score2,clu.4.score3)
  clu.4.simcoeff = t(clu.4.scores)%*%t(clu.4.fpca.loading)+matrix(clu.4.mean.coef,nrow=n.group,ncol=length(clu.4.mean.coef),byrow=TRUE)
  clu.4.sim.df = clu.4.fd
  clu.4.sim.df$coefs = t(clu.4.simcoeff)
  clu.4.sim.data = predict(clu.4.sim.df,mindiff.date.1$date)
  
  name.clu4 = paste0("clu4.",1:n.group)
  colnames(clu.4.sim.data) <- name.clu4
  clu.4.sim.data = data.frame(clu.4.sim.data)
  
  simu.all = cbind(clu.1.sim.data,clu.2.sim.data,clu.3.sim.data,clu.4.sim.data)
  simu.all$t = mindiff.date.1$date
  
  if(plotting == T){
    par(mfrow=c(2,2))
    plot(clu.1.sim.df,ylim = c(0,10))
    plot(clu.2.sim.df,ylim = c(0,10))
    plot(clu.3.sim.df,ylim = c(0,10))
    plot(clu.4.sim.df,ylim = c(0,10))
  }
  
  return (simu.all)
}


######################## Functions for Cluster Number Selection Study #################################
cluster.select.num.real = function(simu.all,n.group,plotting = FALSE){
  t = simu.all$t
  simu.all = simu.all[,-which(names(simu.all) %in% c("t"))]
  simu.all.mat = matrix(unlist(simu.all),ncol = 4*n.group)
  
  bks.1 = as.Date(quantile(unclass(as.Date(mindiff.date.1$date)), seq(0,1,length = 12)), origin = "1970-01-01")
  err1.basis = create.bspline.basis(rangeval = c(min(mindiff.date.1$date),max(mindiff.date.1$date)), breaks = bks.1, norder = 4) 
  simu.fd<-smooth.basis(argvals = t, y = simu.all.mat, err1.basis)$fd
  
  if(plotting == TRUE){plot(simu.fd)}
  
  ## FunFEM
  fem.clu.bic = funFEM(simu.fd,K = 2:6,model = "all",crit = "bic", init = "kmeans",eps = 1e-5)
  fem.num.bic = fem.clu.bic$K
  
  fem.clu.icl = funFEM(simu.fd,K = 2:6,model = "all",crit = "icl", init = "kmeans",eps = 1e-5)
  fem.num.icl = fem.clu.icl$K
  
  ## kmeans on smoothed FPC scores
  range = c(as.numeric(min(t)),as.numeric(max(t)))
  bspFpcPen = fdPar(fdobj=err1.basis, Lfdobj=2, lambda = 2*10^6)
  smooth.fpc = pca.fd(simu.fd,nharm=3,harmfdPar=bspFpcPen)
  fd.scores = smooth.fpc$scores
  # all index  
  fpc.kmeans.best <- NbClust(fd.scores, distance = "euclidean", 
                             min.nc=2, max.nc=6, 
                             method = "kmeans", 
                             index = "all")
  fpc.clu.select =  max(fpc.kmeans.best$Best.partition)
  
  ## kmeans on coef
  fd.coefs = t(simu.fd$coefs)
  # all index  
  coeff.kmeans.best <- NbClust(fd.coefs, distance = "euclidean", 
                               min.nc=2, max.nc=6, 
                               method = "kmeans", 
                               index = "all")
  bsp.clu.select =  max(coeff.kmeans.best$Best.partition)
  
  simu.clunum.compare = rbind(fem.num.bic,fem.num.icl,
                              bsp.clu.select,fpc.clu.select)
  return(simu.clunum.compare)
}

#################Clustering Number Selection Validation Simulation ##########################
simunum.real.20= data.frame(matrix(ncol = 4, nrow = 1))
colnames(simunum.real.20) <- c("method","index","clu.num","n.group")
seed = 888
for (i in 1:52){ 
  set.seed(seed)
  df = generate.curve.real(20)
  clunum.result = data.frame(matrix(ncol = 4, nrow = 4))
  colnames(clunum.result) <- c("method","index","clu.num","n.group")
  clunum.result$method = c("fem","fem","bsp","fpc")
  clunum.result$index = c("bic","icl","all","all")
  clunum.result$clu.num = cluster.select.num.real(simu.all=df,n.group=20)
  clunum.result$n.group = rep(20,each = 4)
  
  simunum.real.20 = rbind(simunum.real.20,clunum.result)
  seed = seed+10000*i
}


seed = 8888
for (i in 1:42){
  set.seed(seed)
  df = generate.curve.real(20)
  clunum.result = data.frame(matrix(ncol = 4, nrow = 4))
  colnames(clunum.result) <- c("method","index","clu.num","n.group")
  clunum.result$method = c("fem","fem","bsp","fpc")
  clunum.result$index = c("bic","icl","all","all")
  clunum.result$clu.num = cluster.select.num.real(simu.all=df,n.group=20)
  clunum.result$n.group = rep(20,each = 4)
  
  simunum.real.20 = rbind(simunum.real.20,clunum.result)
  seed = seed+10000*i
}

seed = 666
for (i in 1:11){
  set.seed(seed)
  df = generate.curve.real(20)
  clunum.result = data.frame(matrix(ncol = 4, nrow = 4))
  colnames(clunum.result) <- c("method","index","clu.num","n.group")
  clunum.result$method = c("fem","fem","bsp","fpc")
  clunum.result$index = c("bic","icl","all","all")
  clunum.result$clu.num = cluster.select.num.real(simu.all=df,n.group=20)
  clunum.result$n.group = rep(20,each = 4)
  
  simunum.real.20 = rbind(simunum.real.20,clunum.result)
  seed = seed+10000*i
}

seed = 6666
for (i in 1:11){
  set.seed(seed)
  df = generate.curve.real(20)
  clunum.result = data.frame(matrix(ncol = 4, nrow = 4))
  colnames(clunum.result) <- c("method","index","clu.num","n.group")
  clunum.result$method = c("fem","fem","bsp","fpc")
  clunum.result$index = c("bic","icl","all","all")
  clunum.result$clu.num = cluster.select.num.real(simu.all=df,n.group=20)
  clunum.result$n.group = rep(20,each = 4)
  
  simunum.real.20 = rbind(simunum.real.20,clunum.result)
  seed = seed+10000*i
}
simunum.real.20 = na.omit(simunum.real.20)

seed = 68
for (i in 1:41){
  set.seed(seed)
  df = generate.curve.real(20)
  clunum.result = data.frame(matrix(ncol = 4, nrow = 4))
  colnames(clunum.result) <- c("method","index","clu.num","n.group")
  clunum.result$method = c("fem","fem","bsp","fpc")
  clunum.result$index = c("bic","icl","all","all")
  clunum.result$clu.num = cluster.select.num.real(simu.all=df,n.group=20)
  clunum.result$n.group = rep(20,each = 4)
  
  simunum.real.20 = rbind(simunum.real.20,clunum.result)
  seed = seed+10000*i
}

seed = 86
for (i in 1:43){
  set.seed(seed)
  df = generate.curve.real(20)
  clunum.result = data.frame(matrix(ncol = 4, nrow = 4))
  colnames(clunum.result) <- c("method","index","clu.num","n.group")
  clunum.result$method = c("fem","fem","bsp","fpc")
  clunum.result$index = c("bic","icl","all","all")
  clunum.result$clu.num = cluster.select.num.real(simu.all=df,n.group=20)
  clunum.result$n.group = rep(20,each = 4)
  
  simunum.real.20 = rbind(simunum.real.20,clunum.result)
  seed = seed+10000*i
}
simunum.real.20 = na.omit(simunum.real.20)
write.csv(simunum.real.20, file = "simunum_real_20.csv", row.names = FALSE)

### Summary analysis of the cluster number detection.
clunum.sim.20 = read.csv("simunum_real_20.csv",header = TRUE)
clunum.sim.20.summary = clunum.sim.20%>%group_by(n.group, method,
                             index, clu.num)%>%summarize(count = n())
clunum.sim.20.summary

############################## Clustering Validation Study Functions #####################################
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
get.df.cluster = function(simu.all,n.group){
  t = simu.all$t
  simu.all = simu.all[,-which(names(simu.all) %in% c("t"))]
  simu.all.mat = matrix(unlist(simu.all),ncol = 4*n.group)
  
  bks.1 = as.Date(quantile(unclass(as.Date(mindiff.date.1$date)), seq(0,1,length = 12)), origin = "1970-01-01")
  simu.basis = create.bspline.basis(rangeval = c(min(mindiff.date.1$date),max(mindiff.date.1$date)), breaks = bks.1, norder = 4) 
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
  bspFpcPen = fdPar(fdobj=simu.basis, Lfdobj=2, lambda = 2e6)
  ## We will explain this lambda selection in the next part
  smooth.fpc = pca.fd(simu.fd,nharm=3,harmfdPar=bspFpcPen)
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


#################Clustering Number Selection Validation Simulation##########################
simuvid.real.20= data.frame(matrix(ncol = 4, nrow = 1))
colnames(simuvid.real.20) <- c("method","index","n.group","accuracy")
seed = 888
for (i in  1:200){
  set.seed(seed)
  df = generate.curve.real(20)
  cluvid.result = data.frame(matrix(ncol = 4, nrow = 4))
  colnames(cluvid.result) <- c("method","index","n.group","accuracy")
  cluvid.result$method = c("bsp","fpc","fem","fem")
  cluvid.result$index = c("all","all","bic","icl")
  cluvid.result$n.group = rep(20,each = 4)
  df.clu = get.df.cluster(df,n.group = 20)
  
  cluvid.result$accuracy = get.clu.accuracy(df.clu, 20)
  
  simuvid.real.20 = rbind(simuvid.real.20,cluvid.result)
  seed = seed+10000*i
}
simuvid.real.20 = na.omit(simuvid.real.20)
write.csv(simuvid.real.20, file = "simuvid_real_20.csv", row.names = FALSE)

### Summary the validation result
simuvid.raw = read.csv("simuvid_real_20.csv",header = TRUE)
simuvid.raw.summary = simuvid.raw%>%
  group_by(n.group, method,index)%>%
  summarize(mean.accuracy = mean(accuracy),
            sd.accuracy = sd(accuracy),
            se.accuracy = sd(accuracy)/sqrt(200))
simuvid.raw.summary


####################################### Kmeans on Raw Data #######################################
##### Select Number of Cluster #####
simunum.realraw.20= data.frame(matrix(ncol = 4, nrow = 1))
colnames(simunum.realraw.20) <- c("method","index","clu.num","n.group")
seed = 888
for (i in 1:52){
  set.seed(seed)
  df = generate.curve.real(20)
  df = df[,-which(names(df) %in% c("t"))]
  clunum.result.raw = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(clunum.result.raw) <- c("method","index","clu.num","n.group")
  clunum.result.raw$method = "kmeans on raw"
  clunum.result.raw$index = "average silhouette width"
  pamk.best <- pamk(df,krange=2:6)
  clunum.result.raw$clu.num = pamk.best$nc
  clunum.result.raw$n.group = 20
  
  simunum.realraw.20 = rbind(simunum.realraw.20,clunum.result.raw)
  seed = seed+10000*i
}

seed = 8888
for (i in 1:42){
  set.seed(seed)
  df = generate.curve.real(20)
  df = df[,-which(names(df) %in% c("t"))]
  clunum.result.raw = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(clunum.result.raw) <- c("method","index","clu.num","n.group")
  clunum.result.raw$method = "kmeans on raw"
  clunum.result.raw$index = "average silhouette width"
  pamk.best <- pamk(df,krange=2:6)
  clunum.result.raw$clu.num = pamk.best$nc
  clunum.result.raw$n.group = 20
  
  simunum.realraw.20 = rbind(simunum.realraw.20,clunum.result.raw)
  seed = seed+10000*i
}

seed = 666
for (i in 1:11){
  set.seed(seed)
  df = generate.curve.real(20)
  df = df[,-which(names(df) %in% c("t"))]
  clunum.result.raw = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(clunum.result.raw) <- c("method","index","clu.num","n.group")
  clunum.result.raw$method = "kmeans on raw"
  clunum.result.raw$index = "average silhouette width"
  pamk.best <- pamk(df,krange=2:6)
  clunum.result.raw$clu.num = pamk.best$nc
  clunum.result.raw$n.group = 20
  
  simunum.realraw.20 = rbind(simunum.realraw.20,clunum.result.raw)
  seed = seed+10000*i
}

seed = 6666
for (i in 1:11){
  set.seed(seed)
  df = generate.curve.real(20)
  df = df[,-which(names(df) %in% c("t"))]
  clunum.result.raw = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(clunum.result.raw) <- c("method","index","clu.num","n.group")
  clunum.result.raw$method = "kmeans on raw"
  clunum.result.raw$index = "average silhouette width"
  pamk.best <- pamk(df,krange=2:6)
  clunum.result.raw$clu.num = pamk.best$nc
  clunum.result.raw$n.group = 20
  
  simunum.realraw.20 = rbind(simunum.realraw.20,clunum.result.raw)
  seed = seed+10000*i
}

seed = 68
for (i in 1:41){
  set.seed(seed)
  df = generate.curve.real(20)
  df = df[,-which(names(df) %in% c("t"))]
  clunum.result.raw = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(clunum.result.raw) <- c("method","index","clu.num","n.group")
  clunum.result.raw$method = "kmeans on raw"
  clunum.result.raw$index = "average silhouette width"
  pamk.best <- pamk(df,krange=2:6)
  clunum.result.raw$clu.num = pamk.best$nc
  clunum.result.raw$n.group = 20
  
  simunum.realraw.20 = rbind(simunum.realraw.20,clunum.result.raw)
  seed = seed+10000*i
}

seed = 86
for (i in 1:43){
  set.seed(seed)
  df = generate.curve.real(20)
  df = df[,-which(names(df) %in% c("t"))]
  clunum.result.raw = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(clunum.result.raw) <- c("method","index","clu.num","n.group")
  clunum.result.raw$method = "kmeans on raw"
  clunum.result.raw$index = "average silhouette width"
  pamk.best <- pamk(df,krange=2:6)
  clunum.result.raw$clu.num = pamk.best$nc
  clunum.result.raw$n.group = 20
  
  simunum.realraw.20 = rbind(simunum.realraw.20,clunum.result.raw)
  seed = seed+10000*i
}
simunum.realraw.20 = na.omit(simunum.realraw.20)
write.csv(simunum.realraw.20, file = "simunum_realraw_20.csv", row.names = FALSE)

### Summary analysis of the cluster number detection.
clunum.sim.20.raw = read.csv("simunum_realraw_20.csv",header = TRUE)
clunum.sim.20.raw.summary = clunum.sim.20.raw%>%group_by(n.group, method,
                                                 index, clu.num)%>%summarize(count = n())
clunum.sim.20.raw.summary

################### Cluster Result Validation #############################
get.clu.accuracy.raw = function(df.clu, n.group){
  # Hungarian Alg
  real = rep(c(1,2,3,4),each = n.group)
  a = minWeightBipartiteMatching(df.clu, real)
  df.clu.1 = a[df.clu]
  
  raw.accu = sum(real == df.clu.1)/(n.group*4) 
  return(raw.accu)
}

#### Start the Simulation
simuvid.realraw.20= data.frame(matrix(ncol = 4, nrow = 1))
colnames(simuvid.realraw.20) <- c("method","index","n.group","accuracy")
seed = 88
for (i in  1:200){
  set.seed(seed)
  df = generate.curve.real(20)
  df = df[,-which(names(df) %in% c("t"))]
  cluvid.raw.result = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(cluvid.raw.result) <- c("method","index","n.group","accuracy")
  cluvid.raw.result$method = c("kmeans on raw")
  cluvid.raw.result$index = c("average silhouette width")
  cluvid.raw.result$n.group = 20
  df.clu <- kmeans(t(as.matrix(df)),4)
  cluvid.raw.clu = df.clu$cluster
  
  cluvid.raw.result$accuracy = get.clu.accuracy.raw(cluvid.raw.clu, 20)
  
  simuvid.realraw.20 = rbind(simuvid.realraw.20,cluvid.raw.result)
  seed = seed+10000*i
}
simuvid.realraw.20 = na.omit(simuvid.realraw.20)
write.csv(simuvid.realraw.20, file = "simuvid_realraw_20.csv", row.names = FALSE)

### Summary the validation result
simuvid.realraw = read.csv("simuvid_realraw_20.csv",header = TRUE)
simuvid.realraw.summary = simuvid.realraw%>%
  group_by(n.group, method,index)%>%
  summarize(mean.accuracy = mean(accuracy),
            sd.accuracy = sd(accuracy),
            se.accuracy = sd(accuracy)/sqrt(200))
simuvid.realraw.summary

############## Kmeans on Raw Data with Noise #######################################
generate.curve.real.noise<-function(n.group, plotting = F){
  ## Simulating 4 clusters based on real application result
  # IN cluster 1
  clu.1 = cluster.result[cluster.result$fem.4 == 1,]$state
  clu.1.raw = mindiff.state.1[,names(mindiff.state.1)%in%clu.1]
  clu.1.rawmat.1 = matrix(unlist(clu.1.raw),ncol = 12)
  bks.1 = as.Date(quantile(unclass(as.Date(mindiff.date.1$date)), seq(0,1,length = 15)), origin = "1970-01-01")
  err1.basis = create.bspline.basis(rangeval = c(min(mindiff.date.1$date),max(mindiff.date.1$date)), breaks = bks.1, norder = 4) 
  clu.1.fd<-smooth.basis(argvals = mindiff.date.1$date, y = clu.1.rawmat.1, err1.basis)$fd
  clu.1.fpca = pca.fd(clu.1.fd,nharm=3) 
  clu.1.fpca.loading = clu.1.fpca$harmonics$coefs
  clu.1.fpca.lambda = clu.1.fpca$values[1:3]
  clu.1.mean.coef = clu.1.fpca$meanfd$coefs
  
  ## Simulate n curves:
  clu.1.score1 = rnorm(n.group,0,sqrt(clu.1.fpca.lambda[1]))
  clu.1.score2 = rnorm(n.group,0,sqrt(clu.1.fpca.lambda[2]))
  clu.1.score3 = rnorm(n.group,0,sqrt(clu.1.fpca.lambda[3]))
  clu.1.scores = rbind(clu.1.score1,clu.1.score2,clu.1.score3)
  clu.1.simcoeff = t(clu.1.scores)%*%t(clu.1.fpca.loading)+matrix(clu.1.mean.coef,nrow=n.group,ncol=length(clu.1.mean.coef),byrow=TRUE)
  clu.1.sim.df = clu.1.fd
  clu.1.sim.df$coefs = t(clu.1.simcoeff)
  clu.1.sim.data = predict(clu.1.sim.df,mindiff.date.1$date)
  
  name.clu1 = paste0("clu1.",1:n.group)
  colnames(clu.1.sim.data) <- name.clu1
  clu.1.sim.data = data.frame(clu.1.sim.data)
  
  # IN cluster 2
  clu.2 = cluster.result[cluster.result$fem.4 == 2,]$state
  clu.2.raw = mindiff.state.1[,names(mindiff.state.1)%in%clu.2]
  clu.2.rawmat.1 = matrix(unlist(clu.2.raw),ncol = 7)
  clu.2.fd<-smooth.basis(argvals = mindiff.date.1$date, y = clu.2.rawmat.1, err1.basis)$fd
  clu.2.fpca = pca.fd(clu.2.fd,nharm=3) 
  clu.2.fpca.loading = clu.2.fpca$harmonics$coefs
  clu.2.fpca.lambda = clu.2.fpca$values[1:3]
  clu.2.mean.coef = clu.2.fpca$meanfd$coefs
  
  ## Simulate n curves:
  clu.2.score1 = rnorm(n.group,0,sqrt(clu.2.fpca.lambda[1]))
  clu.2.score2 = rnorm(n.group,0,sqrt(clu.2.fpca.lambda[2]))
  clu.2.score3 = rnorm(n.group,0,sqrt(clu.2.fpca.lambda[3]))
  clu.2.scores = rbind(clu.2.score1,clu.2.score2,clu.2.score3)
  clu.2.simcoeff = t(clu.2.scores)%*%t(clu.2.fpca.loading)+matrix(clu.2.mean.coef,nrow=n.group,ncol=length(clu.2.mean.coef),byrow=TRUE)
  clu.2.sim.df = clu.2.fd
  clu.2.sim.df$coefs = t(clu.2.simcoeff)
  clu.2.sim.data = predict(clu.2.sim.df,mindiff.date.1$date)
  
  name.clu2 = paste0("clu2.",1:n.group)
  colnames(clu.2.sim.data) <- name.clu2
  clu.2.sim.data = data.frame(clu.2.sim.data)
  
  # IN cluster 3
  clu.3 = cluster.result[cluster.result$fem.4 == 3,]$state
  clu.3.raw = mindiff.state.1[,names(mindiff.state.1)%in%clu.3]
  clu.3.rawmat.1 = matrix(unlist(clu.3.raw),ncol = 5)
  clu.3.fd<-smooth.basis(argvals = mindiff.date.1$date, y = clu.3.rawmat.1, err1.basis)$fd
  clu.3.fpca = pca.fd(clu.3.fd,nharm=3) 
  clu.3.fpca.loading = clu.3.fpca$harmonics$coefs
  clu.3.fpca.lambda = clu.3.fpca$values[1:3]
  clu.3.mean.coef = clu.3.fpca$meanfd$coefs
  
  ## Simulate n curves:
  clu.3.score1 = rnorm(n.group,0,sqrt(clu.3.fpca.lambda[1]))
  clu.3.score2 = rnorm(n.group,0,sqrt(clu.3.fpca.lambda[2]))
  clu.3.score3 = rnorm(n.group,0,sqrt(clu.3.fpca.lambda[3]))
  clu.3.scores = rbind(clu.3.score1,clu.3.score2,clu.3.score3)
  clu.3.simcoeff = t(clu.3.scores)%*%t(clu.3.fpca.loading)+matrix(clu.3.mean.coef,nrow=n.group,ncol=length(clu.3.mean.coef),byrow=TRUE)
  clu.3.sim.df = clu.3.fd
  clu.3.sim.df$coefs = t(clu.3.simcoeff)
  clu.3.sim.data = predict(clu.3.sim.df,mindiff.date.1$date)
  
  name.clu3 = paste0("clu3.",1:n.group)
  colnames(clu.3.sim.data) <- name.clu3
  clu.3.sim.data = data.frame(clu.3.sim.data)
  
  # IN cluster 4
  clu.4 = cluster.result[cluster.result$fem.4 == 4,]$state
  clu.4.raw = mindiff.state.1[,names(mindiff.state.1)%in%clu.4]
  clu.4.rawmat.1 = matrix(unlist(clu.4.raw),ncol = 16)
  clu.4.fd<-smooth.basis(argvals = mindiff.date.1$date, y = clu.4.rawmat.1, err1.basis)$fd
  clu.4.fpca = pca.fd(clu.4.fd,nharm=3) 
  clu.4.fpca.loading = clu.4.fpca$harmonics$coefs
  clu.4.fpca.lambda = clu.4.fpca$values[1:3]
  clu.4.mean.coef = clu.4.fpca$meanfd$coefs
  
  ## Simulate n curves:
  clu.4.score1 = rnorm(n.group,0,sqrt(clu.4.fpca.lambda[1]))
  clu.4.score2 = rnorm(n.group,0,sqrt(clu.4.fpca.lambda[2]))
  clu.4.score3 = rnorm(n.group,0,sqrt(clu.4.fpca.lambda[3]))
  clu.4.scores = rbind(clu.4.score1,clu.4.score2,clu.4.score3)
  clu.4.simcoeff = t(clu.4.scores)%*%t(clu.4.fpca.loading)+matrix(clu.4.mean.coef,nrow=n.group,ncol=length(clu.4.mean.coef),byrow=TRUE)
  clu.4.sim.df = clu.4.fd
  clu.4.sim.df$coefs = t(clu.4.simcoeff)
  clu.4.sim.data = predict(clu.4.sim.df,mindiff.date.1$date)
  
  name.clu4 = paste0("clu4.",1:n.group)
  colnames(clu.4.sim.data) <- name.clu4
  clu.4.sim.data = data.frame(clu.4.sim.data)
  
  simu.all = cbind(clu.1.sim.data,clu.2.sim.data,clu.3.sim.data,clu.4.sim.data) 
  simu.mat = as.matrix(simu.all)
  length = dim(simu.mat)[1] * dim(simu.mat)[2]
  noise = matrix(rnorm(length, mean = 0, sd = 0.5), dim(simu.mat)[1])
  simu.all = simu.all + noise
  
  
  if(plotting == T){
    rawmat.1 = matrix(unlist(simu.all[,1:n.group]),ncol = n.group)
    rawmat.2 = matrix(unlist(simu.all[,(n.group+1):(2*n.group)]),ncol = n.group)
    rawmat.3 = matrix(unlist(simu.all[,(2*n.group+1):(3*n.group)]),ncol = n.group)
    rawmat.4 = matrix(unlist(simu.all[,(3*n.group+1):(4*n.group)]),ncol = n.group)
    rawmat1.fd<-smooth.basis(argvals = mindiff.date.1$date, y = rawmat.1, err1.basis)$fd
    rawmat2.fd<-smooth.basis(argvals = mindiff.date.1$date, y = rawmat.2, err1.basis)$fd
    rawmat3.fd<-smooth.basis(argvals = mindiff.date.1$date, y = rawmat.3, err1.basis)$fd
    rawmat4.fd<-smooth.basis(argvals = mindiff.date.1$date, y = rawmat.4, err1.basis)$fd
    par(mfrow=c(2,2))
    plot(rawmat1.fd,ylim = c(0,10));plot(rawmat2.fd,ylim = c(0,10))
    plot(rawmat3.fd,ylim = c(0,10));plot(rawmat4.fd,ylim = c(0,10))
  }
  
  simu.all$t = mindiff.date.1$date
  
  return (simu.all)
}

try = generate.curve.real.noise(20,plotting = T)

############################ Cluster Number Detection Simulation   #######################################################
simunum.realraw.noise.20= data.frame(matrix(ncol = 4, nrow = 1))
colnames(simunum.realraw.noise.20) <- c("method","index","clu.num","n.group")
seed = 888
for (i in 1:52){
  set.seed(seed)
  df = generate.curve.real.noise(20)
  df = df[,-which(names(df) %in% c("t"))]
  clunum.result.raw = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(clunum.result.raw) <- c("method","index","clu.num","n.group")
  clunum.result.raw$method = "kmeans on raw"
  clunum.result.raw$index = "average silhouette width"
  pamk.best <- pamk(df,krange=2:6)
  clunum.result.raw$clu.num = pamk.best$nc
  clunum.result.raw$n.group = 20
  
  simunum.realraw.noise.20 = rbind(simunum.realraw.noise.20,clunum.result.raw)
  seed = seed+10000*i
}

seed = 8888
for (i in 1:42){
  set.seed(seed)
  df = generate.curve.real.noise(20)
  df = df[,-which(names(df) %in% c("t"))]
  clunum.result.raw = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(clunum.result.raw) <- c("method","index","clu.num","n.group")
  clunum.result.raw$method = "kmeans on raw"
  clunum.result.raw$index = "average silhouette width"
  pamk.best <- pamk(df,krange=2:6)
  clunum.result.raw$clu.num = pamk.best$nc
  clunum.result.raw$n.group = 20
  
  simunum.realraw.noise.20 = rbind(simunum.realraw.noise.20,clunum.result.raw)
  seed = seed+10000*i
}

seed = 666
for (i in 1:11){
  set.seed(seed)
  df = generate.curve.real.noise(20)
  df = df[,-which(names(df) %in% c("t"))]
  clunum.result.raw = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(clunum.result.raw) <- c("method","index","clu.num","n.group")
  clunum.result.raw$method = "kmeans on raw"
  clunum.result.raw$index = "average silhouette width"
  pamk.best <- pamk(df,krange=2:6)
  clunum.result.raw$clu.num = pamk.best$nc
  clunum.result.raw$n.group = 20
  
  simunum.realraw.noise.20 = rbind(simunum.realraw.noise.20,clunum.result.raw)
  seed = seed+10000*i
}

seed = 6666
for (i in 1:11){
  set.seed(seed)
  df = generate.curve.real.noise(20)
  df = df[,-which(names(df) %in% c("t"))]
  clunum.result.raw = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(clunum.result.raw) <- c("method","index","clu.num","n.group")
  clunum.result.raw$method = "kmeans on raw"
  clunum.result.raw$index = "average silhouette width"
  pamk.best <- pamk(df,krange=2:6)
  clunum.result.raw$clu.num = pamk.best$nc
  clunum.result.raw$n.group = 20
  
  simunum.realraw.noise.20 = rbind(simunum.realraw.noise.20,clunum.result.raw)
  seed = seed+10000*i
}

seed = 68
for (i in 1:41){
  set.seed(seed)
  df = generate.curve.real.noise(20)
  df = df[,-which(names(df) %in% c("t"))]
  clunum.result.raw = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(clunum.result.raw) <- c("method","index","clu.num","n.group")
  clunum.result.raw$method = "kmeans on raw"
  clunum.result.raw$index = "average silhouette width"
  pamk.best <- pamk(df,krange=2:6)
  clunum.result.raw$clu.num = pamk.best$nc
  clunum.result.raw$n.group = 20
  
  simunum.realraw.noise.20 = rbind(simunum.realraw.noise.20,clunum.result.raw)
  seed = seed+10000*i
}

seed = 86
for (i in 1:43){
  set.seed(seed)
  df = generate.curve.real.noise(20)
  df = df[,-which(names(df) %in% c("t"))]
  clunum.result.raw = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(clunum.result.raw) <- c("method","index","clu.num","n.group")
  clunum.result.raw$method = "kmeans on raw"
  clunum.result.raw$index = "average silhouette width"
  pamk.best <- pamk(df,krange=2:6)
  clunum.result.raw$clu.num = pamk.best$nc
  clunum.result.raw$n.group = 20
  
  simunum.realraw.noise.20 = rbind(simunum.realraw.noise.20,clunum.result.raw)
  seed = seed+10000*i
}
simunum.realraw.noise.20 = na.omit(simunum.realraw.noise.20)
write.csv(simunum.realraw.noise.20, file = "simunum_realraw_noise_20.csv", row.names = FALSE)

### Summary analysis of the cluster number detection.
clunum.sim.noise.20.raw = read.csv("simunum_realraw_noise_20.csv",header = TRUE)
clunum.sim.noise.20.raw.summary = clunum.sim.noise.20.raw%>%group_by(n.group, method,
                                                         index, clu.num)%>%summarize(count = n())
clunum.sim.noise.20.raw.summary


################### Cluster Result Validation on Raw Data #############################
get.clu.accuracy.raw = function(df.clu, n.group){
  # Hungarian Alg
  real = rep(c(1,2,3,4),each = n.group)
  a = minWeightBipartiteMatching(df.clu, real)
  df.clu.1 = a[df.clu]
  
  raw.accu = sum(real == df.clu.1)/(n.group*4) 
  return(raw.accu)
}

#### Start the Simulation
simuvid.realraw.noise.20= data.frame(matrix(ncol = 4, nrow = 1))
colnames(simuvid.realraw.noise.20) <- c("method","index","n.group","accuracy")
seed = 888
for (i in  1:200){
  set.seed(seed)
  df = generate.curve.real.noise(20)
  df = df[,-which(names(df) %in% c("t"))]
  cluvid.raw.result = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(cluvid.raw.result) <- c("method","index","n.group","accuracy")
  cluvid.raw.result$method = c("kmeans on raw")
  cluvid.raw.result$index = c("average silhouette width")
  cluvid.raw.result$n.group = 20
  df.clu <- kmeans(t(as.matrix(df)),4)
  cluvid.raw.clu = df.clu$cluster
  
  cluvid.raw.result$accuracy = get.clu.accuracy.raw(cluvid.raw.clu, 20)
  
  simuvid.realraw.noise.20 = rbind(simuvid.realraw.noise.20,cluvid.raw.result)
  seed = seed+10000*i
}
simuvid.realraw.noise.20 = na.omit(simuvid.realraw.noise.20)
write.csv(simuvid.realraw.noise.20, file = "simuvid_realraw_noise_20.csv", row.names = FALSE)

### Summary the validation result
simuvid.realraw.noise = read.csv("simuvid_realraw_noise_20.csv",header = TRUE)
simuvid.realraw.noise.summary = simuvid.realraw%>%
  group_by(n.group, method,index)%>%
  summarize(mean.accuracy = mean(accuracy),
            sd.accuracy = sd(accuracy),
            se.accuracy = sd(accuracy)/sqrt(200))
simuvid.realraw.noise.summary
