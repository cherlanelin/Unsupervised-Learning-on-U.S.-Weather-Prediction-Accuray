## This is the R file including all the content related to the simulated data based on Null case clustering.
## Which means, detecting number of cluster
## Please run all of it from the beginning to the end by order
## The main content of this file includes:
##     0. Installing required packages 
##     1. Functions for Data simulation and Analysis
##     2. Clustering Number Selection Validation
######################################################################################


################################## Installing required packages ######################
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



# Common Noise Simulation
generate.curve.null<-function(t.min,t.max,t.by,n.group,err.sd){
  t = seq(t.min, t.max, by = t.by)
  
  mean.clu1 = sin(2*t) 
  
  simu.clu1 = data.frame(matrix(ncol = n.group, nrow = length(t)))
  name.clu1 = paste0("clu1.",1:n.group)
  colnames(simu.clu1) <- name.clu1
  
  for(column in name.clu1){
    simu.clu1[,column] = mean.clu1 + rnorm(length(t),mean = 0, sd = err.sd)
  }
  
  simu.all = simu.clu1
  simu.all$t = t
  
  return (simu.all)
}



# Value-Proportion Noise Simulation
generate.curve.ctvardiff.null<-function(t.min,t.max,t.by,n.group,err.sd){
  t = seq(t.min, t.max, by = t.by)
  
  mean.clu1 = sin(2*t) 
 
  simu.clu1 = data.frame(matrix(ncol = n.group, nrow = length(t)))
  name.clu1 = paste0("clu1.",1:n.group)
  colnames(simu.clu1) <- name.clu1
  
  for(column in name.clu1){
    simu.clu1[,column] = mean.clu1 + rnorm(length(t),mean = 0, sd = err.sd)*abs(sin(2*t))
  }
  
  simu.all = simu.clu1
  simu.all$t = t
  
  return (simu.all)
}

#Cluster Number Selection
cluster.select.num = function(simu.all,n.group,knots.num,nharm,plotting = FALSE){
  t = simu.all$t
  simu.all = simu.all[,-which(names(simu.all) %in% c("t"))]
  simu.all.mat = matrix(unlist(simu.all),ncol = n.group)
  
  bks = quantile(t, seq(0,1,length = knots.num+1))
  simu.basis = create.bspline.basis(rangeval = c(min(t),max(t)), breaks = bks, norder = 4) 
  simu.fd<-smooth.basis(argvals = t, y = simu.all.mat, simu.basis)$fd
  
  if(plotting == TRUE){plot(simu.fd)}
  
  #Kmeans on Raw Data
  
  
  #FEM
  fem.clu.bic = funFEM(simu.fd,K = 1:6,model = "all",crit = "bic", init = "kmeans",eps = 1e-4)
  fem.num.bic = fem.clu.bic$K
  
  fem.clu.icl = funFEM(simu.fd,K = 1:6,model = "all",crit = "icl", init = "kmeans",eps = 1e-4)
  fem.num.icl = fem.clu.icl$K
  
  
  # kmeans on coef
  fd.coefs = t(simu.fd$coefs)
  #average silhouette width  
  coeff.kmeans.best <- pamk(fd.coefs,krange=1:6)
  bsp.clu.select =  coeff.kmeans.best$nc

  # kmeans on smoothed FPC scores
  range = c(as.numeric(min(t)),as.numeric(max(t)))
  bspFpcPen = fdPar(fdobj=simu.basis, Lfdobj=2, lambda = 10^(-2))
  smooth.fpc = pca.fd(simu.fd,nharm=nharm,harmfdPar=bspFpcPen)
  fd.scores = smooth.fpc$scores
  
  #all index  
  fpc.kmeans.best <- pamk(fd.scores,krange=1:6)
  fpc.clu.select =  fpc.kmeans.best$nc
  
  simu.clunum.compare = rbind(fem.num.bic,fem.num.icl,
                              bsp.clu.select,fpc.clu.select)
  return(simu.clunum.compare)
}


#################Clustering Number Selection Validation Simulation ##########################

# Round 1, common noise
simunum.null.comnoise.40 = data.frame(matrix(ncol = 5, nrow = 1))
colnames(simunum.null.comnoise.40) <- c("method","index","clu.num","n.group","round")
seed = 888
for (i in  1:98){
  set.seed(seed)
  df = generate.curve.null(t.min=0,t.max=10,t.by=0.01,
                           n.group=40,err.sd=0.5)
  cluster.result = data.frame(matrix(ncol = 5, nrow = 4))
  colnames(cluster.result) <- c("method","index","clu.num","n.group","round")
  cluster.result$method = c("fem","fem","bsp","fpc")
  cluster.result$index = c("bic","icl","average silhouette width","average silhouette width")
  cluster.result$n.group = rep(40,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 0.5, common noise",each = 4)
  cluster.result$clu.num = cluster.select.num(simu.all=df,n.group=40,knots.num=10,nharm = 4,plotting = F)
  
  simunum.null.comnoise.40 = rbind(simunum.null.comnoise.40,cluster.result)
  seed = seed+10000*i
}

seed = 8888
for (i in  1:2){
  set.seed(seed)
  df = generate.curve.null(t.min=0,t.max=10,t.by=0.01,
                           n.group=40,err.sd=0.5)
  cluster.result = data.frame(matrix(ncol = 5, nrow = 4))
  colnames(cluster.result) <- c("method","index","clu.num","n.group","round")
  cluster.result$method = c("fem","fem","bsp","fpc")
  cluster.result$index = c("bic","icl","average silhouette width","average silhouette width")
  cluster.result$n.group = rep(40,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 0.5, common noise",each = 4)
  cluster.result$clu.num = cluster.select.num(simu.all=df,n.group=40,knots.num=10,nharm = 3,plotting = F)
  
  simunum.null.comnoise.40 = rbind(simunum.null.comnoise.40,cluster.result)
  seed = seed+10000*i
}
simunum.null.comnoise.40 = na.omit(simunum.null.comnoise.40)
write.csv(simunum.null.comnoise.40, file = "simunum_null_comnoise_40.csv", row.names = FALSE)

# Round 2, varied noise
simunum.null.varnoise.40 = data.frame(matrix(ncol = 5, nrow = 1))
colnames(simunum.null.varnoise.40) <- c("method","index","clu.num","n.group","round")
seed = 888
for (i in  1:100){
  set.seed(seed)
  df = generate.curve.ctvardiff.null(t.min=0,t.max=10,t.by=0.01,
                           n.group=40,err.sd=0.5)
  cluster.result = data.frame(matrix(ncol = 5, nrow = 4))
  colnames(cluster.result) <- c("method","index","clu.num","n.group","round")
  cluster.result$method = c("fem","fem","bsp","fpc")
  cluster.result$index = c("bic","icl","average silhouette width","average silhouette width")
  cluster.result$n.group = rep(40,each = 4)
  cluster.result$round = rep("by = 0.01, sd = 0.5, varied noise",each = 4)
  cluster.result$clu.num = cluster.select.num(simu.all=df,n.group=40,knots.num=10,nharm = 3,plotting = F)
  
  simunum.null.varnoise.40 = rbind(simunum.null.varnoise.40,cluster.result)
  seed = seed+10000*i
}
simunum.null.varnoise.40 = na.omit(simunum.null.varnoise.40)
write.csv(simunum.null.varnoise.40, file = "simunum_null_varnoise_40.csv", row.names = FALSE)

### Summary analysis of the cluster number detection.
simunum.null.comnoise.40.df = read.csv("simunum_null_comnoise_40.csv",header = TRUE)
simunum.null.varnoise.40.df = read.csv("simunum_null_varnoise_40.csv",header = TRUE)


#common noise
simunum.null.comnoise.40.summary = simunum.null.comnoise.40.df%>%group_by(round, n.group, method,
                             index, clu.num)%>%summarize(count = n())
simunum.null.comnoise.40.summary

#varied noise
simunum.null.varnoise.40.summary = simunum.null.varnoise.40.df%>%group_by(round, n.group, method,
                             index, clu.num)%>%summarize(count = n())
simunum.null.varnoise.40.summary
