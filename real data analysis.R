## This is the R file including all the content related to real data analysis.
## Please run all of it from the beginning to the end by order
## The main content of this file includes:
##     0. Installing required packages 
##     1. Data importing and dataset merging
##     2. Smooth FPCA on B-spline Non-parametric Regression
##     3. Clustering
##     4. Final Result Visualization
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
              "NbClust","funFEM","fda.usc","tclust","dtwclust","mapproj")
ipak(packages)
## The package fiftystater is for the U.S map visualization, and it may be unavailable on CRAN
## So its achieved file is attached in the Github repository with the path 


################################# Data importing, dataset merging and data preprocessing ##################
## 1. Importing each dataset
#Location Dataset
loc <- read.csv("locations.csv",header = T, stringsAsFactors = F)
loc$city_index <- as.numeric(rownames(loc))

# Historical Weather Dataset
hist <- read.csv("histWeather.csv",header = T, stringsAsFactors = F)
hist$Date = as.Date(hist$Date,format = "%Y-%m-%d")
hist <- unique(hist)
index=which.max(hist$Max_TemperatureF)
hist=hist[-c(index),]
 
# Forecast Weather Dataset
forecast <- read.table("forecast.dat", header=F, stringsAsFactors = F)
colnames(forecast) = c("city_index","forecast_date","forecast_value","forecast_goal",
                       "forecast_made_on")
forecast <- forecast %>% 
  mutate(
    forecast_date = as.Date(forecast_date,format = "%Y-%m-%d"),
    forecast_made_on = as.Date(forecast_made_on,format = "%Y-%m-%d"),
    # double check the meaning of "M" in forecast_value 
    forecast_value = as.numeric(forecast_value)
  )

## 2. Joining the dataset
# merge MinTemp with MaxTemp
fore_min_max <- full_join(
  x = forecast %>% filter(forecast_goal=="MinTemp") %>%  dplyr::select(-forecast_goal) %>% na.omit() %>% unique(),
  y = forecast %>% filter(forecast_goal=="MaxTemp") %>%  dplyr::select(-forecast_goal) %>% na.omit() %>% unique(),
  by = c("city_index","forecast_date","forecast_made_on"),
  suffix = c("_min", "_max")) %>% unique()

# merge MinTemp, MaxTemp, ProbPrecip and location to get fore_full_loc data 
fore_full_loc <- 
  # merge fore_min_max with ProbPrecip 
  full_join(
    x = fore_min_max,
    y = forecast %>% filter(forecast_goal=="ProbPrecip") %>% dplyr::select(-forecast_goal),
    by = c("city_index", "forecast_date", "forecast_made_on")
  ) %>% unique() %>% 
  dplyr::rename(  forecast_MinTemp = forecast_value_min,
                  forecast_MaxTemp = forecast_value_max,
                  forecast_ProbPre = forecast_value) %>%
  # merge fore_min_max_precip with location data
  left_join(y = loc, by = "city_index")

# master data with all infomation
df <- left_join( x = fore_full_loc, y = hist, 
                 by = c("forecast_date"= "Date", "AirPtCd" = "AirPtCd"))

# Daily Forecast Final Dataset
df.min.max <- df %>% dplyr::group_by(city_index,forecast_date,forecast_made_on,city,state,AirPtCd) %>% 
  dplyr::summarise(forecast_MinTemp = mean(forecast_MinTemp, na.rm=T),
                   Min_TemperatureF = mean(Min_TemperatureF, na.rm=T),
                   Mean_TemperatureF = mean(Mean_TemperatureF, na.rm=T),
                   forecast_MaxTemp = mean(forecast_MaxTemp, na.rm=T),
                   Max_TemperatureF = mean(Max_TemperatureF, na.rm=T)
  ) %>% ungroup() %>% 
  mutate(dayDiff = forecast_date - forecast_made_on,
         diffMin = forecast_MinTemp - Min_TemperatureF,
         diffMax = forecast_MaxTemp - Max_TemperatureF,
         year = lubridate::year(forecast_date),
         month = month(forecast_date),
         week = week(forecast_date),
         day = wday(forecast_date,label = T))

df.all <- left_join( x = df.min.max, y = hist %>% dplyr::select(-`Min_TemperatureF`,
                                                         Mean_TemperatureF,
                                                         Max_TemperatureF), 
                     by = c("forecast_date"= "Date", "AirPtCd" = "AirPtCd"))

## We load the final dataset to a RData for further convenience
save(df.all, file = "df.min.max.RData")


## 3. Data preprocessing
load("df.min.max.RData")
df.all = df.all %>% mutate(city = factor(city), state = factor(state))
df.all$city_state = paste(df.all$city,df.all$state)
region <- read.csv("region-to-state.new.csv",header = T, stringsAsFactors = F)

df.full = df.all %>% left_join(region, by="state") # number of [diffMin == NA's] = 89911

# Remove outliers
outlier = df.full %>% filter(diffMin< -50 | diffMin> 50) # 22 outliers
df = df.full %>% filter(diffMin>= -50)%>% filter(diffMin <= 50) # 726714 = 816647-89911-22

# Replace [Events == ""] to Normal and change the format of some variables
df$Events=as.character(df$Events)
df$Events[df$Events==""]<-"Normal"
df$CloudCover = factor(df$CloudCover)
df$PrecipitationIn = as.numeric(df$PrecipitationIn)
df$season = 4
df[df$month<12 & df$month>8,]$season = 3
df[df$month<9 & df$month>5,]$season = 2
df[df$month<6 & df$month>2,]$season = 1
df$season = factor(df$season)
df$city_state = paste(df$city,df$state,sep = "_") ## Some states share some same city names (3 cities have this problem)

# Filtering the data related to the weather prediction with day difference  = 1
df$diffMin = abs(df$diffMin)
df$state = tolower(df$state)
state.mindiff.df.1=df%>%filter(dayDiff == 1) %>% 
  dplyr::select(diffMin,forecast_date,state) %>% 
  dplyr::group_by(forecast_date, state) %>% 
  dplyr::summarise(avg=mean(diffMin,na.rm=T))

# Filling the missing value using weekly average for every state
try.1=state.mindiff.df.1 %>% spread(state, avg, fill=NA)

try.1$week=week(try.1$forecast_date)
try.1$year=lubridate::year(try.1$forecast_date)
try_new.1 <- try.1
try_new.1$forecast_date <- NULL

week.mean.err0.1 <- try_new.1 %>% dplyr::group_by(year,week) %>% 
  dplyr::summarise_all(dplyr::funs(mean(.,na.rm = T)))

for (row in 1: nrow(try.1)) {
  for (col in 1: ncol(try.1)) {
    if(is.na(try.1[row,col])) {
      week=as.numeric(try.1[row,"week"])
      year=as.numeric(try.1[row,"year"])
      try.1[row,col] = as.numeric(week.mean.err0.1[week.mean.err0.1$week == week & week.mean.err0.1$year == year,col+1])
    } 
  }
}
sum(is.na(try.1)) # summation  = 0, no missing value now

mindiff.state.1 = try.1[,-which(names(try.1) %in% c("forecast_date","week","year"))] 
names(mindiff.state.1) = tolower(names(mindiff.state.1))
mindiff.date.1 = try.1$forecast_date
mindiff.state.mat.1 = matrix(unlist(mindiff.state.1),ncol = 50)

write.csv(mindiff.date.1,"mindiff_day1_nona.csv", row.names = F)
write.csv(mindiff.state.1,"mindiff_state1_nona.csv", row.names = F)
################### Smooth FPCA on B-spline Non-parametric Regression ######################

## 1. Fitting Longitudinal Data to Functional Data Using B-spline
# Generalized Cross-validation fo Knots Selection
TableGCV<-c()
for (j in 13:37) { # the interval width is select from every 4 months to every 1 month
  bks.cv = as.Date(quantile(unclass(mindiff.date.1), seq(0,1,length = j)), origin = "1970-01-01")
  err.basis.cv = create.bspline.basis(rangeval = c(min(mindiff.date.1),max(mindiff.date.1)), breaks = bks.cv, norder = 4) 
  err.fit.cv = smooth.basis(argvals = mindiff.date.1, y = mindiff.state.mat.1, err.basis.cv)
  gcv.max <- max(err.fit.cv$gcv) # under each basis, record the maximum of GCV in 50 fitted curves for all states
  
  gcv.mu <- mean(err.fit.cv$gcv) # record the average GCV in 50 fitted curves for all states
  gcv.med <- median(err.fit.cv$gcv) # record the median of GCV in 50 fitted curves for all states
  TableGCV <- rbind(TableGCV,c(j,gcv.max,gcv.mu,gcv.med))
}

plot(TableGCV[,1], TableGCV[,2], type = "l")
plot(TableGCV[,1], TableGCV[,3], type = "l")
plot(TableGCV[,1], TableGCV[,4], type = "l")

bks.num.tab = TableGCV[order(TableGCV[,2]),]
bks.num.tab[1:5,] 

bks.num.tab = TableGCV[order(TableGCV[,3]),]
bks.num.tab[1:5,] 

bks.num.tab = TableGCV[order(TableGCV[,4]),]
bks.num.tab[1:5,] 
## All three plots shows the first local minimum at around 15, 
## so we finally select 15 seperating points (2 boundary and 13 interiors) 
## producing 14 intervals with same data amount.

# Fitting and Visualize the Function Object
bks.1 = as.Date(quantile(unclass(mindiff.date.1), seq(0,1,length = 15)), origin = "1970-01-01")
err1.basis = create.bspline.basis(rangeval = c(min(mindiff.date.1),max(mindiff.date.1)), breaks = bks.1, norder = 4) 
err1.fd<-smooth.basis(argvals = mindiff.date.1, y = mindiff.state.mat.1, err1.basis)$fd
plot(err1.fd,ylab = "Absolute Prediction Error (F)", xlab = "date")


## 2. Using One-curve-out method to select the penalty size of the PCA Smoothing
# Step 1, determine the number of PC Functions through Non-smoothing PCA 
# (at least 90% of Variation should be explained in the raw PC functions)
raw.fpca = pca.fd(err1.fd,nharm=5)
sum(raw.fpca$varprop) # Reach 90% Variation Excatly with 5 PC functions

# Step 2, Tuning Lambda
# Function for calculating the CV error for a given lambda 
get.mcv.smoothfpca = function(lambda){
  range = c(as.numeric(min(mindiff.date.1)),as.numeric(max(mindiff.date.1)))
  bspFpcPen = fdPar(fdobj=err1.basis, Lfdobj=2, lambda = lambda)
  smooth.fpc.cv = c()
  for (i in 1:50){
    smooth.fpc.raw = err1.fd
    smooth.fpc.raw$coefs = err1.fd$coefs[,-i]
    smooth.fpc = pca.fd(smooth.fpc.raw,nharm=5,harmfdPar=bspFpcPen)
    smooth.fpc.coefs = as.matrix(smooth.fpc$harmonics$coefs) #17*5
    smooth.fpc.score = as.matrix(smooth.fpc$scores) #49*5
    
    smooth.fpc.xi.fd = err1.fd
    smooth.fpc.xi.fd$coefs = err1.fd$coefs[,i]
    smooth.fpc.xi = eval.fd(mindiff.date.1,smooth.fpc.xi.fd) - eval.fd(mindiff.date.1,smooth.fpc$meanfd) #1033*1
    smooth.fpc.basismat = eval.fd(mindiff.date.1,smooth.fpc$harmonics) #1033*5
    
    # estimate the SSE = y'y - y'X(X'X)X'y
    smooth.fpc.ri = t(smooth.fpc.xi)%*%smooth.fpc.xi - t(smooth.fpc.xi)%*%
      smooth.fpc.basismat%*%ginv(t(smooth.fpc.basismat)%*%smooth.fpc.basismat)%*%
      t(smooth.fpc.basismat)%*%smooth.fpc.xi
    smooth.fpc.cv = c(smooth.fpc.cv,smooth.fpc.ri)
  }
  return(sum(smooth.fpc.cv))
}

# Choosing lambda from (0:500)*5e2 and visualize the result
table.lambda.mcv<-c()
for (j in 1:501) { 
  smooth.fpc.cv  =  get.mcv.smoothfpca((j-1)*500)
  table.lambda.mcv <- rbind(table.lambda.mcv,c(j,smooth.fpc.cv))
}

table.lambda.mcv = table.lambda.mcv[order(table.lambda.mcv[,1]),]
plot(table.lambda.mcv[1:100,1]*5e2, table.lambda.mcv[1:100,2], type = "l", 
     main = "CV plot Lambda Selection", 
     xlab = "Lambda Value", ylab = "CV")

table.lambda.mcv[,1] = table.lambda.mcv[,1]*500
table.lambda.mcv = table.lambda.mcv[order(table.lambda.mcv[,2]),]
table.lambda.mcv[1:5,] 
# We can see that the local minimum of CV error is around the lambda = 18500


################################# Clustering ####################################
### Applying the Smooth FPCA on the Time series data with the tuning parameter lambda = 18500
### However, to visualize the FPC functions with more smoothness to show the main pattern, we use 2,000,000
### The Clustering result are exactly the same between lambda = 18500 and lambda = 2,000,000
range = c(as.numeric(min(mindiff.date.1)),as.numeric(max(mindiff.date.1)))
bspFpcPen = fdPar(fdobj=err1.basis, Lfdobj=2, lambda = 2e6) ## You can change to 18500
                                                            ## The following clustering results are the same

smooth.fpc = pca.fd(err1.fd,nharm=5,harmfdPar=bspFpcPen)
smooth.fpc$varprop
sum(smooth.fpc$varprop)

# PC Functions Visualization
op <- par(mfrow=c(2,3))
plot(smooth.fpc$harmonics[1],main = paste("Smoothed PC1 function"),
     ylab = "PC 1 Function Value", xlab = "date")
plot(smooth.fpc$harmonics[2],main = paste("Smoothed PC2 function"),
     ylab = "PC 2 Function Value", xlab = "date")
plot(smooth.fpc$harmonics[3],main = paste("Smoothed PC3 function"),
     ylab = "PC 3 Function Value", xlab = "date")
plot(smooth.fpc$harmonics[4],main = paste("Smoothed PC4 function"),
     ylab = "PC 4 Function Value", xlab = "date")
plot(smooth.fpc$harmonics[5],main = paste("Smoothed PC5 function"),
     ylab = "PC 5 Function Value", xlab = "date")
plot(smooth.fpc$meanfd,main = "Mean Function",
     ylab = "Absolute Prediction Error (F)", xlab = "date")

############ Clustering Goal 1: Select Number of Clusters #############
####### Method 1: K-means on Smoothed FPC Scores
fpca.smooth <- data.frame(state = colnames(mindiff.state.1),
                          pc1 = smooth.fpc$scores[,1],
                          pc2 = smooth.fpc$scores[,2],
                          pc3 = smooth.fpc$scores[,3],
                          pc4 = smooth.fpc$scores[,4],
                          pc5 = smooth.fpc$scores[,5],
                          var.explain = sum(smooth.fpc$varprop),
                          daydiff = 1)

fpca.score = fpca.smooth%>%dplyr::select(pc1,pc2,pc3,pc4,pc5)
# Select number of cluster
set.seed(888)
fpc.kmeans.num <- NbClust(fpca.score, distance = "euclidean", 
                          min.nc=2, max.nc=8, 
                          method = "kmeans", 
                          index = "all")
max(fpc.kmeans.num$Best.partition) ## suggest 2 

####### Method 2: K-means on B-spline COefficients
err1.fd.coefs = t(err1.fd$coefs)
# Select number of cluster
set.seed(888)
bsp.kmeans.num <- NbClust(err1.fd.coefs, distance = "euclidean", 
                          min.nc=2, max.nc=8, 
                          method = "kmeans", 
                          index = "all")
max(bsp.kmeans.num$Best.partition) ## suggest 2 

####### Method 3: FunFEM Clustering Using BIC Criterions
set.seed(888)
fem.bic = funFEM(err1.fd,K = 2:8,model = "all", crit = "bic", init = "kmeans",eps = 1e-5)
fem.bic$K ##Suggest 4

####### Method 4: FunFEM Clustering Using ICL Criterions
set.seed(888)
fem.icl = funFEM(err1.fd,K = 2:8,model = "all", crit = "icl", init = "kmeans",eps = 1e-5)
fem.icl$K ##Suggest 4

### According to the Simulation Study, FunFEM has outstanding performance compared to K-means on Cluster Number Selection
### So we believe the Cluster Number Selection Result in FunFEM, and choose cluster number K = 4



############### Clustering Goal 2: Final Clustering Result with Cluster Number = 4
####### Method 1: K-means on Smoothed FPC Scores 
# Tuning Start Point: 8,88,888,8888,88888,888888,8888888,88888888
sd.fpc.clu4  <- matrix(NA, nrow=8, ncol=2)
rr=1
seed.candidate = c(8,88,888,8888,88888,888888,8888888,88888888)
sd.fpc.clu4[,1] = seed.candidate
for (rr in 1:8){
  set.seed(seed.candidate[rr])
  
  fpca.smooth.temp = fpca.smooth 
  fpc4.temp <- kmeans(fpca.score , 4) # 4 cluster solution
  fpca.smooth.temp$fpca4 <- fpc4.temp$cluster
  
  ## 4 clusters
  # 1
  kfpc4.state.1.temp = fpca.smooth.temp[fpca.smooth.temp$fpca4 == 1,]$state
  kfpc4.state.1.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.1.temp]),ncol= table(fpca.smooth.temp$fpca4)[1])
  kfpc4.state.1.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.1.ts.temp, err1.basis)$fd
  
  # 2
  kfpc4.state.2.temp = fpca.smooth.temp[fpca.smooth.temp$fpca4 == 2,]$state
  kfpc4.state.2.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.2.temp]),ncol= table(fpca.smooth.temp$fpca4)[2])
  kfpc4.state.2.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.2.ts.temp, err1.basis)$fd
  
  # 3
  kfpc4.state.3.temp = fpca.smooth.temp[fpca.smooth.temp$fpca4 == 3,]$state
  kfpc4.state.3.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.3.temp]),ncol= table(fpca.smooth.temp$fpca4)[3])
  kfpc4.state.3.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.3.ts.temp, err1.basis)$fd
  
  # 4
  kfpc4.state.4.temp = fpca.smooth.temp[fpca.smooth.temp$fpca4 == 4,]$state
  kfpc4.state.4.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.4.temp]),ncol= table(fpca.smooth.temp$fpca4)[4])
  kfpc4.state.4.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.4.ts.temp, err1.basis)$fd
  
  sd.all.temp = c(int.simpson(sd.fd(kfpc4.state.1.fd.temp),method="TRAPZ"),
                  int.simpson(sd.fd(kfpc4.state.2.fd.temp),method="TRAPZ"),
                  int.simpson(sd.fd(kfpc4.state.3.fd.temp),method="TRAPZ"),
                  int.simpson(sd.fd(kfpc4.state.4.fd.temp),method="TRAPZ"))
  sd.fpc.clu4[rr,2] <- mean(sd.all.temp) 
}

best.fpc.seed<- sd.fpc.clu4[which.min(sd.fpc.clu4[,2]),1] ## Choose 8
set.seed(best.fpc.seed)
kmeans.fpca.4 <- kmeans(fpca.score, 4) # Final 4 cluster solution
# Final clustering Result 
fpca.smooth <- data.frame(fpca.smooth, fpc.4.large = kmeans.fpca.4$cluster)

####### Method 2: K-means on B-spline Coefficients with Tuning Start Point 8,88,888,8888,88888,888888,8888888,88888888
sd.bsp.clu4 <- matrix(NA, nrow=8, ncol=2)
rr=1
seed.candidate = c(8,88,888,8888,88888,888888,8888888,88888888)
sd.bsp.clu4[,1] = seed.candidate
for (rr in 1:8){
  set.seed(seed.candidate[rr])
  
  fpca.ts.1.temp = fpca.smooth
  fpca.temp <- kmeans(err1.fd.coefs, 4) # 4 cluster solution
  fpca.ts.1.temp$bsp4 <- fpca.temp$cluster
  
  ## 4 clusters
  # 1
  kfpc4.state.1.temp = fpca.ts.1.temp[fpca.ts.1.temp$bsp4 == 1,]$state
  kfpc4.state.1.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.1.temp]),ncol= table(fpca.ts.1.temp$bsp4)[1])
  kfpc4.state.1.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.1.ts.temp, err1.basis)$fd
  
  # 2
  kfpc4.state.2.temp = fpca.ts.1.temp[fpca.ts.1.temp$bsp4 == 2,]$state
  kfpc4.state.2.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.2.temp]),ncol= table(fpca.ts.1.temp$bsp4)[2])
  kfpc4.state.2.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.2.ts.temp, err1.basis)$fd
  
  # 3
  kfpc4.state.3.temp = fpca.ts.1.temp[fpca.ts.1.temp$bsp4 == 3,]$state
  kfpc4.state.3.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.3.temp]),ncol= table(fpca.ts.1.temp$bsp4)[3])
  kfpc4.state.3.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.3.ts.temp, err1.basis)$fd
  
  # 4
  kfpc4.state.4.temp = fpca.ts.1.temp[fpca.ts.1.temp$bsp4 == 4,]$state
  kfpc4.state.4.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.4.temp]),ncol= table(fpca.ts.1.temp$bsp4)[4])
  kfpc4.state.4.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.4.ts.temp, err1.basis)$fd
  
  sd.all.temp = c(int.simpson(sd.fd(kfpc4.state.1.fd.temp),method="TRAPZ"),
                  int.simpson(sd.fd(kfpc4.state.2.fd.temp),method="TRAPZ"),
                  int.simpson(sd.fd(kfpc4.state.3.fd.temp),method="TRAPZ"),
                  int.simpson(sd.fd(kfpc4.state.4.fd.temp),method="TRAPZ"))
  sd.bsp.clu4[rr,2] <- mean(sd.all.temp) 
}

best.bsp.seed<- sd.bsp.clu4[which.min(sd.bsp.clu4[,2]),1] ## Choose 88888
set.seed(best.bsp.seed)
kmeans.bsp.4 <- kmeans(err1.fd.coefs, 4) # Final 4 cluster solution
# Final clustering Result 
fpca.smooth <- data.frame(fpca.smooth, bsp.4 = kmeans.bsp.4$cluster)

####### Method 3: FunFEM Clustering Using BIC Criterions with Tuning Start Point 8,88,888,8888,88888,888888,8888888,88888888
sd.fem.clu4 <- matrix(NA, nrow=8, ncol=2)
rr=1
seed.candidate = c(8,88,888,8888,88888,888888,8888888,88888888)
sd.fem.clu4[,1] = seed.candidate
for (rr in 1:8){
  set.seed(seed.candidate[rr])
  
  fpca.ts.1.temp = fpca.smooth
  fpca.temp <- funFEM(err1.fd,K = 4,model = "all",crit = "bic", init = "kmeans",eps = 1e-5) # 4 cluster solution
  fpca.ts.1.temp$fem4 <- fpca.temp$cls
  
  ## 4 clusters
  # 1
  kfpc4.state.1.temp = fpca.ts.1.temp[fpca.ts.1.temp$fem4 == 1,]$state
  kfpc4.state.1.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.1.temp]),ncol= table(fpca.ts.1.temp$fem4)[1])
  kfpc4.state.1.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.1.ts.temp, err1.basis)$fd
  
  # 2
  kfpc4.state.2.temp = fpca.ts.1.temp[fpca.ts.1.temp$fem4 == 2,]$state
  kfpc4.state.2.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.2.temp]),ncol= table(fpca.ts.1.temp$fem4)[2])
  kfpc4.state.2.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.2.ts.temp, err1.basis)$fd
  
  # 3
  kfpc4.state.3.temp = fpca.ts.1.temp[fpca.ts.1.temp$fem4 == 3,]$state
  kfpc4.state.3.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.3.temp]),ncol= table(fpca.ts.1.temp$fem4)[3])
  kfpc4.state.3.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.3.ts.temp, err1.basis)$fd
  
  # 4
  kfpc4.state.4.temp = fpca.ts.1.temp[fpca.ts.1.temp$fem4 == 4,]$state
  kfpc4.state.4.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.4.temp]),ncol= table(fpca.ts.1.temp$fem4)[4])
  kfpc4.state.4.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.4.ts.temp, err1.basis)$fd
  
  sd.all.temp = c(int.simpson(sd.fd(kfpc4.state.1.fd.temp),method="TRAPZ"),
                  int.simpson(sd.fd(kfpc4.state.2.fd.temp),method="TRAPZ"),
                  int.simpson(sd.fd(kfpc4.state.3.fd.temp),method="TRAPZ"),
                  int.simpson(sd.fd(kfpc4.state.4.fd.temp),method="TRAPZ"))
  sd.fem.clu4[rr,2] <- mean(sd.all.temp) 
}

best.fem.seed<- sd.fem.clu4[which.min(sd.fem.clu4[,2]),1] ## Choose 88888
set.seed(best.fem.seed)
kmeans.fem.4 <- funFEM(err1.fd,K = 4,model = "all",crit = "bic", init = "kmeans",eps = 1e-5) # Final 4 cluster solution
# 4 cluster solution
fpca.smooth <- data.frame(fpca.smooth, fem.4 = kmeans.fem.4$cls)

####### Method 4: FunFEM Clustering Using ICL Criterions with Tuning Start Point 8,88,888,8888,88888,888888,8888888,88888888
sd.fem.clu4 <- matrix(NA, nrow=8, ncol=2)
rr=1
seed.candidate = c(8,88,888,8888,88888,888888,8888888,88888888)
sd.fem.clu4[,1] = seed.candidate
for (rr in 1:8){
  set.seed(seed.candidate[rr])
  
  fpca.ts.1.temp = fpca.smooth
  fpca.temp <- funFEM(err1.fd,K = 4,model = "all",crit = "icl", init = "kmeans",eps = 1e-5) # 4 cluster solution
  fpca.ts.1.temp$fem4 <- fpca.temp$cls
  
  ## 4 clusters
  # 1
  kfpc4.state.1.temp = fpca.ts.1.temp[fpca.ts.1.temp$fem4 == 1,]$state
  kfpc4.state.1.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.1.temp]),ncol= table(fpca.ts.1.temp$fem4)[1])
  kfpc4.state.1.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.1.ts.temp, err1.basis)$fd
  
  # 2
  kfpc4.state.2.temp = fpca.ts.1.temp[fpca.ts.1.temp$fem4 == 2,]$state
  kfpc4.state.2.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.2.temp]),ncol= table(fpca.ts.1.temp$fem4)[2])
  kfpc4.state.2.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.2.ts.temp, err1.basis)$fd
  
  # 3
  kfpc4.state.3.temp = fpca.ts.1.temp[fpca.ts.1.temp$fem4 == 3,]$state
  kfpc4.state.3.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.3.temp]),ncol= table(fpca.ts.1.temp$fem4)[3])
  kfpc4.state.3.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.3.ts.temp, err1.basis)$fd
  
  # 4
  kfpc4.state.4.temp = fpca.ts.1.temp[fpca.ts.1.temp$fem4 == 4,]$state
  kfpc4.state.4.ts.temp = matrix(unlist(mindiff.state.1[,names(mindiff.state.1)%in%kfpc4.state.4.temp]),ncol= table(fpca.ts.1.temp$fem4)[4])
  kfpc4.state.4.fd.temp = smooth.basis(argvals = mindiff.date.1, y = kfpc4.state.4.ts.temp, err1.basis)$fd
  
  sd.all.temp = c(int.simpson(sd.fd(kfpc4.state.1.fd.temp),method="TRAPZ"),
                  int.simpson(sd.fd(kfpc4.state.2.fd.temp),method="TRAPZ"),
                  int.simpson(sd.fd(kfpc4.state.3.fd.temp),method="TRAPZ"),
                  int.simpson(sd.fd(kfpc4.state.4.fd.temp),method="TRAPZ"))
  sd.fem.clu4[rr,2] <- mean(sd.all.temp) 
}

best.fem.seed<- sd.fem.clu4[which.min(sd.fem.clu4[,2]),1] ## Choose 88888
set.seed(best.fem.seed)
kmeans.fem.4 <- funFEM(err1.fd,K = 4,model = "all",crit = "icl", init = "kmeans",eps = 1e-5) # Final 4 cluster solution
# Final clustering Result 
fpca.smooth <- data.frame(fpca.smooth, fem.4.icl = kmeans.fem.4$cls)


####### Matching the Cluster Results Among 4 Clustering Methods Using Hungarian Algorithm For Cluster Mapping
# labels from cluster A will be matched on the labels from cluster B
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

minWeightBipartiteMatching(fpca.smooth$bsp.4,fpca.smooth$fpc.4.large)
minWeightBipartiteMatching(fpca.smooth$bsp.4,fpca.smooth$fem.4)
minWeightBipartiteMatching(fpca.smooth$fem.4,fpca.smooth$fem.4.icl)

############################## Final Result Visualization ###################################
######## Map Visualization
region <- read.csv("region-to-state.new.csv",header = T, stringsAsFactors = F)
region$state = tolower(region$state)
region_noremoted = region[!(region$state %in% c("alaska","hawaii")),]

# K-means Clustering on Smoothed FPC Scores with Start Point Tuning
fpca.smooth$state = tolower(fpca.smooth$state)
fpc4.smooth.tuned.large = ggplot(fpca.smooth, aes(map_id = state)) +
  geom_map(aes(fill = factor(fpc.4.large)),map = fifty_states) +
  expand_limits(x = fifty_states$long, y = fifty_states$lat) +
  coord_map()  +
  labs(x = "",y = "",  fill = "Cluster")+
  geom_text(data = region_noremoted, aes(long,lag,label = state.short), size=3)+ 
  fifty_states_inset_boxes()+ borders("state")
# If there is an error when visualizing the plot, please run
#dev.off() 
fpc4.smooth.tuned.large

# B-spline Coefficients with Start Point Tuning
fpca.smooth$state = tolower(fpca.smooth$state)
bsp4.tuned = ggplot(fpca.smooth, aes(map_id = state)) +
  geom_map(aes(fill = factor(bsp.4)),map = fifty_states) +
  expand_limits(x = fifty_states$long, y = fifty_states$lat) +
  coord_map()  +
  labs(x = "",y = "",  fill = "Cluster")+
  geom_text(data = region_noremoted, aes(long,lag,label = state.short), size=3)+ 
  fifty_states_inset_boxes()+ borders("state")
# If there is an error when visualizing the plot, please run
#dev.off() 
bsp4.tuned

# FunFEM BIC
fpca.smooth$state = tolower(fpca.smooth$state)
fem4.tuned = ggplot(fpca.smooth, aes(map_id = state)) +
  geom_map(aes(fill = factor(fem.4)),map = fifty_states) +
  expand_limits(x = fifty_states$long, y = fifty_states$lat) +
  coord_map()  +
  labs(x = "",y = "",  fill = "Cluster")+
  geom_text(data = region_noremoted, aes(long,lag,label = state.short), size=3)+ 
  fifty_states_inset_boxes()+ borders("state")
# If there is an error when visualizing the plot, please run
#dev.off() 
fem4.tuned

# FunFEM ICL
fpca.smooth$state = tolower(fpca.smooth$state)
fem4.tuned.icl = ggplot(fpca.smooth, aes(map_id = state)) +
  geom_map(aes(fill = factor(fem.4.icl)),map = fifty_states) +
  expand_limits(x = fifty_states$long, y = fifty_states$lat) +
  coord_map()  +
  labs(x = "",y = "",  fill = "Cluster")+
  geom_text(data = region_noremoted, aes(long,lag,label = state.short), size=3)+ 
  fifty_states_inset_boxes()+ borders("state")
# If there is an error when visualizing the plot, please run
#dev.off() 
fem4.tuned.icl

######## Curves Visualization in each Cluster with 95% Confidence Interval
plot.ci = function(result.df, original.ts, method){
  temp.df = result.df
  if (method == "fpc.large"){
    temp.df$temp.clu = temp.df$fpc.4.large
  }else if (method == "bsp"){
    temp.df$temp.clu = temp.df$bsp.4
  }else if (method == "fem"){
    temp.df$temp.clu = temp.df$fem.4
  }else{
    temp.df$temp.clu = temp.df$fem.4.icl
  }
  
  par(mfrow=c(2,4))
  
  ## 4 clusters
  temp.state.1 = temp.df[temp.df$temp.clu == 1,]$state
  temp.state.1.ts = matrix(unlist(original.ts[,names(original.ts)%in%temp.state.1]),ncol= table(temp.df$temp.clu)[1])
  temp.state.1.fd = smooth.basis(argvals = mindiff.date.1, y = temp.state.1.ts, err1.basis)$fd
  plot(temp.state.1.fd,ylim = c(1,12),ylab = "Absolute Prediction Error (F)", xlab = "date",cex.lab=1.65,cex.axis=1.4)
  
  temp.state.2 = temp.df[temp.df$temp.clu == 2,]$state
  temp.state.2.ts = matrix(unlist(original.ts[,names(original.ts)%in%temp.state.2]),ncol= table(temp.df$temp.clu)[2])
  temp.state.2.fd = smooth.basis(argvals = mindiff.date.1, y = temp.state.2.ts, err1.basis)$fd
  plot(temp.state.2.fd,ylim = c(1,12),ylab = "Absolute Prediction Error (F)", xlab = "date",cex.lab=1.65,cex.axis=1.4)
  
  temp.state.3 = temp.df[temp.df$temp.clu == 3,]$state
  temp.state.3.ts = matrix(unlist(original.ts[,names(original.ts)%in%temp.state.3]),ncol= table(temp.df$temp.clu)[3])
  temp.state.3.fd = smooth.basis(argvals = mindiff.date.1, y = temp.state.3.ts, err1.basis)$fd
  plot(temp.state.3.fd,ylim = c(1,12),ylab = "Absolute Prediction Error (F)", xlab = "date",cex.lab=1.65,cex.axis=1.4)
  
  temp.state.4 = temp.df[temp.df$temp.clu == 4,]$state
  temp.state.4.ts = matrix(unlist(original.ts[,names(original.ts)%in%temp.state.4]),ncol= table(temp.df$temp.clu)[4])
  temp.state.4.fd = smooth.basis(argvals = mindiff.date.1, y = temp.state.4.ts, err1.basis)$fd
  plot(temp.state.4.fd,ylim = c(1,12),ylab = "Absolute Prediction Error (F)", xlab = "date",cex.lab=1.65,cex.axis=1.4)
  
  # 95% CI
  muhat = mean.fd(temp.state.1.fd)
  sdhat = sd.fd(temp.state.1.fd)
  n = length(temp.state.1.fd)
  se_hat_u = fd(basisobj = err1.basis)
  se_hat_u$coefs = 1.96*sdhat$coefs/sqrt(n) + muhat$coefs
  se_hat_l = fd(basisobj = err1.basis)
  se_hat_l$coefs = -1.96*sdhat$coefs/sqrt(n) + muhat$coefs
  plot.fd(muhat,ylim = c(1,12),ylab = "Absolute Prediction Error (F)", xlab = "date",cex.lab=1.65,cex.axis=1.4)
  plot.fd(se_hat_l,add = TRUE,col = "red",lty = 2)
  plot.fd(se_hat_u,add = TRUE,col = "red",lty = 2)
  
  muhat = mean.fd(temp.state.2.fd)
  sdhat = sd.fd(temp.state.2.fd)
  n = length(temp.state.2.fd)
  se_hat_u = fd(basisobj = err1.basis)
  se_hat_u$coefs = 1.96*sdhat$coefs/sqrt(n) + muhat$coefs
  se_hat_l = fd(basisobj = err1.basis)
  se_hat_l$coefs = -1.96*sdhat$coefs/sqrt(n) + muhat$coefs
  plot.fd(muhat,ylim = c(1,12),ylab = "Absolute Prediction Error (F)", xlab = "date",cex.lab=1.65,cex.axis=1.4)
  plot.fd(se_hat_l,add = TRUE,col = "red",lty = 2)
  plot.fd(se_hat_u,add = TRUE,col = "red",lty = 2)
  
  muhat = mean.fd(temp.state.3.fd)
  sdhat = sd.fd(temp.state.3.fd)
  n = length(temp.state.3.fd)
  se_hat_u = fd(basisobj = err1.basis)
  se_hat_u$coefs = 1.96*sdhat$coefs/sqrt(n) + muhat$coefs
  se_hat_l = fd(basisobj = err1.basis)
  se_hat_l$coefs = -1.96*sdhat$coefs/sqrt(n) + muhat$coefs
  plot.fd(muhat, ylim = c(1,12),ylab = "Absolute Prediction Error (F)", xlab = "date",cex.lab=1.65,cex.axis=1.4)
  plot.fd(se_hat_l,add = TRUE,col = "red",lty = 2)
  plot.fd(se_hat_u,add = TRUE,col = "red",lty = 2)
  
  muhat = mean.fd(temp.state.4.fd)
  sdhat = sd.fd(temp.state.4.fd)
  n = length(temp.state.4.fd)
  se_hat_u = fd(basisobj = err1.basis)
  se_hat_u$coefs = 1.96*sdhat$coefs/sqrt(n) + muhat$coefs
  se_hat_l = fd(basisobj = err1.basis)
  se_hat_l$coefs = -1.96*sdhat$coefs/sqrt(n) + muhat$coefs
  plot.fd(muhat, ylim = c(1,12),ylab = "Absolute Prediction Error (F)", xlab = "date",cex.lab=1.65,cex.axis=1.4)
  plot.fd(se_hat_l,add = TRUE,col = "red",lty = 2)
  plot.fd(se_hat_u,add = TRUE,col = "red",lty = 2)
  
  sd.all = c(int.simpson(sd.fd(temp.state.1.fd),method="TRAPZ"),
             int.simpson(sd.fd(temp.state.2.fd), method="TRAPZ"),
             int.simpson(sd.fd(temp.state.3.fd), method="TRAPZ"),
             int.simpson(sd.fd(temp.state.4.fd), method="TRAPZ"))
  
  rank.df = data.frame(sd.mu = mean(sd.all),
                       mu.clu1 = int.simpson(mean.fd(temp.state.1.fd), method="TRAPZ"),
                       mu.clu2 = int.simpson(mean.fd(temp.state.2.fd), method="TRAPZ"),
                       mu.clu3 = int.simpson(mean.fd(temp.state.3.fd), method="TRAPZ"),
                       mu.clu4 = int.simpson(mean.fd(temp.state.4.fd), method="TRAPZ"))
  
  return(rank.df)
}

plot.ci(fpca.smooth, mindiff.state.1, "fpc.large")
plot.ci(fpca.smooth, mindiff.state.1, "bsp")
plot.ci(fpca.smooth, mindiff.state.1, "fem")
plot.ci(fpca.smooth, mindiff.state.1, "fem.icl")

write.csv(fpca.smooth, file = 'real_analysis_result.csv',row.names = F)
