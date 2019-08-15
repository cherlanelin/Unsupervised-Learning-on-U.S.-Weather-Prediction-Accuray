# Unsupervised-Learning-on-U.S.-Weather-Prediction-Accuray
Using functional data unsupervised learning techniques, FPCA and clustering, to evaluate the performance of U.S weather prediction.

#### This is the reproducible code package for the paper Unsupervised Learning on U.S. Weather Prediction Accuray. The repository includes the real and raw U.S. weather data, R code of real data, and R code for simulation study.

### Programming language of the project
R and Rstudio

### Real and Raw Datasets
1. forecast.zip

The file contains all the forecast records. Please get the forecast.dat file from it first.

2. location.csv

The file contains all the geoinformation of the weather forecast locations.

3. histWeather.csv

The file contains all the historical records of the real weather features measurements.

4. region-to-state.new.csv

The file contains all the short names of states and the geolocation of the center of the states.

5. Description.txt

This the introduction the contents of the real data files 1. forecast.zip, 2. location.csv, and 3. histWeather.csv.

### real_data_analysis.R
This is the R files including all the content related to real data analysis. It contains the following sections:

1. Data importing and dataset merging
2. Smooth FPCA on B-spline Non-parametric Regression
3. Clustering

#### The package "fiftystater" is used for clustering result visualization on U.S. map. It may be unavailable on CRAN, so I attached the achieved fold of the package, the "fiftystater_1.0.1.tar.gz" in the Github repository.
 

### simulation study analysis.R
This is the R files including all the content related to simulation study. The main content of this file includes:
1. Functions for Data simulation and Analysis
2. Smooth FPCA on B-spline Non-parametric Regression
3. Clustering Number Selection Validation
4. Clustering Validation Study

