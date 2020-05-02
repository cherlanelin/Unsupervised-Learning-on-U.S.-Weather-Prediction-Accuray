# Unsupervised-Learning-on-U.S.-Weather-Prediction-Accuray
Using functional data unsupervised learning techniques, FPCA and clustering, to evaluate the performance of U.S weather prediction.

#### This is the reproducible code package for the paper Unsupervised Learning on U.S. Weather Prediction Accuray. The repository includes the real and raw U.S. weather data, R code of real data, and R code for simulation study.

### Programming language of the project and execution system
R and Rstudio, Windows 10

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

### real_data_analysis.R and related output CSV files
This is the R files including all the content related to real data analysis. It contains the following sections:

1. Data importing and dataset merging
2. Smooth FPCA on B-spline Non-parametric Regression
3. Clustering

From this R files execultion, 3 more CSV files will be generated to save the output of the real application:

1. mindiff_day1_nona.csv: The date list with collected temperature difference data.
2. mindiff_state1_nona.csv: The unsmoothed raw data of daily absoluted temperature difference of 50 U.S. states.
3. real_analysis_result.csv: Final clustering result from real data analysis

#### The package "fiftystater" is used for clustering result visualization on U.S. map. It may be unavailable on CRAN, so I attached the achieved fold of the package, the "fiftystater_1.0.1.tar.gz" in the Github repository.
 

### simulation study analysis.R
This is the R files including all the content related to simulation study. The main content of this file includes:
1. Functions for Data simulation and Analysis
2. Smooth FPCA on B-spline Non-parametric Regression
3. Clustering Number Selection Validation
4. Clustering Validation Study

### simulation on real analysis.R
This is the R files including all the content related to the analysis from the simulation based on real data. **Please run the code real_data_analysis.R first** for real application so as to get the output files from the real data application. The main content of this file includes:
1. Functions for Data simulation and Analysis
2. Clustering Number Selection Validation
3. Clustering Validation Study

### null case simulation.R
This is the R files including all the content related to the null case simulation study (detecting cluster number K = 1). The main content of this file includes:
1. Functions for Data simulation and Analysis
2. Clustering Number Selection Validation
