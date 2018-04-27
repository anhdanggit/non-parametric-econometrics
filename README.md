# Non-parametric Methods

I provide the R code for 02 take-home projects in non-parametric course (Master 2, EEE, TSE). This repository includes the following code files, corresponding to `[Mai-Anh Dang] Report-2.pdf` and `[Mai-Anh Dang] Report-1.pdf`

Report                | Code files
----------------------| ----------------------------
Report-1              | `mean_regression.R`
Report-1                         | `simulation_kernel.R`
Report-1                         | `density_estimate.R`
----------------------| ----------------------------                      
Report-2              | `dWADE.R`
Report-2                       | `nonparam_bootstraps.R`
Report-2                       | `qtile_reg_BG90.R`     

### Theoretical frameworks
The theoretical and formal equations behind each methods and code files could be found in the associated reports. The methods covered include:

* Kernel Estimator
* Mean Regression Functions: Local Constant, Local Linear
* Density-weighted Average Derivative Estimator (dWADE)
* Bootstraps to construct Confident Interval for non-parametric estimates
* Bhattacharya and Gangopadhyay (1990) estimator
* Quantile Regressions

### Applications
These mentioned methods are conducted and assessed through both simulations and pratical applications, using the below data:

* GDP 2005 and 2016 `data/GDP.xlsx`
* Annual Household Income and Food Expenditure in Belgium `data/Engel.dta`
* House Price and Other Charactersitics `data/anglin.gencay.1996.csv`

### Simulation: Non-parametric Density Estimate
The nonparametric kernel estimators do not fit the true density perfectly, but performing quite well, even for the small sample of n = 100. When n increase, order of error decrease. For the large sample n = 1000, the estimated density is closer to the true curve.
![alt tag](https://github.com/maianhdang/non_parametric/blob/master/results/Grahs_MC/MC-Y.png)

### Density Estimate and CI by Bootstraps
We use the Pivotal Bootstraps approach to construct the CIs for density estimates on `GDP.xlsx`. Based on these interval, we test the null hypotheses visually.  
![alt tag](https://github.com/maianhdang/non_parametric/blob/master/results/Null_Btrap.png)

### Mean Regression: Local Linear and Local Constant
To estimate the expected Food Expenditure at a given level of Income, on the data set `Engel.dta`
![alt tag](https://github.com/maianhdang/non_parametric/blob/master/results/Graphs_mean_reg/exercise2.png)

### BG90 Estimator and Quantile Regression
To estimate the expected Income at a given level of Food Expenditure, on the data set `Engel.dta`
![alt tag](https://github.com/maianhdang/non_parametric/blob/master/results/Ex3_BG90_npqreg.png)

### dWADE Estimator 
This method is applied n the data set `anglin.gencay.1996.csv` for a hedonic analysis, describing the relationship between housing price and observed characteristics. 
