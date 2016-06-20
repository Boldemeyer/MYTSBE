# MYTSBE
Multi Year Time Stratified Bayesian Estimator

`MYTSBE` is an R package for estimating juvenile salmonid abundances using capture-mark-recapture collected at rotary screw traps. `MYTSBE` contains several functions for formatting juvenile salmonid data into the appropriate capture-mark-recapture format needed for the `MYTSBE` model (`Chinook_Format()` and `Steelhead_Format()`). The package also contains three functions that run the same hierarchical Bayesian model but produce different life-stage and cohort summaries dependent on the species (`MYTSBE_Cohort()`, `MYTSBE_Cohort_Fry()`, and `MYTSBE_Calendar()`). 


## Getting Started

To install `MYTSBE` you can use Hadley Wickham's `devtools` package. To install and load the `devtools` package use:
```
install.packages("devtools")
library(devtools)
```
NOTE: To use `devtools`, you may also have to download and install Rtools (although you shouldn't). The latest version on Rtools can be found at
https://cran.r-project.org/bin/windows/Rtools/

Once `devtools` is successfully installed, use the following to install SCOBI:
```
devtools::install_github("Boldemeyer/MYTSBE")
```
If you are interested in making contributions to `MYTSBE`, consider getting a GitHub account, fork this repository, clone to a local directory, modify, and send me a pull request. I can then review any changes and merge.

For further information email me or, see:

Oldemeyer, B.N., Copeland, T.S., and B.P. Kennedy. (in review. 2016). A multi-year hierarchical Bayesian mark-recapture model using recurring salmonid behavior to account for sparse or missing data. Canadian Journal of Fisheries and Aquatic Sciences. 

* Will update upon publication
