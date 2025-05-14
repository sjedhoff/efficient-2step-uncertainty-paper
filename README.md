
# Efficient Uncertainty Propagation in Bayesian Two-Step Procedures

__Authors:__  Svenja Jedhoff, Hadi Kutabi, Anne Meyer, Paul-Christian Bürkner

------------------------------------------------------------------------

## Overview

This repository contains code to replicate key experiments from our
[paper]() 'Efficient Uncertainty Propagation in Bayesian Two-Step Procedures' as well
as running the iterative algorithm for the multiple imputation setup for new datasets.

<img src="images/algorithm_multiple_imputation.png" width="100%"/>

------------------------------------------------------------------------


## Replicating

### Multiple Imputation


### Surrogate Models


------------------------------------------------------------------------

## Running the algorithm on your dataset

The iterative algorithm in the special case of regression with missing data
is presented in 
```
iterative_algo_imputation.R
```
The algorithm needs a dataset as input and a formula for specifying the regression model.
Additional specifications for the model, which will be calculated in 
[brms](https://paulbuerkner.com/brms/) can be passed through the ... operator. 
The algorithm return two lists, where the first one 'draws' contains a brms 
draws matrix corresponding to each of the imputed datasets, which are saved in the second list.

------------------------------------------------------------------------

## Citation


