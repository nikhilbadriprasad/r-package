# r-package

A new R package called climr which allows for:
1) Loading of global, northern hemisphere, and southern hemisphere temperature data sets
2) Fitting of models to said data 
3) Nice plots of the fitted models

A new method for the class “climr” called gp_fit. This function takes in input the class object and a string argument, which should be one of c(“Nelder-Mead”,
“BFGS”, “SANN”, “Brent”). The method gp_fit should should optimise the hyperparameters of the Gaussian process using optim with the optimisation method given by the string argument above. The method should return an object of class climr_gp_fit.

A plot method for the class climr_gp_fit. Similarly to the plot method for a climr_fit object, this method will use ggplot2 to illustrate the data as a scatterplot, and will add the smoothed Guassian process regression line.
