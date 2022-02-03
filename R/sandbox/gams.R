## learning about GAMs
library(tidyverse)
library(mgcv)
library(itsadug)


## additive models fit a curve through data to minimize residuals, but control how wiggly the line can get

## generate data 
set.seed(10)
n = 250
x = runif(n, 0, 5)
y_model = 3 * x/(1 + 2 * x)
y_obs = rnorm(n, y_model, 0.1)

## plot data 
data_plot = qplot(x, y_obs) + geom_line(aes(y = y_model)) + theme_bw()
data_plot

## use GAM to model orginary least squares linear regresson 
linear_model = gam(y_obs ~ x)
model_summary = summary(linear_model)
print(model_summary)

## plot the regression
data_plot = data_plot + geom_line(colour = "red", aes(y = fitted(linear_model)))
print(data_plot)

## plot the residuals 
plot(residuals.gam(linear_model))

## fit a gam
gam_model = gam(y_obs ~ s(x)) ## s(x) where x is the covariate which the smooth is a function of
summary(gam_model)
plot(gam_model)

data_plot = data_plot + geom_line(colour = "blue", aes(y = fitted(gam_model)))
print(data_plot)


## test whether smoothing improves model fit
## set our smoothed model so that it is nested in our linear model
linear_model = gam(y_obs ~ x)
nested_gam_model = gam(y_obs ~ s(x) + x)

print(anova(linear_model, nested_gam_model, test = "Chisq"))
## non-linear term is significant 



## new data
n <- 250
x_test <- runif(n,-5,5)
y_test_fit <- 4*dnorm(x_test)
y_test_obs <- rnorm(n,y_test_fit, 0.2)

data_plot <- qplot(x_test, y_test_obs) + geom_line(aes(y = y_test_fit)) + theme_bw()

## use gam to model linearly
linear_model <- gam(y_test_obs ~ x_test)
summary(linear_model) ## R2 = 0.3211

## plot the regression
data_plot = data_plot + geom_line(colour = "red", aes(y = fitted(linear_model)))
print(data_plot)

## fit a gam
gam_model = gam(y_test_obs ~ s(x_test)) ## s(x) where x is the covariate which the smooth is a function of
summary(gam_model) ## R2 = 0.864
plot(gam_model)

data_plot = data_plot + geom_line(colour = "blue", aes(y = fitted(gam_model)))
print(data_plot)

## test whether smoothing improves model fit
## set our smoothed model so that it is nested in our linear model
linear_model = gam(y_test_obs ~ x_test)
nested_gam_model = gam(y_test_obs ~ s(x_test) + x_test)


print(anova(linear_model, nested_gam_model, test = "Chisq"))
## non-linear term is significant 
summary(nested_gam_model)$s.table

## generate dataset 
gam_data = gamSim(eg = 5)
head(gam_data)

## start with mode that groups by x0 (categorical, effect on intercept) and smoothed relationship with x1
basic_model = gam(y ~ x0 + s(x1), data = gam_data)
basic_summary = summary(basic_model)
print(basic_summary$p.table) ## significance of parametric terms

print(basic_summary$s.table) ## significance of smooth terms
## edf parameter (effective degrees of freedom) indicates wiggliness - higher = more wiggly
## 1 is a straight line, 8-10 is highly non-linear 

plot(basic_model)

ggplot(gam_data, aes(x = x1, y = y, colour = x0)) + geom_point()  +
  geom_line(aes(y = fitted(basic_model), group = x0))


## add a linear term x2 
two_term_model <- gam(y ~ x0 + s(x1) + x2, data = gam_data)
two_term_summary <- summary(two_term_model)
print(two_term_summary$p.table)
print(two_term_summary$s.table)
plot(two_term_model, page = 1)

## try making it non-linear 
two_smooth_model <- gam(y ~ x0 + s(x1) + s(x2), data = gam_data)
two_smooth_summary <- summary(two_smooth_model)
print(two_smooth_summary$p.table)
print(two_smooth_summary$s.table)
plot(two_smooth_model, page = 1)

# When more than one covariable is included in the model, as above, the fitted 
# response can be partitioned into the contributions of each variable as shown. 
# Here we can appreciate the varying magnitude of the effect of each variable; 
# where the y-axis represents the contribution (effect) of each covariate to the
# fitted response, centered on 0. If the confidence intervals had overlapped with 
# zero for certain values of x (or throughout the entire range), this would imply 
# a non-significant effect at those x values (or of x in entirety). When the contribution 
# for an individual covariate changes along the range of x-axis, the change in that 
# covariate is associated with a change in the response.

anova(basic_model, two_term_model, two_smooth_model, test = "Chisq")  ## best fit model = 3

## challenge 2
# Create two new models, with x3 as a linear and smoothed term. Use plots, coefficient 
# tables and the anova function to determine if x3 is an important term to include.

new_model_linear <- gam(y ~ x0 + s(x1) + s(x2) + x3, data = gam_data)
new_model_linear_summary <- summary(new_model_linear)
print(new_model_linear_summary$p.table)
print(new_model_linear_summary$s.table)
plot(new_model_linear, page = 1)

new_model_smooth <- gam(y ~ x0 + s(x1) + s(x2) + s(x3), data = gam_data)
new_model_smooth_summary <- summary(new_model_smooth)
print(new_model_smooth_summary$p.table)
print(new_model_smooth_summary$s.table)
## it's linear (edf = 1)
plot(new_model_smooth, page = 1)

anova(two_smooth_model, new_model_smooth, test = "Chisq")  ## s(x3) is not significant 


## INTERACTIONS
# if one variable is smoothed, and the other isn’t: use the by argument
# by argument lets you have a smooth term vary between different levels of a factor

# Let's ask whether the non-linear smoother s(x2) varies across different levels of x0. 
# To determine whether smooth terms differ significantly among factor levels, 
# we will use an ANOVA test on the interaction term.

categorical_interact <- gam(y ~ x0 + s(x1) + s(x2, by = x0),
                            data = gam_data)

categorical_interact_summary <- summary(categorical_interact)

print(categorical_interact_summary$s.table)

plot(categorical_interact, page = 1)


# or alternatively: plot using vis.gam function, where
# theta is the degree rotation on the x-y plane
vis.gam(categorical_interact, 
        view = c("x2", "x0"), 
        theta = 40,
        n.grid = 500, 
        border = NA)

categorical_interact <- gam(y ~ x0 + s(x1) + s(x2) + s(x2, by = x0), data = gam_data)

anova(two_smooth_model, categorical_interact, test = "Chisq")
## anova and plot confirm that we do not need interaction 


smooth_interact <- gam(y ~ x0 + s(x1, x2), data = gam_data)
smooth_interact_summary <- summary(smooth_interact)
print(smooth_interact_summary$s.table)
plot(smooth_interact, page = 1, scheme = 3)
plot(smooth_interact,page=1, scheme=1) # will give a similar plot to the vis.gam()

vis.gam(smooth_interact, 
        view = c("x1", "x2"), 
        theta = 40, 
        n.grid = 500,
        border = NA)

smooth_interact <- gam(y ~ x0 + s(x1) + s(x2) + s(x1, x2), data = gam_data)

anova(two_smooth_model, smooth_interact, test = "Chisq")







