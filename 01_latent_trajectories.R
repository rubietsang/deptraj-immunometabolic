# Latent class trajectory modelling

library(tidyverse)
library(lcmm)
library(LCTMtools) # devtools::install_github("hlennon/LCTMtools")
library(gtools)

# Class enumeration
# Run 1-class models to get starting values
m1 <- hlme(smfq ~ ageyr + I(ageyr^2),
           subject = 'id',
           data = dep_long_inc,
           maxiter = 500,
           nproc = 10)

m1i <- hlme(smfq ~ ageyr + I(ageyr^2),
            random = ~ 1,
            subject = 'id',
            data = dep_long_inc,
            maxiter = 500,
            nproc = 10)

m1s <- hlme(smfq ~ ageyr + I(ageyr^2),
            random = ~ ageyr + I(ageyr^2),
            subject = 'id',
            data = dep_long_inc,
            maxiter = 500,
            nproc = 10)

# Run GMMs with increasing number of classes
gmmis2 <- gridsearch(rep = 50, maxiter = 10, minit = m1i, cl = 10,
                     hlme(smfq ~ ageyr + I(ageyr^2),
                          random = ~ 1,
                          mixture = ~ ageyr + I(ageyr^2),
                          ng = 2, subject = 'id', nwg = TRUE,
                          data = dep_long_inc, maxiter = 500, nproc = 10,
                          verbose = TRUE))

gmmis3 <- gridsearch(rep = 50, maxiter = 10, minit = m1i, cl = 10,
                     hlme(smfq ~ ageyr + I(ageyr^2),
                          random = ~ 1,
                          mixture = ~ ageyr + I(ageyr^2),
                          ng = 3, subject = 'id', nwg = TRUE,
                          data = dep_long_inc, maxiter = 500, nproc = 10,
                          verbose = TRUE))

gmmis4 <- gridsearch(rep = 50, maxiter = 10, minit = m1i, cl = 10,
                     hlme(smfq ~ ageyr + I(ageyr^2),
                          random = ~ 1,
                          mixture = ~ ageyr + I(ageyr^2),
                          ng = 4, subject = 'id', nwg = TRUE,
                          data = dep_long_inc, maxiter = 500, nproc = 10,
                          verbose = TRUE))

gmmis5 <- gridsearch(rep = 50, maxiter = 10, minit = m1i, cl = 10,
                     hlme(smfq ~ ageyr + I(ageyr^2),
                          random = ~ 1,
                          mixture = ~ ageyr + I(ageyr^2),
                          ng = 5, subject = 'id', nwg = TRUE,
                          data = dep_long_inc, maxiter = 500, nproc = 10,
                          verbose = TRUE))

gmmis6 <- gridsearch(rep = 50, maxiter = 10, minit = m1i, cl = 10,
                     hlme(smfq ~ ageyr + I(ageyr^2),
                          random = ~ 1,
                          mixture = ~ ageyr + I(ageyr^2),
                          ng = 6, subject = 'id', nwg = TRUE,
                          data = dep_long_inc, maxiter = 500, nproc = 10,
                          verbose = TRUE))

# Compare model fit and assess model adequacy
summarytable(m1i, gmmis2, gmmis3, gmmis4, gmmis5, gmmis6,
             which = c("conv", "G", "npm", "loglik",
                       "BIC", "SABIC", "ICL",
                       "%class", "entropy"))

summaryplot(gmmis2, gmmis3, gmmis4, gmmis5)

LCTMtoolkit(gmmis4)

# Test alternative model specifications
gmmi4 <- gridsearch(rep = 50, maxiter = 10, minit = m1i, cl = 10,
                    hlme(smfq ~ ageyr + I(ageyr^2),
                         random = ~ 1,
                         mixture = ~ ageyr + I(ageyr^2),
                         ng = 4, subject = 'id', nwg = FALSE,
                         data = dep_long_inc, maxiter = 500, nproc = 10,
                         verbose = TRUE))

gbtm4 <- gridsearch(rep = 50, maxiter = 10, minit = m1, cl = 10,
                    hlme(smfq ~ ageyr + I(ageyr^2),
                         mixture = ~ ageyr + I(ageyr^2),
                         ng = 4, subject = 'id', nwg = FALSE,
                         data = dep_long_inc, maxiter = 500, nproc = 10,
                         verbose = TRUE))

summarytable(gbtm4, gmmi4, gmmis4,
             which = c("conv", "G", "npm", "loglik",
                       "BIC", "SABIC", "ICL",
                       "%class", "entropy"))
LCTMtoolkit(gbtm4)
LCTMtoolkit(gmmi4)
LCTMtoolkit(gmmis4)