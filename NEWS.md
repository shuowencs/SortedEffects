NEWS
================
Shuowen Chen
5/5/2021

# SortedEffects 1.4.0.
Changes from SortedEffects 1.3.0. For nonparametric bootstrap, we use the multinomial weight resampling to produce more stable bootstrap estimates. 

Changes from SortedEffects 1.2.0. The spe function in the 1.2.0 version produces far wider nonparametric bootstrap confidence bands for some specifications. We fix that issue. 

Changes from SortedEffects 1.1.0. The 1.1.0 version has a check problem due to a recent upgrade of package 'tibble' to version 3.0. The package doesn't really need tibble package, so we change the data type. 
