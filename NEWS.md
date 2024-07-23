# ppgam 1.0.2

* Removed set.seed(1) from `ppgam()`, which clearly shouldn't have been left there! (Thanks to Patrik Bohlinger for spotting this.)

* Added argument `weight.non.numeric` to function `ppgam()`, which allows weights to be calculated for non-numeric variables according to their frequency of occurrence.

* Added functionality for S3 methods, `coef()`, `logLik()`, `simulate()` and `fitted()`.

# ppgam 1.0.1

* `ppgam()` has had a couple of lines changed to account for changes in `evgam:::.joinSmooth()`

# ppgam 1.0.0

* Initial release
