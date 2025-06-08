
<!-- README.md is generated from README.Rmd. Please edit that file -->

# argminCS

Welcome! Here we have source code to perform argmin hypothesis test.

## Overview

The goal of argminCS is to produce confidence set of argmin from iid
samples with a valid type 1 control, while exhibiting desirable
statistical power. In particular, the method ‘softmin.LOO’ is the main
innovative component in our paper *Winners with Confidence: Argmin
Inference over a High-Dimensional Discrete Candidate Set*. Several other
methods are also implemented within the package to ease method
comparison and simulations.

## Citation

If you use the **argminCS** package for your research or any experiment,
please cite our paper “Winners with Confidence: Argmin Inference over a
High-Dimensional Discrete Candidate Set”.

## Installation

You can install the development version of argminCS from this
[GitHub](https://github.com/) webpage with:

``` r
# install.packages("devtools")
devtools::install_github("xu3cl4/argminCS")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(argminCS)
dimension <- 4
sample.size <- 200
p <- 20
mu <- (1:p)/p
cov <- diag(length(mu))
set.seed(108)
data <- MASS::mvrnorm(sample.size, mu, cov)

## to test if 'dimension' is likely to be argmin with (default) softmin.LOO
difference.matrix <- matrix(rep(data[, dimension], p-1), 
                            ncol = p-1, 
                            byrow = FALSE) - data[, -dimension]
argmin.HT(difference.matrix, method='SML')
#> $test.stat.scale
#> [1] -0.7101986
#> 
#> $critical.value
#> [1] 1.644854
#> 
#> $std
#> [1] 0.7251399
#> 
#> $ans
#> [1] "Accept"
#> 
#> $lambda
#> [1] 5.656854
#> 
#> $lambda.capped
#> [1] FALSE
#> 
#> $residual.slepian
#> [1] 0.04882178
#> 
#> $variance.bound
#> [1] 0.03146688
#> 
#> $test.stat.centered
#> NULL

## rather than perform a hypothesis testing for a specific dimension, 
## one can directly generate a discrete confidence set by 
CS.argmin(data, method='SML')
#> [1] 1 2 5 6
```

## Detailed Tutorial

Here is a detailed
[tutorial](https://xu3cl4.github.io/argminCS/demo_CSargmin.html).

Regarding the details of methods and their associated tuning parameters,
we encourage users to install the package and check function
documentation.
<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->

## Key References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-cck.many.moments" class="csl-entry">

Chernozhukov, V., D. Chetverikov, and K. Kato. 2013. “Testing Many
Moment Inequalities.” IDEAS Working Paper Series from RePEc.

</div>

<div id="ref-dey.2024" class="csl-entry">

Dey, N., M. Ryan, and J. P. Williams. 2024. “Anytime-Valid Generalized
Universal Inference on Risk Minimizers.” *arXiv.org*.
<https://doi.org/10.48550/arxiv.2402.00202>.

</div>

<div id="ref-futschik.1995" class="csl-entry">

Futschik, Andreas, and Georg Pflug. 1995. “Confidence Sets for Discrete
Stochastic Optimization.” *Annals of Operations Research* 56 (1):
95–108. <https://doi.org/10.1007/BF02031702>.

</div>

<div id="ref-gupta.1965" class="csl-entry">

Gupta, Shanti S. 1965. “On Some Multiple Decision (Selection and
Ranking) Rules.” *Technometrics* 7 (2): 225–45.
<https://doi.org/10.1080/00401706.1965.10490251>.

</div>

<div id="ref-lei.cvc" class="csl-entry">

Lei, Jing. 2020. “Cross-Validation with Confidence.” *Journal of the
American Statistical Association* 115 (532): 1978–97.
<https://doi.org/10.1080/01621459.2019.1672556>.

</div>

<div id="ref-mogstad2024inference" class="csl-entry">

Mogstad, Magne, Joseph P Romano, Azeem M Shaikh, and Daniel Wilhelm.
2024. “Inference for Ranks with Applications to Mobility Across
Neighbourhoods and Academic Achievement Across Countries.” *Review of
Economic Studies* 91 (1): 476–518.

</div>

</div>
