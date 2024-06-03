
<!-- README.md is generated from README.Rmd. Please edit that file -->

# argminCS

Welcome! Here we have source code to perform argmin hypothesis test.
<!-- badges: start --> <!-- badges: end -->

## Overview

The goal of argminCS is to produce confidence set of argmin from iid
samples with a valid type 1 control, while exhibiting desirable
statistical power. In particular, the method ‘softmin.LOO’ is the main
innovative component in the paper *Winners with Confidence: Argmin
Inference over a High-Dimensional Discrete Candidate Set* by Tianyu
Zhang, Hao Lee and Jing Lei. Several other methods are also implemented
within the package to ease method comparison and simulations.

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
mu <- (1:20)/20
cov <- diag(length(mu))
set.seed(108)
data <- MASS::mvrnorm(sample.size, mu, cov)
sample.mean <- colMeans(data)

## to test if 'dimension' is likely to be argmin with (default) softmin.LOO
argmin.HT(data, dimension, method='SML')
#> $test.stat.scale
#> [1] 0.2384774
#> 
#> $critical.value
#> [1] 1.644854
#> 
#> $std
#> [1] 1.034184
#> 
#> $ans
#> [1] "Accept"

## rather than perform a hypothesis testing for a specific dimension, 
## one can directly generate a discrete confidence set by 
CS.argmin(data, method='SML')
#> [1] 1 2 3 4
```

Regarding the details of methods and their associated tuning parameters,
we encourage users to install the package and check function
documentation. Note that the method in the NC state paper is not
supported for now.
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

<div id="ref-lei.cvc" class="csl-entry">

Lei, Jing. 2020. “Cross-Validation with Confidence.” *Journal of the
American Statistical Association* 115 (532): 1978–97.
<https://doi.org/10.1080/01621459.2019.1672556>.

</div>

</div>
