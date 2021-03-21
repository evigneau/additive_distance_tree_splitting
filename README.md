
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AddDistTreeSplit

<!-- badges: start -->

<!-- badges: end -->

AddDistTreeSplit includes R functions for recursive splitting of an
addirive distance tree

## Installation

You can install the released version of AddDistTreeSplit from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("AddDistTreeSplit")
```

## Example

``` r
library(AddDistTreeSplit)
data(DistAnimal)
tree=ape::nj(DistAnimal)
ape::plot.phylo(tree,type="unrooted")
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
recursive_partitioning(d=DistAnimal)
#> [[1]]
#> [[1]][[1]]
#> [1] "Cat"   "Dog"   "Lion"  "Mouse"
#> 
#> [[1]][[2]]
#> [1] "Bear"   "Cow"    "Deer"   "Goat"   " Horse" "Pig"    "Rabbit" "Sheep" 
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] "Cat"   "Dog"   "Lion"  "Mouse"
#> 
#> [[2]][[2]]
#> [1] "Deer"   "Rabbit"
#> 
#> [[2]][[3]]
#> [1] "Bear"   "Cow"    "Goat"   " Horse" "Pig"    "Sheep"
```

## How to cite

Koenig, L., Cariou, V., Symoneaux, R., Coulon-Leroy, C., Vigneau, E.
(2021). Food Quality and Preference, 89, 104137(2021).
