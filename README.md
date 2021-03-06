# gphmm
Generalized Pair Hidden Markov Model (GPHMM)

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Travis-CI Build Status](https://travis-ci.org/fperraudeau/gphmm.svg?branch=master)](https://travis-ci.org/fperraudeau/gphmm)

This package trains a GPHMM and computes GPHMM probabilities.

The model will is described in details in https://www.biorxiv.org/content/early/2017/12/04/228619.

## Installation

```{r}
install.packages('gphmm')
```

Note that `gphmm` package need compilation of C++ code.

Installing from CRAN is recommended. However, if you want to install the version from GitHub, you can do so with the following lines of code

```{r}
library(devtools)
install_github("fperraudeau/gphmm")
```
