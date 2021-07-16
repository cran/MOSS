## MOSS: Multi-omic integration via sparse singular value decomposition.

[![Travis build status](https://travis-ci.com/agugonrey/MOSS.svg?branch=master)](https://travis-ci.com/agugonrey/MOSS)

Agustin Gonzalez-Reymundez, Alexander Grueneberg, and Ana I. Vazquez.

### Installing and loading MOSS from CRAN.

```
install.packages("MOSS")
library("MOSS")
```
### Installing and loading MOSS from GitHub.

```
if (require("remotes") == FALSE) install.packages("remotes")
install_github("agugonrey/MOSS")
library("MOSS")
```

### Documentation.

  For a description of the package's main function. 

```
help(moss)
```

  For more documentation, see the package's [vignette](https://github.com/agugonrey/MOSS/blob/master/inst/MOSS_working_example.pdf). An example of using MOSS on a multi-omic "big" pan-cancer data can be found [here](https://github.com/agugonrey/MOSS/blob/master/inst/MOSS_pancancer_example.pdf). The data can be found [here](https://data.mendeley.com/datasets/r8p67nfjc8/1).
