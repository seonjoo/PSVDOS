# PSVDOS

## CRAN
To install from github:
```{r}
library(devtools)
install_github("seonjoo/PSVDOS")
library(PSVDOS)
```

To install the most recent package from CRAN type:
```{r}
install.packages("PSVDOS")
library(PSVDOS)
```

## Usage
The PSVDOS function conducts Poisson singular value decomposition with offset estimation. 

```{r}
data(demo)
heatmap(Y)
psvdfit = PSVDOS(Y,K=4,verbose=1,err = 0.0001,niter = 300)
names(psvdfit)
```
