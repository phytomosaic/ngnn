# ngnn

Nonlinear Gradient Nearest Neighbors, a method for predicting multivariate species compositions from environmental variables based on individualistic but possibly coordinated species responses. Nonlinear responses are possible, which improves upon existing gradient nearest neighbors methods. Package also includes methods for other ecological tasks, including in-sample NMS ordination, out-of-sample NMS prediction, Nonparametric Constrained Ordination (NCO), Nonparametric Multiplicative Regression (NPMR), and simplified Gradient Nearest Neighbors (GNN).

## Motivation

Given a set of sample units where species abundances and corresponding predictor values are both known, how does one infer which species should appear in 'new' sample units where only the predictors are known?  NGNN (nonlinear gradient nearest neighbors) approaches the problem of species imputation in the following way:

Regress species individualistically on predictors -> Feed fitted values to NMS ordination -> Find nearest neighbors in ordination space, and assign those species.  This retains realistic communities of co-occurring species, since they've already been observed in at least one other sample unit.  

## Installation

Install the package from github as follows:
```r
install.packages('devtools')
devtools::install_github('phytomosaic/foggy')
devtools::install_github('phytomosaic/ecole')
```

## Examples

Create datasets (here split one artificially)
```r
# set up
set.seed(978)
require(vegan)
data(varespec, varechem)
spe <- varespec ; id  <- varechem
i   <- sample(1:nrow(spe), size=floor(0.75*nrow(spe))) # sample
spe <- spe[i, ]          # in-sample species
idi <- id[i, ]           # in-sample predictors
ido <- id[-i, ]          # out-of-sample predictors
nm  <- c('Al', 'K')      # select 1 or 2 gradients of interest
```

Basic usage
```r
# perform NGNN
res <- ngnn(spe, idi, ido, nm, nmulti=5, method='bray',
            thresh=0.90, neighb=5, maxits=999, k=1)
summary(res)
str(res, 1)

# plot the species response curves
ngnn_plot_spp(res, pick=1:9, nm=nm)

# plot the NCO gradient space
ngnn_plot_nco(res)

# predicted (imputed) species composition for out-of-sample sites
ngnn_get_spp(res)

# how close were predicted species composition to 'true' values?
spe_append <- rbind(spe, res$spp_imputed)   # append to existing
heatmap(t(as.matrix(spe_append)), Rowv=NA, Colv=NA)

# check composition of 'hold-out' data
heatmap(t(as.matrix(varespec[-i,])), Rowv=NA, Colv=NA)
# ... vs new species from NGNN
heatmap(t(as.matrix(res$spp_imputed)), Rowv=NA, Colv=NA)

# Prediction error: Root Mean Square Error
`rmse` <- function(y, ypred, ...){
     sqrt(mean((y-ypred)^2, ...))
}
rmse(varespec[-i,], res$spp_imputed)
```

Can do this manually, avoiding the 'ngnn' wrapper function:
```r
# step 1: NPMR
res_npmr <- npmr(spe, idi, ido, nm, nmulti=5)

# step 2: NCO (NMS)
res_nco  <- nco(res_npmr, method='bray', thresh=0.90)

# step 3: NCOpredict (NMSpredict)
res_nmsp <- nco_predict(res_nco, method='bray', neighb=5, maxits=999)

# step 4: GNN
res_gnn  <- gnn(obj=res_nmsp, k=1)
summary(res_gnn)
```

## References

Kruskal, J. B. 1964. Multidimensional scaling by optimizing goodness of fit to a nonmetric hypothesis. Psychometrika 29: 1-27.
McCune, B. 2006. Non-parametric habitat models with automatic interactions. Journal of Vegetation Science 17(6):819-830. 
McCune, B., and H. T. Root. 2012. Nonparametric constrained ordination. 97th ESA Annual Meeting. Ecological Society of America, Portland, OR.
McCune, B., and H. T. Root. 2017. Nonparametric constrained ordination to describe community and species relationships to environment. Unpublished ms.
Ohmann, J.L., and M.J. Gregory. 2002. Predictive mapping of forest composition and structure with direct gradient analysis and nearest-neighbor imputation in coastal Oregon, U.S.A. Canadian Journal of Forest Research 32:725-741.
