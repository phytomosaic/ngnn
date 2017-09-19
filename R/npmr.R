#' @title NPMR, Nonparametric multiplicative regression
#'
#' @description
#' Nonparametric kernel regression with automatic interactions among
#'     multiple predictors, adapted here for in-sample estimation AND
#'     out-of-sample prediction.  Also known as 'generalized product
#'     kernel' regression in econometrics.
#'
#' @param spe species dataframe, rows = sample units and columns =
#'     species
#'
#' @param idi in-sample predictor dataframe, rows must match 'spe'
#'
#' @param ido out-of-sample predictor dataframe, where rows = new
#'     sample units
#'
#' @param nm  string vector specifying predictors to include (max 2)
#'
#' @param nmulti number of random starts in nonparametric regression
#'
#' @param pa logical, convert to presence/absence?
#'
#' @param pr logical, use 'beals' for probs of joint occurrence?
#'
#' @param obj object of class 'npmr' from call to \code{npmr}
#'
#' @param pick variable to query
#'
#' @param zlim vector of length 2, giving vertical limits for plots
#'
#' @param ... additional arguments passed to function
#'
#' @return
#' List of class 'npmr' with elements:
#' \itemize{
#'   \item spe = original species matrix
#'   \item id_i = in-sample predictors used in NPMR
#'   \item nm = which predictors were used
#'   \item nm_len = their length
#'   \item iYhat = in-sample fitted values from NPMR
#'   \item oYhat = out-of-sample fitted values from NPMR
#'   \item np_stat = fit, tolerances and results of signif tests
#'   \item np_mods = list of NPMR regression models for every species
#'   \item np_bw = list of NPMR  bandwidths for every species
#'   }
#'
#' @details
#' NPMR is nonparametric multiplicative regression (McCune 2006). This
#' implementation uses a variant \code{\link[np]{npreg}} from package
#' `np`, where it is known as generalized product kernel regression
#' (Li and Racine 2007). The current function is used in a predictive
#' capacity for estimating environmentally constrained fitted values
#' for NGNN, and therefore requires a predictor matrix for
#' out-of-sample observations.
#'
#' Sensitivity analysis is available, following Eqn 9 in McCune
#' (2006:825), aka Sensitivity Formula 1 in Hyperniche software.
#' Sensitivity (Q) is the mean absolute difference in the response
#' resulting from nudging the predictors +/- 5%, expressed as a
#' proportion of the range of the response variable. Q = 1 means that
#' nudging a predictor results in a change in response of equal
#' magnitude, and Q = 0 means that nudging a predictor has no
#' detectable effect on the response.
#'
#' @examples
#' # set up
#' set.seed(978)
#' require(vegan)
#' data(varespec, varechem)
#' spe <- varespec ; id  <- varechem
#' i   <- sample(1:nrow(spe), size=floor(0.75*nrow(spe))) # sample
#' spe <- spe[i, ]          # in-sample species
#' idi <- id[i, ]           # in-sample predictors
#' ido <- id[-i, ]          # out-of-sample predictors
#' nm  <- c('Al', 'K')      # select 1 or 2 gradients of interest
#'
#'
#' # NPMR basic usage
#' res_npmr <- npmr(spe, idi, ido, nm, nmulti=5)
#' summary(res_npmr)
#'
#' # sensitivity analysis
#' Q <- npmr_sens(obj=res_npmr, pick = 'Vaccmyrt', nm)
#' Q
#'
#' # plot NPMR species response curves
#' plot(res_npmr, pick=4, nm)   # alternatively: pick='Vaccmyrt'
#' plot(res_npmr, pick=1:9, nm)
#'
#' @family npmr functions
#'
#' @references
#' Li, Q. and J. S. Racine. 2007. Nonparametric Econometrics: Theory
#'   and Practice, Princeton University Press.
#'
#' McCune, B. 2006. Non-parametric habitat models with automatic
#'   interactions. Journal of Vegetation Science 17(6):819-830.
#'
#' @export
#' @rdname npmr
### NPMR, nonparametric multiplicative regression
# TODO: parallelize
`npmr` <- function(spe, idi, ido, nm, nmulti, pa, pr, ...){
     # coerce, then subset predictors of interest
     spe <- data.frame(spe)
     idi <- data.frame(idi)
     ido <- data.frame(ido)
     nm_len <- length(nm)
     id_i <- idi[ ,nm]
     id_o <- ido[ ,nm]
     stopifnot(identical(row.names(spe), row.names(id_i)))
     stopifnot(identical(names(id_i), names(id_o)))
     ## TODO: admit factors in generalized product kernel
     #     but for now, coerce factors to numeric
     `unfac` <- function(f){ as.numeric(levels(f))[f] }
     isfac   <- sapply(id_i, is.factor)
     id_i[isfac] <- lapply(id_i[isfac], FUN=unfac)
     id_o[isfac] <- lapply(id_o[isfac], FUN=unfac)
     # convert presence/absence
     if(missing(pa)) pa <- FALSE
     if(pa) spe[spe>0] <- 1
     # Beals smoothing for _pr_obabilities of joint occurrence
     if(missing(pr)) pr <- FALSE
     if(pr) spe <- vegan::beals(x=spe, type=1)
     # fit NPMR models for in-sample sample units (TIME WARN ! ! ! !)
     #  bw0 heuristic for bandwidth selection equal to 5% the range:
     # bw0 <- sapply(id_i, FUN=function(x){diff(range(x))*0.05})
     `fn1` <- function(x){
          bw  <- np::npregbw(
               stats::as.formula(paste('abund ~',
                                       paste(nm,collapse='+'))),
               data=x, nmulti=nmulti)
     }
     `fn2` <- function(x){
          mod <- np::npreg(bw=x)
     }
     cat('Calculating kernel bandwidths, may take a moment...\n')
     tmp   <- reshape2::melt(cbind(spe, id_i), id.vars=nm,
                             variable.name='spp', value.name='abund')
     bw    <- plyr::dlply(tmp, plyr::.(spp), fn1, .progress='text')
     cat('Performing np regressions...\n')
     mod   <- plyr::llply(bw, .fun=fn2, .progress = 'text')
     iYhat <- plyr::ldply(mod, .fun=stats::fitted)
     row.names(iYhat) <- iYhat[,1]
     iYhat <- iYhat[,-1]
     iYhat <- data.frame(t(iYhat))
     row.names(iYhat) <- row.names(id_i)
     # extract summaries
     `fn3` <- function(x) {
          c(round(x$R2,3), round(x$MSE, 2), round(x$MAPE,2),
            format(x$bw, scientific=FALSE))
     }
     cat('Getting in-sample regression summaries...\n')
     np_stat <- plyr::ldply(mod, fn3, .progress = 'text')
     names(np_stat) <- c('SPE','R2','MSE','MAPE',
                         paste0('tol_',1:nm_len))
     # NPMR predictions for out-of-sample sample units
     cat('Predictions for out-of-sample SUs, may take a moment...\n')
     `fn4` <- function(x){
          preds <- stats::predict(x, newdata = id_o)
     }
     oYhat <- plyr::ldply(mod, fn4, .progress = 'text')
     sppnames <- oYhat[,1]
     oYhat <- oYhat[,-1]
     oYhat <- data.frame(t(oYhat))
     names(oYhat) <- as.character(sppnames)
     row.names(oYhat) <- row.names(ido)
     cat('Remember to inspect species response curves...\n')
     out_npmr <- list(spe     = spe,
                      id_i    = id_i,
                      nm      = nm,
                      nm_len  = nm_len,
                      iYhat   = iYhat,
                      oYhat   = oYhat,
                      np_stat = np_stat,
                      np_mods = mod,
                      np_bw   = bw)
     class(out_npmr) <- 'npmr'
     out_npmr
}
#' @export
#' @rdname npmr
# sensitivity analysis # TO DO: need to iterate large number of times
`npmr_sens` <- function(obj, pick, nm, ...){
     stopifnot(class(obj)=='npmr')
     stopifnot(obj$nm_len==2)
     if(length(pick)!=1) stop('Choose only 1 species')
     if(is.numeric(pick)){ pick <- names(spe)[pick] }
     Q   <- rep(NA, length=obj$nm_len) # initialize vector
     spe <- obj$spe
     id  <- obj$id_i
     mod <- obj$np_mods
     # nudge 1st predictor 5% its range
     nudge <- 0.05*diff(range(id[,obj$nm[[1]]]))
     eval <- data.frame(
          c( id[,obj$nm[[1]]] + nudge, id[,obj$nm[[1]]] - nudge),
          c( id[,obj$nm[[2]]] , id[,obj$nm[[2]]] ),
          rep(pick, nrow(id)*2 ))
     names(eval) <- c(paste0(obj$nm[[1]]), paste0(obj$nm[[2]]), 'spp')
     cat('Fitting nudged values for var 1, may take a moment...\n\n')
     nudgpr <- stats::predict(obj = mod[[pick]],
                              data = obj$spe[,pick],
                              newdata = eval)
     numer <- sum(abs( nudgpr - c(spe[, pick], spe[, pick])))
     denom <- 2*nrow(id)*diff(range( spe[, pick] ))*0.05
     Q[1]  <- numer/denom
     # nudge 2nd predictor 5% its range
     nudge <- 0.05*diff(range(id[,obj$nm[[2]]] ))
     eval <- data.frame(
          c( id[,obj$nm[[1]]] , id[,obj$nm[[1]]]),
          c( id[,obj$nm[[2]]] + nudge , id[,obj$nm[[2]]] - nudge ),
          rep(pick, nrow(id)*2 ))
     names(eval) <- c(paste0(obj$nm[[1]]), paste0(obj$nm[[2]]), 'spp')
     cat('Fitting nudged values for var 2, may take a moment...\n\n')
     nudgpr <- stats::predict(obj = mod[[pick]],
                              data = obj$spe[,pick],
                              newdata = eval)
     numer <- sum(abs( nudgpr - c(spe[, pick], spe[, pick])))
     denom <- 2*nrow(id)*diff(range( spe[, pick] ))*0.05
     Q[2] <- numer/denom
     Q <- round(Q, 3)
     names(Q) <- nm
     Q
}
#' @export
#' @rdname npmr
# summary method
`summary.npmr` <- function(obj, ...){
     stopifnot(class(obj)=='npmr')
     out <- list(nm        = obj$nm,
                 np_stat   = obj$np_stat)
     out
}
#' @export
#' @rdname npmr
# inspect species response curves from NPMR
`plot.npmr` <- function(obj, pick=NULL, nm, zlim, ...){
     stopifnot(class(obj)=='npmr')
     spe <- obj$spe
     obj <- obj$np_mods
     ev <- 30
     wasmissing <- missing(zlim)
     if(is.null(pick)) pick <- c(1:ncol(spe))
     if(is.character(pick)) pick <- which(names(spe) %in% pick)
     `fn5` <- function (nn=length(pick)) {
          if (nn <= 3)  c(1, nn)
          else if (nn <= 6)  c(2, (nn + 1)%/%2)
          else if (nn <= 12) c(3, (nn + 2)%/%3)
          else c(ceiling( nn / (nr <- ceiling(sqrt(nn)))), nr)
     }
     x1 <- seq(min(obj[[1]]$eval[[1]]),max(obj[[1]]$eval[[1]]),len=ev)
     x2 <- seq(min(obj[[1]]$eval[[2]]),max(obj[[1]]$eval[[2]]),len=ev)
     dd <- expand.grid(x1, x2)
     names(dd) <- nm
     graphics::par(mfrow=c(fn5()[1], fn5()[2]),
                   oma=c(0,0,0,0), mar=c(0,0,.9,0))
     cat('Plotting species response curves, just a moment...\n')
     for (i in pick){
          f <- matrix(stats::predict(obj[[i]], newdata=dd), ev, ev)
          if(wasmissing) {
               if(identical(min(f), max(f))) mn <- 0 else mn <- min(f)
               zlim <- c(mn*0.9, max(f)*1.1)
          }
          graphics::persp(x1, x2, f, col='lightblue',
                          main=names(obj)[[i]],
                          xlab=nm[[1]], ylab=nm[[2]], zlab='',
                          theta=125, phi=35, d=1.5, ltheta = -30,
                          lphi = 55, shade=0.9, ticktype='d',
                          expand=0.7, zlim=zlim)
     }
}
###   end   ###