#' @title NCO, Nonparametric constrained ordination
#'
#' @description
#' Constrained ordination based on nonparametric regression and NMS.
#'
#' @param obj object of class 'npmr' from call to \code{\link{npmr}}
#'
#' @param method distance measure for all ordinations
#'
#' @param thresh numeric threshold for stepacross dissimilarities
#'
#' @param type either 'points' or 'text' for plotting
#'
#' @param cexn expansion factor for points and text
#'
#' @param ... additional arguments passed to function
#'
#' @return
#' List of class 'nco' with elements:
#' \itemize{
#'   \item scr_i = environmentally constrained scores from NCO
#'   \item NCO_model = the NCO model itself
#'   \item R2_internal = squared correlation of Dhat and Dz
#'   \item R2_enviro = squared correlation of D and Dz
#'   \item R2_partial = squared correlation of Dz and each predictor
#'   \item Axis_tau = rank correlation of each axis and predictor
#'   }
#'
#' @details
#' Combines existing algorithms in multivariate workflow:
#'
#'         NPMR + NMS = NCO
#'
#' NCO (McCune and Root 2012; McCune and Root 2017) is nonmetric
#' multidimensional scaling ordination (NMS; Kruskal 1964) of sample
#' units based on environmentally constrained fitted values from a set
#' of nonparametric multiplicative regressions (NPMR; McCune 2006) for
#' each species. This implementation uses \code{\link[np]{npreg}} from
#' package `np` (NPMR stage) and \code{\link[vegan]{metaMDS}} from
#' package 'vegan' (NMS stage).
#'
#' Variance explained \code{R2_enviro} measures the strength of the
#' relationship between the original community data and the ordination
#' as constrained to the *measured* predictors.  Like all constrained
#' ordinations, NCO does not account for variation in species
#' composition related to *unmeasured* predictors.
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
#' # NPMR basic usage
#' res_npmr <- npmr(spe, idi, ido, nm, nmulti=5)
#' summary(res_npmr)
#' plot(res_npmr, pick=1:9, nm=nm)
#'
#' # NCO basic usage
#' res_nco <- nco(res_npmr, method='bray', thresh=0.90)
#' summary(res_nco)
#'
#' # plot the NCO gradient space
#' plot(res_nco)
#' plot(res_nco, type='text')
#'
#' @family nco functions
#'
#' @references
#' Kruskal, J. B. 1964. Multidimensional scaling by optimizing
#'   goodness of fit to a nonmetric hypothesis. Psychometrika 29:
#'   1-27.
#'
#' McCune, B. 2006. Non-parametric habitat models with automatic
#'   interactions. Journal of Vegetation Science 17(6):819-830.
#'
#' McCune, B., and H. T. Root. 2012. Nonparametric constrained
#'   ordination. 97th ESA Annual Meeting. Ecological Society of
#'   America, Portland, OR.
#'
#' McCune, B., and H. T. Root. 2017. Nonparametric constrained
#'   ordination to describe community and species relationships to
#'   environment. Unpublished ms.
#'
#' @export
#' @rdname nco
`nco` <- function(obj, method, thresh, ...){
     stopifnot(class(obj)=='npmr')
     # unconstrained AND constrained distances, optionally stepacross
     cat('Calculating dissimilarities...\n\n')
     if (!missing(thresh)){
          D    <- vegan::stepacross(
               vegan::vegdist(obj$spe, method=method),
               path='shortest', toolong=thresh)
          Dhat <- vegan::stepacross(
               vegan::vegdist(obj$iYhat, method=method),
               path='shortest', toolong=thresh)
     } else {
          D    <- vegan::vegdist(obj$spe, method=method)
          Dhat <- vegan::vegdist(obj$iYhat, method=method)
     }
     # ordinate constrained values
     cat('\nPerforming NCO, may take a moment...\n\n')
     nax <- max(obj$nm_len, 2)  # constrain n axes = n predictors
     ord <- vegan::metaMDS(Dhat, k=nax, try=99, trymax=99, maxit=500,
                           autotransform=FALSE, weakties=TRUE, ...)
     row.names(ord$points) <- row.names(obj$spe)
     names(ord$points)     <- c(paste0('NCO', 1:nax))
     scr_i <- ord$points                # NCO calibration scores
     Dz  <- stats::dist(scr_i, 'euc')   # NCO ordination distances
     R2_internal <- stats::cor(Dhat, Dz)^2 # R2 btwn constrnd Dhat,Dz
     R2_enviro <- stats::cor(D, Dz)^2   # R2 btwn orig spp dist and Dz
     R2_partial  <- matrix(rep(NA,nax), nrow=nax, ncol=1, byrow=T,
                           dimnames=list(obj$nm[1:obj$nm_len],
                                         'R2_partial'))
     for (i in 1:obj$nm_len){
          R2_partial[i,] <-
               stats::cor(Dz, stats::dist(obj$id_i[obj$nm[[i]]]))^2
     }
     Axis_tau <- stats::cor(obj$id_i, scr_i, method='kendall')
     out_nco <- c(obj,
                  list(scr_i        = scr_i,
                       NCO_model    = ord,
                       R2_internal  = R2_internal,
                       R2_enviro    = R2_enviro,
                       R2_partial   = R2_partial,
                       Axis_tau     = Axis_tau))
     class(out_nco) <- 'nco'
     out_nco
}
#' @export
#' @rdname nco
# summary method
`summary.nco` <- function(obj, ...){
     stopifnot(class(obj)=='nco')
     out <- list(scr_i        = obj$scr_i,
                 NCO_model    = obj$NCO_model,
                 R2_internal  = obj$R2_internal,
                 R2_enviro    = obj$R2_enviro,
                 R2_partial   = obj$R2_partial,
                 Axis_tau     = obj$Axis_tau)
     out
}
#' @export
#' @rdname nco
# plot method
`plot.nco` <- function(obj, type='points', cexn=NULL, ...){
     stopifnot(class(obj)=='nco')
     type <- match.arg(type, c('points', 'text', 'none'))
     ylim <- c(min(obj$scr_i[,2])*1.1, max(obj$scr_i[,2])*1.1)
     xlim <- c(min(obj$scr_i[,1])*1.1, max(obj$scr_i[,1])*1.1)
     if(!is.null(cexn)){
          cexn <- normalize(obj$id_i[cexn])
     } else { cexn <- 0.7 }
     vegan::ordiplot(obj$scr_i, type=type, display='sites',
                     xlim=xlim, ylim=ylim, cex=cexn, las=1, ...)
}
# automatic scale for cex (this function is also in package `ecole`)
`normalize` <- function(x, ...){
     (x - min(x, ...))/(max(x, ...) - min(x, ...))
}
###   end   ###