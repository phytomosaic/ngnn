#' @title Out-of-sample prediction for Nonparametric Constrained
#'      Ordination
#'
#' @description
#' Calculates new scores for out-of-sample observations in a
#'      multivariate NCO gradient space.
#'
#' @param obj object of class 'nco' from call to \code{\link{nco}}
#'
#' @param method distance measure for all ordinations
#'
#' @param neighb number of adjacent distances considered in prediction
#'
#' @param maxits number of iterations
#'
#' @param type either 'points' or 'text' for plotting
#'
#' @param cexn expansion factor for points and text
#'
#' @param ocol color value or vector for out-of-sample points
#'
#' @param ... additional arguments passed to function
#'
#' @return
#' List of class 'ncopredict' including new out-of-sample predicted
#' scores, and flags indicating which (if any) were beyond the range
#' of original scores.  Specifically, the first item \code{nmsp} is a
#' list of 5 items describing predicted NCO scores:
#' \itemize{
#'   \item scr_both = combined new and old NCO scores
#'   \item scr_o = new NCO scores only
#'   \item stress = overall stress of the solution, with new scores
#'   \item iters = number of iterations performed
#'   \item cause = ICAUSE from FORTRAN, reason for termination of
#'       iterations: 1 = max iterations used, 2 = stress fell below
#'       STRMIN, 3 = stress ratio exceeded SRATMX, 4 = scale factor of
#'       gradient fell below SFGRMN
#'   \item R2_enviro = squared correlation of D and Dz
#'   \item R2_partial = squared correlation of Dz and each predictor
#'   \item Axis_tau = rank correlation of each axis and predictor
#'   }
#'
#' @details
#' Combines existing algorithms in multivariate workflow:
#'
#'         NPMR + NMS = NCO[in-sample] \cr
#'         NCO + NCOpredict  = NCO[out-of-sample]
#'
#' The main function is a wrapper for \code{add.points()}, written by
#' \href{http://ecology.msu.montana.edu/labdsv/R/labs/lab9/lab9.html}{
#' Dave Roberts}, who deserves all credit. The intellectual pedigree
#' of the NMS prediction concept includes Dave Roberts, Peter Minchin,
#' Jari Oksanen, and Bruce McCune, although they each differ in
#' algorithmic approach. McCune et al. (1997a, 1997b) gave its first
#' use in the literature, while McMurray et al. (2015) gave the first
#' results from the R implementation.  PC-ORD (McCune and Mefford
#' 2016) seems to be the only other software that offers NMS
#' prediction, and does so with a different algorithm that also flags
#' poorly fit new points.
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
#' # NPMR
#' res_npmr <- npmr(spe, idi, ido, nm, nmulti=5)
#' summary(res_npmr)
#'
#' # NCO (NMS)
#' res_nco  <- nco(res_npmr, method='bray', thresh=0.90)
#' summary(res_nco)
#'
#' # NCOpredict (NMSpredict)
#' res_nmsp <- nco_predict(res_nco, method='bray', neighb=5,
#'                         maxits=999)
#' summary(res_nmsp)
#'
#' # plot the NCO gradient space
#' par(mfrow=c(1,2))
#' plot(res_nco)       # original in-sample points
#' plot(res_nmsp)      # add out-of-sample points
#' par(mfrow=c(1,1))
#'
#' @family ncopredict functions
#'
#' @references
#' McCune, B., J.P. Dey, J.E. Peck, K. Heiman, and S. Will-Wolf. 1997a.
#'   Regional gradients in lichen communities of the southeast United
#'   States. Bryologist 100:145-158.
#'
#' McCune, B., J.P. Dey, J.E. Peck, D. Cassell, K. Heiman, S.
#'   Will-Wolf, and P.N. Neitlich. 1997b. Repeatability of community
#'   data: species richness versus gradient scores in large-scale lichen
#'   studies. Bryologist 100:40-46.
#'
#' McCune, B., and M. J. Mefford. 2016. PC-ORD. Multivariate Analysis
#'     of Ecological Data. Version 7. MjM Software Design, Gleneden
#'     Beach, OR.
#'
#' McMurray, J.A., D.W. Roberts, and L.H. Geiser. 2015. Epiphytic
#'   lichen indication of nitrogen deposition and climate in the
#'   northern rocky mountains, USA. Ecological Indicators 49:154-161.
#'
#' Roberts, D.W. 2017. LabDSV: Non-metric Multidimensional Scaling.
#'   URL http://ecology.msu.montana.edu/labdsv/R/labs/lab9/lab9.html
#'   [original code, as \code{add.points()}]
#'
#' @export
#' @rdname nco_predict
### wrapper for NMSpredict, get scores for out-of-sample sample units
`nco_predict` <- function(obj, method, neighb, maxits, ...){
     stopifnot(class(obj)=='nco')
     if(!identical(names(obj$iYhat), names(obj$oYhat))){
          stop('Species must agree between old and new Yhat')
     }
     cat('Getting distances of combined in/out-of-sample data...\n\n')
     D0 <- vegan::vegdist(rbind(obj$iYhat, obj$oYhat), method)
     cat('Predicting new NCO scores...\n\n')
     nmsp <- NMSpredict(scr=obj$scr_i, dis=D0, neighb, maxits)
     names(nmsp)[names(nmsp)=='newpoints'] <- 'scr_both'
     names(nmsp)[names(nmsp)=='newpoints'] <- 'scr_o'
     # flag if beyond range of existing scores
     flagax1 <- nmsp$scr_o[
          which(nmsp$scr_o[,1] > max(obj$scr_i[,1])|
                     nmsp$scr_o[,1] < min(obj$scr_i[,1])),]
     flagax2 <- nmsp$scr_o[
          which(nmsp$scr_o[,2] > max(obj$scr_i[,2]) |
                     nmsp$scr_o[,2] < min(obj$scr_i[,2])),]
     out <- c(obj,
              list(nmsp=nmsp,
                   flagax1=flagax1,
                   flagax2=flagax2))
     class(out) <- 'ncopredict'
     out
}

### NMSpredict core algorithm lightly adapted from add.points()
`NMSpredict` <- function(scr, dis, neighb, maxits){
     # convert original scores to local matrix, and get sizes
     points  <- list(points=scr)
     class(points) <- 'nmds'
     points <- points$points
     oldn <- nrow(points)
     ndim <- ncol(points)
     totn <- attr(dis,'Size')
     newn <- totn - oldn
     # TODO: allow 1 point instead of many - until then, force > 1
     if (newn < 1) stop('come on, try more than just 1 point!')
     # test correspondence
     if (!identical(dimnames(points)[[1]],attr(dis,'Labels')[1:oldn]))
          stop('ordination and dissimilarity matrix do not match')
     # decompose dissimilarity object to isolate new values
     diss <- as.matrix(dis)[1:oldn,(oldn+1):totn]
     # set up initial conditions
     ndis <- oldn * newn
     tmp <- matrix(rep(0,newn*ndim),ncol=ndim)
     for (i in 1:newn) {
          pnt <- seq(1,oldn)[order(diss[,i])][1:neighb]
          weight <- 1-diss[pnt,i]
          for (j in 1:ncol(points)) {
               tmp[i,j] <- weighted.mean(points[pnt,j],w=weight)
          }
     }
     xinit <- rbind(points,tmp)
     dimnames(xinit)[[1]] <- attr(dis,'Labels')
     # set up indices
     iidx <- rep((1:oldn),newn)
     jidx <- NULL
     for (i in (oldn+1):totn) jidx <- c(jidx,rep(i,oldn))
     # set up ordination
     nfix <- oldn
     ngrp <- istart <- 1
     isform <- 2
     ities <- 1
     iregn <- 1
     iscal <- 0
     sratmx <- 0.99999
     strmin <- 1e-07
     sfgrmn <- 1e-05
     dist <- rep(0,ndis)
     dhat <- rep(0,ndis)
     x <- matrix(0,nrow=totn,ncol=ndim)
     stress <- 1
     strs <- ngrp
     iters <- 1
     icause <- 1
     maxits <- maxits
     iwork <- rep(0,ndis)
     grad <- matrix(0,nrow=totn,ncol=ndim)
     grlast <- matrix(0,nrow=totn,ncol=ndim)
     out <- .Fortran('monoMDS',
                     nobj=as.integer(totn),
                     nfix=as.integer(nfix),
                     ndim=as.integer(ndim),
                     ndis=as.integer(ndis),
                     ngrp=as.integer(ngrp),
                     diss=as.double(diss),
                     iidx=as.integer(iidx),
                     jidx=as.integer(jidx),
                     xinit=as.double(xinit),
                     istart=as.integer(istart),
                     isform=as.integer(isform),
                     ities=as.integer(ities),
                     iregn=as.integer(iregn),
                     iscal=as.integer(iscal),
                     maxits=as.integer(maxits),
                     sratmx=as.double(sratmx),
                     strmin=as.double(strmin),
                     sfgrmn=as.double(sfgrmn),
                     dist=as.double(dist),
                     dhat=as.double(dhat),
                     points=as.double(x),
                     stress=as.double(stress),
                     strs=as.double(strs),
                     iters=as.integer(iters),
                     cause=as.integer(icause))
     res <- list(points=matrix(out$points,ncol=ndim),
                 newpoints=matrix(out$points,ncol=ndim)[(oldn+1):totn,],
                 stress=out$stress,
                 iters=out$iters,
                 cause=out$cause)
     dimnames(res$points)[[1]] <- attr(dis,'Labels')
     dimnames(res$newpoints)[[1]] <- attr(dis,'Labels')[(oldn+1):totn]
     class(res) <- 'nmds'
     res
}
#' @export
#' @rdname nco
# summary method
`summary.ncopredict` <- function(obj, ...){
     stopifnot(class(obj)=='ncopredict')
     out <- list(nmsp       = obj$nmsp,
                 flagax1    = obj$flagax1,
                 flagax2    = obj$flagax2)
     out
}
#' @export
#' @rdname nco_predict
# plot both old and new scores
`plot.ncopredict` <- function(obj, type='points', ocol=2, cexn=NULL,
                              ...){
     stopifnot(class(obj)=='ncopredict')
     type <- match.arg(type, c('points', 'text', 'none'))
     ylim <- c(min(obj$scr_i[,2],obj$nmsp$scr_o[,2])*1.1,
               max(obj$scr_i[,2],obj$nmsp$scr_o[,2])*1.1)
     xlim <- c(min(obj$scr_i[,1],obj$nmsp$scr_o[,1])*1.1,
               max(obj$scr_i[,1],obj$nmsp$scr_o[,1])*1.1)
     if(!is.null(cexn)){
          cexn <- normalize(obj$id_i[cexn])
     } else { cexn <- 0.7 }
     # par(oma=rep(0,4), mar=c(4,4,1,1))
     ordiplot(obj$scr_i, type=type, display='sites', xlim=xlim,
              ylim=ylim, cex=cexn, las=1)
     if(type=='points'){
          points(obj$nmsp$scr_o, col=ocol, pch=16, cex=0.7)
     }
     if(type=='text'){
          text(obj$nmsp$scr_o, labels=row.names(obj$nmsp$scr_o),
               col=ocol, cex=0.7)
     }
}
###   end   ###