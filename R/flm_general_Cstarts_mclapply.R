pflm <- function(formula, nsd = 1, mc.cores = 1, npf, tbounds, nruns, startd = NULL, dx, knotsx, 
                 pars, db, knotsb = NULL, criterion = c("A","D"), lambda = 0,
                 dlbound = -1, dubound = 1, tol = 0.0001, progress = FALSE){
  
  if (mc.cores < 1) stop("The number of cores must be at least one")
  if (nsd < 1) stop("The number of starting designs must be at least one")
  if (tbounds[1] != 0) stop("The time line should start from zero")
  tu <- max(tbounds)
  if (tu <= 0) stop("Argument t should be the time interval from 0 to a positive time value")
  if (npf <= 0) stop("The number of profile factors should be positive")
  if (any(dx < 0)) stop("Degree entries in dx cannot be negative")
  if (any(db < 0)) stop("Degree entries in db cannot be negative")
  if (npf != length(dx)) stop("The length of dx should be equal to the number of profile factors")
  if (npf != length(knotsx)) stop("The length of knotsx should be equal to the number of profile factors")
  if (!(all(pars %in% c("power","bspline")))) stop("Parameters should be power or bspline")
  if (!(criterion %in% c("A","D"))) stop("Criterion should be A or D")
  if (lambda < 0) stop("Lambda cannot be negative")
  if (tol <= 0) stop("Tolerance should be a positive value")
  if (dubound <= dlbound) stop("The upper bound of the design should be greater than the lower bound")
  
  ptm <- proc.time()[3]
  
  for (ii in 1:npf) {
    if (!is.null(knotsx[[ii]])) {
      if (min(knotsx[[ii]]) < 0 | max(knotsx[[ii]]) > tu) stop("Knots for profile factors are allowed only in the time line")
    }
    knotsx[[ii]] <- as.numeric(knotsx[[ii]])
  }
  
  no.terms <- length(attr(terms(formula), "term.labels"))
  
  for (jj in 1:no.terms) {
    if (!is.null(knotsb[[jj]])){
      if (min(knotsb[[jj]]) < 0 | max(knotsb[[jj]]) > tu) stop("Knots for parameters are allowed only in the time line")
    }
    knotsb[[jj]] <- as.numeric(knotsb[[jj]])
  }
  
  if (no.terms < npf) stop("The number of terms in the formula must be at least equal to npf")
  if (length(knotsb) != no.terms) stop("The length of knotsb should be equal to the number of terms in the formula")
  if (length(db) != no.terms) stop("The length of db should be equal to the number of terms in the formula")
  if (length(pars) != no.terms) stop("The length of pars should be equal to the number of terms in the formula")
  
  nb <- rep(0, no.terms)
  for (g in 1:no.terms) {
    nb[g] <- length(knotsb[[g]]) + db[[g]] + 1
  }
  
  inter <- attr(terms(formula), "intercept")
  dimnb <- inter + sum(nb)
  if (nruns < dimnb & identical(criterion, "A")) stop("For invertible information matrix, the number of runs must be greater or equal than the number of basis functions for the functional parameters")
  
  if (identical(criterion, "A")) {
    flinear <- linearcpp.aopt
  } else {
    flinear <- linearcpp.dopt
  }
  
  nx <- rep(0, npf)
  for (r in 1:npf) {
    nx[r] <- length(knotsx[[r]]) + dx[[r]] + 1
    if (nx[r] < nb[r] & identical(criterion, "A")) stop("For invertible information matrix, the number of basis functions for a function of a profile factor must be greater or equal to the number of basis functions for a functional parameter")
  }
  
  if (lambda == 0) {
    v <- matrix(0, nrow=dimnb, ncol=dimnb)
  } else {
    vtemp <- list()
    for (i in 1:no.terms) {
      if (pars[i]=="power") {
        vtemp[[i]] <- Vpowercpp(db=db[i], tu=tu)
      } else {
        vtemp[[i]] <- Vbsplinecpp(db=db[i], tu=tu, knots=knotsb[[i]])
      }
    }
    v <- as.matrix((rbind(0, cbind(0, bdiag(vtemp))) + 1 - 1) * lambda)
    if (inter == 1) {
      v <- v
    } else {
      v <- v[,-1]
    }
  }
  
  if (!is.null(startd)) {
    if (!is.list(startd)) stop("Argument startd should be NULL or a list for the starting design of length npf")
    if (length(startd) != nsd) stop("Argument startd should be NULL or a list of length nsd, each corresponding to a starting design of length npf")
    d <- startd
    for(w in 1:nsd) {
      for (i in 1:npf) {
        if (nrow(d[[w]][[i]]) != nruns) stop("The number of rows in objects in startd should be equal to nruns")
        if (ncol(d[[w]][[i]]) != nx[i]) stop("The number of columns in objects in startd should be equal to the number of basis of each profile factor")
      }
    }
  } else {
    dd <- list()
    d <- list()
    for (c in 1:nsd) {
      for (i in 1:npf) {
        dd[[i]] <- matrix(runif(nruns * nx[i], dlbound, dubound), nrow=nruns, ncol=nx[i])
        names(dd)[i] <- paste0("x", i, sep="")
      }
      d[[c]] <- dd
    }
    startd <- d
  }
  
  if (no.terms == npf) {
    result <- mclapply(X=d, FUN=flm.ce.me, mc.cores=mc.cores, t=tbounds, n=nruns, npf=npf, 
                       db=db, knotsb=knotsb, tol=tol, dx=dx, knotsx=knotsx, pars=pars,
                       dlbound=dlbound, dubound=dubound, progress=progress, nb=nb, nx=nx,
                       flinear=flinear, v=v, inter=inter, ptm=ptm)
  } else {
    #no.ints <- no.terms - nrow(attr(terms(formula), "factors"))
    nrfactors <- nrow(attr(terms(formula), "factors"))
    no.ints <- no.terms - nrfactors
    no.pols <- nrfactors - npf
    if (nrfactors == npf) {
      result <- mclapply(X=d, FUN=flm.ce.ints, mc.cores=mc.cores, formula=formula, t=tbounds, 
                         n=nruns, npf=npf, dx=dx, knotsx=knotsx, pars=pars, flinear=flinear,
                         no.terms=no.terms, no.ints=no.ints, nb=nb, nx=nx, v=v,
                         db=db, knotsb=knotsb, tol=tol, inter=inter,
                         dlbound=dlbound, dubound=dubound, progress=progress, ptm=ptm)
    } else {
      result <- mclapply(X=d, FUN=flm.ce.intspols, mc.cores=mc.cores, formula=formula, t=tbounds, 
                         n=nruns, npf=npf, dx=dx, knotsx=knotsx, pars=pars, flinear=flinear,
                         no.terms=no.terms, no.pols=no.pols, no.ints=no.ints, nb=nb, nx=nx,
                         db=db, knotsb=knotsb, v=v, tol=tol, inter=inter,
                         dlbound=dlbound, dubound=dubound, progress=progress, ptm=ptm)
    }
  }
  
  allobjvals <- rep(0, nsd)
  alldesigns <- list()
  for (c in 1:nsd) {
    allobjvals[c] <- result[[c]]$objval
    alldesigns[[c]] <- result[[c]]$design
  }
  bestrep <- which.min(allobjvals)
  result1 <- result[[bestrep]]
  
  result2 <- list(startd=startd[[bestrep]], tbounds=tbounds, npf=npf, criterion=criterion, 
                nruns=nruns, formula=formula, dx=dx, knotsx=knotsx, 
                lambda=lambda, dbounds=c(dlbound, dubound), bestrep=bestrep,
                allobjvals=allobjvals, alldesigns=alldesigns, allstartd=startd)
  
  output <- c(result1, result2)
  class(output) <- "flm"
  output
}


# resall.mc <- pflm(formula=~x1+x2+x1:x2+P(x2,2), t=c(0,1), nsd=3, n=20, npf=2, 
#             dx=c(2,3), knotsx=list(c(0.33,0.66), c(0.25,0.50,0.75)), 
#               pars=c("bspline","bspline","bspline","bspline"), startd=NULL,
#                 db=c(0,1,2,2), knotsb=list(c(0.5), c(0.5),c(0.25,0.75), c(0.25,0.75)), 
#                   criterion="D", tol=0.00001, dlbound=-1, dubound=1, lambda=0, progress=TRUE)
# 
# 
# res <- pflm(formula=~x1, t=c(0,1), nsd=2, n=4, npf=1, 
#          dx=c(0), knotsx=list(c(0.5)), 
#          pars=c("power"), startd=NULL,
#          db=c(1), knotsb=list(c()), 
#          criterion="A", tol=0.00001, dlbound=-1, dubound=1, lambda=0, progress=TRUE)

# compd.func <- function(allres){
#   C <- length(allres)
#   objvals <- rep(0,C)
#   for (c in 1:C) {
#     objvals[c] <- allres[[c]]$objective.value
#   }
#   objmin <- which.min(objvals)
#   allres[[objmin]]
# }
# 
# compd.func(resall)
# compd.func(resall.mc)


