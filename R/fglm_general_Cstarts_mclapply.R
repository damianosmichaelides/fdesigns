pfglm <- function(formula, nsd = 1, mc.cores = 1, npf, tbounds, nruns, startd = NULL, 
                  dx, knotsx, pars, db, knotsb = NULL, lambda = 0, criterion = c("A","D"),
                  family, method = c("quadrature", "MC"), level = NULL, B = NULL, prior, 
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
  if (!(method %in% c("quadrature","MC"))) stop("Method should be quadrature or MC")
  if (identical(method, "MC") & !is.function(prior)) stop("For MC method, argument prior must be a function.")
  if (identical(method, "quadrature") & !is.list(prior)) stop("For quadrature method, argument prior must be a list with specific names; see the help file")
  if (!(criterion %in% c("A","D"))) stop("Criterion should be A or D")
  if (lambda < 0) stop("Lambda cannot be negative")
  if (tol <= 0) stop("Tolerance should be a positive value")
  if (dubound <= dlbound) stop("The upper bound of the design should be greater than the lower bound")

  for (ii in 1:npf) {
    if (!is.null(knotsx[[ii]])) {
      if (min(knotsx[[ii]]) < 0 | max(knotsx[[ii]]) > tu) stop("Knots for profile factors are allowed only in the time line")
    }
    knotsx[[ii]] <- as.numeric(knotsx[[ii]])
  }

  no.terms <- length(attr(terms(formula), "term.labels"))

  for (jj in 1:no.terms) {
    if (!is.null(knotsb[[jj]])){if (min(knotsb[[jj]]) < 0 | max(knotsb[[jj]]) > tu) stop("Knots for parameters are allowed only in the time line")}
    knotsb[[jj]] <- as.numeric(knotsb[[jj]])
  }
  
  if(no.terms < npf) stop("The number of terms in the formula must be at least equal to npf")
  if(length(knotsb) != no.terms) stop("The length of knotsb should be equal to the number of terms in the formula")
  if(length(db) != no.terms) stop("The length of db should be equal to the number of terms in the formula")
  if(length(pars) != no.terms) stop("The length of pars should be equal to the number of terms in the formula")

  nb <- rep(0, no.terms)
  for (g in 1:no.terms) {
    nb[g] <- length(knotsb[[g]]) + db[[g]] + 1
  }

  inter <- attr(terms(formula), "intercept")
  notinter <- no.terms + inter
  dimnb <- inter + sum(nb)
  if (nruns < dimnb & identical(criterion, "A")) stop("For invertible information matrix, the number of runs must be greater or equal than the number of basis functions for the parameters plus the intercept")

  nx <- rep(0, npf)
  for (r in 1:npf) {
    nx[r] <- length(knotsx[[r]]) + dx[[r]] + 1
    if (nx[r] < nb[r] & identical(criterion, "A")) stop("For invertible information matrix, the number of basis functions for a function of a profile factor must be greater or equal to the number of basis functions for a functional parameter")
  }
  
  if (!is.function(family)) {
    if (is.list(family)) {
      famlist <- family
    }
    else {
      famlist.get <- get(family, mode = "function", envir = parent.frame())
      famlist <- famlist.get()
    }
  } else {
    famlist <- family()
  }
  
  if(famlist$family != "binomial" & famlist$family != "poisson") stop("Family or link not implemented for negative squared error loss (NSEL) utility yet")
  if(identical(famlist$family, "binomial") & !identical(famlist$link, "logit")) stop("The link specified for the binomial family is not available yet; see the help file for the available families and links")
  if(identical(famlist$family, "poisson") & !identical(famlist$link, "log")) stop("The link specified for the poisson family is not available yet; see the help file for the available families and links")

  if (identical(criterion, "A")) {
    if (identical(famlist$family, "binomial")) {
      fgeneralised <- logisticcpp.aopt
    } else {
      fgeneralised <- poissoncpp.aopt
    }
  } else {
    if (identical(famlist$family, "binomial")) {
      fgeneralised <- logisticcpp.dopt
    } else {
      fgeneralised <- poissoncpp.dopt
    } 
  }

  ptm <- proc.time()[3]
  
  if (is.null(B) & identical(method, "MC")) {
    B <- 10000
  } else {
    B <- B
  }
  
  if (is.null(level) & identical(method, "quadrature")) {
    level <- 5
  } else {
    level <- level
  }
 
  if (identical(method, "quadrature")) {
    if (identical(names(prior)[1:2], c("mu", "sigma2"))) { 
      if (is.vector(prior$mu)) {
        if (identical(length(prior$mu), as.integer(1))) {
          priormu <- rep(prior$mu, dimnb)
        } else {
          if (identical(length(prior$mu), as.integer(notinter))) {
            priormu <- list()
            if (inter == 1) {nb1 <- c(1, nb)} else {nb1 <- nb}
            for (kk in 1:notinter) { 
              priormu[[kk]] <- rep(prior$mu[kk], nb1[kk])
            }
            priormu <- unlist(priormu)
          } else {
            stop("prior$mu should be scalar or a vector of the length of the number of terms in the formula")
          }
        }
      } else {
        stop("prior$mu should be scalar or a vector")
      }
      if (is.vector(prior$sigma2)) {
        lps <- length(prior$sigma2)
        if (identical(lps, as.integer(1))) {
          psig2 <- prior$sigma2
          priorvarcovar <- diag(psig2, nrow=dimnb)
        } else {
          if (identical(lps, as.integer(notinter))) {
            psig2 <- prior$sigma2
            priorsigma2 <- list()
            if (inter == 1) {nb1 <- c(1, nb)} else {nb1 <- nb}
            for (ll in 1:notinter){
              priorsigma2[[ll]] <- psig2[ll] * diag(nb1[ll])
            }
            priorvarcovar <- as.matrix(bdiag(priorsigma2)) + 1 - 1
          } else {
            stop("prior$sigma2 should be a scalar for a scalar common variance,
                  a vector of the length of the number of terms in the formula,
                  or a square matrix with number of rows and columns equal to 
                  the number of terms in the formula")
          }
        }
      } else {
          if (is.matrix(prior$sigma2)) {
            nrps <- nrow(prior$sigma2)
            ncps <- ncol(prior$sigma2)
            if (nrps != notinter | ncps != notinter) stop("prior$sigma2 should be a scalar 
                  for a scalar common variance, a vector of the length of the 
                  number of terms in the formula, or a square matrix with number 
                  of rows and columns equal to the number of terms in the formula")
            psig2 <- prior$sigma2
            priorsigma2 <- list()
            if (inter == 1) {nb1 <- c(1, nb)} else {nb1 <- nb}
            for (ll in 1:notinter) {
              priorsigma2[[ll]] <- diag(psig2)[ll] * diag(nb1[ll])
            }
            priorvarcovar <- as.matrix(bdiag(priorsigma2)) + 1 - 1
          } else {
            stop("prior$sigma2 should be a scalar for a scalar common variance,
                  a vector of the length of the number of terms in the formula,
                  or a square matrix with number of rows and columns equal to 
                  the number of terms in the formula")
          }
      }
      grid.nig <- createNIGrid(dim=dimnb, type="GHe", level=level) 
      rescale(grid.nig, m=priormu, C=priorvarcovar, dec.type=2)
      abc <- getNodes(grid.nig)
      qw <- as.vector(getWeights(grid.nig))
      wei <- qw * dmvnorm(abc, priormu, priorvarcovar)
    } else { 
      if (identical(names(prior)[1], c("unifbound"))) {
        punifb <- prior$unifbound
        if (is.vector(punifb)) {
          if (!identical(length(punifb), as.integer(2))) stop("If the bounds for a Uniform
                  prior are specified as a vector, the length should be 2")
          if (punifb[2] < punifb[1]) {
            punifb <- sort(punifb)
            warning("The upper bound of the Uniform bound was smaller than the lower bound.
                    The bounds were reversed.")
          }
          punifb.mat <- matrix(rep(punifb, dimnb), nrow=2, ncol=dimnb)
        } else {
          if (is.matrix(punifb)) {
            if (!identical(dim(punifb), as.integer(c(2, notinter)))) stop("If 
                 the bounds for a Uniform prior are specified as a matrix, the dimension 
                 should be c(2, number of terms in the formula)")
            punifb.mat <- NULL
            if (inter == 1) {nb1 <- c(1, nb)} else {nb1 <- nb}
            for (ll in 1:notinter) {
              punifb.mat <- cbind(punifb.mat, matrix(rep(c(punifb[,ll]), nb1[ll]), nrow=2, ncol=nb1[ll]))
            }
          }
          grid.nig <- createNIGrid(dim=dimnb, type="GLe", level=level)
          rescale(grid.nig, domain=t(punifb.mat), dec.type=2)
          abc <- getNodes(grid.nig)
          wei <- as.vector(getWeights(grid.nig))
        }
      } else {
        stop("Argument prior is required to specify a uniform or a normal distribution
             as a list of certain names; see the help file")
      }
    }
  } 

  if (identical(method, "MC")){ 
    prMCnb <- prior(B=B, Q=dimnb)
    if (!identical(dim(prMCnb), as.integer(c(B, dimnb)))) stop("The prior function must return a matrix of dimensions c(B, Q)")
    abc <- prMCnb
    wei <- rep(1/B, B)
  }

  if (lambda == 0){
    v <- matrix(0, nrow=dimnb, ncol=dimnb)
  } else {
      vtemp <- list()
      for (i in 1:no.terms) {
        if (pars[i]=="power") {
          vtemp[[i]] <- Vpowercpp(db=db[i], tu=tu)
        } else {
          vtemp[[i]] <- Vbsplinespp(db=db[i], tu=tu, knots=knotsb[[i]])
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
    result <- mclapply(X=d, FUN=fglm.me, mc.cores=mc.cores, tbounds=tbounds, formula=formula,  
                n=nruns, npf=npf, dx=dx, knotsx=knotsx, pars=pars, db=db, knotsb=knotsb,
                v=v, nb=nb, nx=nx, abc=abc, wei=wei, fgeneralised=fgeneralised, inter=inter, 
                ptm=ptm, tol=tol, dlbound=dlbound, dubound=dubound, progress=progress)
  } else {
      nrfactors <- nrow(attr(terms(formula), "factors"))
      no.ints <- no.terms - nrfactors
      no.pols <- nrfactors - npf
      if (nrfactors == npf) {
        result <- mclapply(X=d, FUN=fglm.ints, mc.cores=mc.cores, tbounds=tbounds, formula=formula,  
                           n=nruns, npf=npf, dx=dx, knotsx=knotsx, pars=pars, db=db, knotsb=knotsb,
                           v=v, nb=nb, nx=nx, abc=abc, wei=wei, fgeneralised=fgeneralised, inter=inter, 
                           ptm=ptm, tol=tol, dlbound=dlbound, dubound=dubound, progress=progress,
                           no.terms=no.terms, no.ints=no.ints, criterion=criterion)
      } else {
        result <- mclapply(X=d, FUN=fglm.intspols, mc.cores=mc.cores, tbounds=tbounds, formula=formula,  
                           n=nruns, npf=npf, dx=dx, knotsx=knotsx, pars=pars, db=db, knotsb=knotsb,
                           v=v, nb=nb, nx=nx, abc=abc, wei=wei, fgeneralised=fgeneralised, inter=inter, 
                           ptm=ptm, tol=tol, dlbound=dlbound, dubound=dubound, progress=progress,
                           no.terms=no.terms, no.ints=no.ints, no.pols=no.pols, criterion=criterion)
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
              nruns=nruns, formula=formula, family=c(famlist$family, famlist$link), 
              method=method, B=B, prior=prior, dx=dx, knotsx=knotsx, lambda=lambda,
              dbounds=c(dlbound, dubound), bestrep=bestrep, allobjvals=allobjvals, 
              alldesigns=alldesigns, allstartd=startd)
  
  output <- c(result1, result2)
  class(output) <- "fglm"
  output
}

# set.seed(100)
# pr.mc <- function(B,Q){
#   matrix(rnorm(B*Q, mean=0, sd=sqrt(2)), nrow=B, ncol=Q)
# }
# 
# res <- pfglm(formula = ~ 1 + x1, nsd = 10, mc.cores = 1, npf = 1, t=c(0,1), n=12, startd=NULL, dx=c(0),
#              knotsx=list(c(0.25,0.50,0.75)), pars=c("power"),
#              db=c(1), knotsb=list(c()), lambda=0, family=binomial,
#              method=c("quadrature"), B=10000, prior=list(mu=c(0,0), sigma2=2*diag(2)),
#              criterion="A", tol=0.0001, dlbound=-1, dubound=1, progress=TRUE)
# 
# resMC <- pfglm(formula=~1+x1+x2+P(x1,2)+x1:x2, nsd=2, mc.cores=1, npf=2, t=c(0,1), n=12, startd=NULL, dx=c(1,0),
#              knotsx=list(c(0.25,0.50,0.75), c()), pars=c("power","power","power", "bspline"),
#              db=c(1,0,0,1), knotsb=list(c(),c(),c(),c()), lambda=0, family=binomial,
#              method=c("MC"), B=10000, prior=pr.mc,
#              criterion="A", tol=0.0001, dlbound=-1, dubound=1, progress=TRUE)

# tidy_dir("R")
# 
# formula=~1+x1
# C=1
# npf=1
# t=c(0,1)
# n=12
# startd=NULL
# dx=c(1)
# knotsx=list(c(0.5))
# pars=c("power")
# db=c(0)
# knotsb=list(c())
# lambda=0
# family=binomial
# method=c("MC")
# B=10000
# prior=pr.mc
# criterion="A"
# tol=0.0001
# dlbound=-1
# dubound=1
# progress=TRUE

# compd.func <- function(allres){
#   C <- length(allres)
#   objvals <- rep(0,C)
#   for (c in 1:C) {
#     objvals[c] <- allres[[c]]$objective.value
#   }
#   objmin <- which.min(objvals)
#   allres[[objmin]]
# }


# Linear models:
# formula - formula to create the model equation and need matching list names for startd
# npf - number of profile factors
# t - time vector sequence from 0 to T, OR vector of 2 values boundaries of time
# n - number of runs
# startd - starting design but if NULL then a random design is generated
# C - number of random starts
# dx - degree of profile factors (length as value of npf)
# knotsx - list of knot vectors of profile factors (length as value of npf)
# pars - basis for parameters "power" or "bspline" length equal to npf
# db - degree of parameters, for power series db=1,2 linear and then no interior knots
# knotsb - list of knots with length equal to number of terms, NULL for power pars
# criterion - choose A- or D- optimality objective function
# tol - tolerance value change to stop optimisation algorithm
# dlbounds - design lower bounds, so function's lower bounds (length 1)
# dubounds - design upper bounds, so function's upper bounds (length 1)
# progress - shown iterations, default is false

