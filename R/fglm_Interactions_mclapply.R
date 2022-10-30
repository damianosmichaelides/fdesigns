fglm.ints <- function(d, tbounds, formula, n, npf, dx, knotsx, pars, db, knotsb, v, nb,
                      nx, abc, wei, fgeneralised, inter, ptm, tol, dlbound, dubound, progress,
                      no.terms, no.ints, criterion){
  
  tu <- max(tbounds)
  eps <- tol 
  diff <- eps + 1
  
  jmat <- list()
  djmat <- list()
  for (i in 1:npf) {
    if (pars[i] == "power") {
      jmat[[i]] <- jcbfunccpp(db=db[i], degree=dx[i], knots=knotsx[[i]], t=tbounds)
    } else {
      jmat[[i]] <- bs2jcbcpp(t=tbounds, dbeta=db[i], dx=dx[i], knots_beta=knotsb[[i]], knots_x=knotsx[[i]])
    }
  }
  
  modmat <- model.matrix(object=formula, data=d)
  modframe.ints <- attr(terms(formula), "factors")
  allknotsxb <- list
  
  for (q in (npf+1):(npf+no.ints)) {
    d[[q]] <- modmat[,which(attr(modmat, "assign")==q)]
    colnames(d[[q]]) <- NULL
    names(d)[q] <- attr(terms(formula), "term.labels")[q]
    dx.ints <- modframe.ints[,q]
    allknotsxb <- knotsx[dx.ints==1]
    nx.ints <- nx[dx.ints==1]
    if (prod(nx.ints) < nb[q] & identical(criterion, "A")) stop("For A-optimality, the product of basis functions for the profile factors involved in an interaction must be greater or equal to the number of basis functions for the parameter")
    lkxi <- length(allknotsxb)
    lnxvec <- list()
    for (m in 1:lkxi) {
      lnxvec[[m]] <- 1:nx.ints[m]
    }  
    aposx <- as.matrix(expand.grid(lnxvec))
    if (pars[q] == "power") {
      jmat[[q]] <- bsalljcbpowercpp(t=tbounds, j=0, nbs=lkxi, aposx=aposx, alldeg=dx[dx.ints==1], allknotsx=allknotsxb, nx=nx.ints, nb=nb[q])
    } else {
      allknotsxb[[lkxi+1]] <- knotsb[[q]]
      jmat[[q]] <- bsalljcbcpp(t=tbounds, j=0, nbs=lkxi+1, aposx=aposx, alldeg=c(dx[dx.ints==1],db[q]), allknotsxb=allknotsxb, nx= nx.ints,  nb=nb[q])
    }
  } 
  
  djmat <- mlist(d, jmat)
  if (inter == 1) {
    z <- cbind(1, djmat)
  } else {
    z <- djmat
  }
  
  curr.eval <- fgeneralised(z=z, v=v, a=abc, w=wei) 
  it <- 1
  if (progress == TRUE) {
    cat("Starting design", c("( Current Value =", curr.eval,")"), "\n")
  }
  
  while (diff > eps) {
    
    for (i in 1:n) {
      for (p in 1:npf) {
        for (j in 1:nx[p]) {
          fnobj <- function(k) {
            f <- d
            f[[p]][i,j] <- k
            modmat <- model.matrix(object=formula, data=f)
            for (w in (npf+1):(npf+no.ints)) {
              f[[w]] <- modmat[,which(attr(modmat, "assign")==w)]
              colnames(f[[w]]) <- NULL
              names(f)[w] <- attr(terms(formula), "term.labels")[w]
            }
            if (inter == 1) {
              zopt <- cbind(1, mlist(f,jmat))
            } else {
              zopt <- mlist(f,jmat)
            }
            result <- fgeneralised(z=zopt, v=v, a=abc, w=wei)
            result
          }
          d[[p]][i,j] <- optim(par=d[[p]][i,j], fn=fnobj, lower=dlbound, 
                               upper=dubound, method="L-BFGS-B")$par
        }
      }
    }
    
    modmat <- model.matrix(object=formula, data=d)
    for (q in (npf+1):(npf+no.ints)) {
      d[[q]] <- modmat[,which(attr(modmat, "assign")==q)]
      colnames(d[[q]]) <- NULL
      names(d)[q] <- attr(terms(formula), "term.labels")[q]
    }
    
    if (inter == 1) {
      z <- cbind(1, mlist(d,jmat))
    } else {
      z <- mlist(d,jmat)
    }
    final.opt <- fgeneralised(z=z, v=v, a=abc, w=wei)
    diff <- abs(curr.eval - final.opt)
    curr.eval <- final.opt
    
    if (progress == TRUE) {
      cat("Iteration", c(it,"( Current Value = ", curr.eval,")"), "\n")
    }
    it <- it+1
  }

  ptm <- proc.time()[3] - ptm
  
  output <- list(objval=final.opt, design=d[][1:npf], nits=it, time=ptm)
  output
}


