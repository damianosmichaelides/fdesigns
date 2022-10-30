fglm.intspols <- function(d, tbounds, formula, n, npf, dx, knotsx, pars, db, knotsb, v, nb,
                          nx, abc, wei, fgeneralised, inter, ptm, tol, dlbound, dubound, progress,
                          no.terms, no.ints, no.pols=no.pols, criterion){
  
  tu <- max(tbounds)
  eps <- tol 
  diff <- eps + 1
  
  jmat <- list()
  djmat <- list()
  for (i in 1:npf) {
    if(pars[i] == "power"){
      jmat[[i]] <- jcbfunccpp(db=db[i], degree=dx[i], knots=knotsx[[i]], t=tbounds)
    } else{
      jmat[[i]] <- bs2jcbcpp(t=tbounds, dbeta=db[i], dx=dx[i], knots_beta=knotsb[[i]], knots_x=knotsx[[i]])
    }
  }
  
  modmat <- model.matrix(object=formula, data=d)
  modframe.intspols <- attr(terms(formula), "factors")[1:npf,]
  
  for (h in (npf+1):(npf+no.pols)) {
    mm.pol <- attr(terms(formula), "term.labels")[h]
    attach(d)
    P_a <- attr(eval(parse(text=mm.pol), envir=.GlobalEnv), "x")
    P_b <- attr(eval(parse(text=mm.pol), envir=.GlobalEnv), "deg")
    detach(d)
    for (p in 1:npf) {
      if (npf == 1) {
        modframe.intspols[h] <- ifelse(P_a == names(modframe.intspols)[p], P_b, 0)
      } else {
        if (npf > 1) {
          modframe.intspols[p,h] <- ifelse(P_a == rownames(modframe.intspols)[p], P_b, 0)
        }
      }
    }
    d[[h]] <- modmat[,which(attr(modmat, "assign")==h)]
    colnames(d[[h]]) <- NULL
    names(d)[h] <- attr(terms(formula), "term.labels")[h]
    if (npf == 1) {
      dx.pols <- modframe.intspols[h]
    } else {
      dx.pols <- modframe.intspols[,h]
    }
    dxpols.coord <- which(dx.pols > 0)
    allknotsxb <- replicate(dx.pols[dxpols.coord], unlist(knotsx[dxpols.coord]), simplify=FALSE)
    nx.ints <- rep(nx[dxpols.coord], dx.pols[dxpols.coord])
    if (prod(nx.ints) < nb[h] & identical(criterion, "A")) stop("For A-optimality, the product of basis functions for the profile factors involved in a polynomial must be greater or equal to the number of basis functions for the parameter")
    lkxi <- length(allknotsxb)
    lnxvec <- list()
    for (u in 1:lkxi) {
      lnxvec[[u]] <- 1:nx.ints[u]
    } 
    aposx <- as.matrix(expand.grid(lnxvec))
    if (pars[h] == "power") {
      jmat[[h]] <- bsalljcbpowercpp(t=tbounds, j=0, nbs=lkxi, aposx=aposx, alldeg=c(rep(dx[dxpols.coord], dx.pols[dxpols.coord])), allknotsx=allknotsxb, nx=nx.ints, nb=nb[h])
    } else {
      allknotsxb[[lkxi+1]] <- knotsb[[h]]
      jmat[[h]] <- bsalljcbcpp(t=tbounds, j=0, nbs=lkxi+1, aposx=aposx, alldeg=c(rep(dx[dxpols.coord], dx.pols[dxpols.coord]),db[h]), allknotsxb=allknotsxb, nx= nx.ints,  nb=nb[h])
    }
  }
  
  if (no.ints > 0) {
    for (q in (npf+no.pols+1):(npf+no.pols+no.ints)) {
      d[[q]] <- modmat[,which(attr(modmat, "assign")==q)]
      colnames(d[[q]]) <- NULL
      names(d)[q] <- attr(terms(formula), "term.labels")[q]
      dx.ints <- modframe.intspols[,q]
      allknotsxb <- knotsx[dx.ints==1]
      nx.ints <- nx[dx.ints==1]
      if(prod(nx.ints) < nb[q] & identical(criterion, "A")) stop("For A-optimality, the product of basis functions for the profile factors involved in an interaction must be greater or equal to the number of basis functions for the parameter")
      lkxi <- length(allknotsxb)
      lnxvec <- list()
      for (m in 1:lkxi) {
        lnxvec[[m]] <- 1:nx.ints[m]
      }  
      aposx <- as.matrix(expand.grid(lnxvec))
      if(pars[q] == "power"){
        jmat[[q]] <- bsalljcbpowercpp(t=tbounds, j=0, nbs=lkxi, aposx=aposx, alldeg=dx[dx.ints==1], allknotsx=allknotsxb, nx=nx.ints, nb=nb[q])
      } else{
        allknotsxb[[lkxi+1]] <- knotsb[[q]]
        jmat[[q]] <- bsalljcbcpp(t=tbounds, j=0, nbs=lkxi+1, aposx=aposx, alldeg=c(dx[dx.ints==1],db[q]), allknotsxb=allknotsxb, nx= nx.ints,  nb=nb[q])
      }
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
  if(progress == TRUE){
    cat("Starting design", c("( Current Value =", curr.eval,")"), "\n")
  }
  
  while(diff > eps) {
    
    for(i in 1:n) {
      for (p in 1:npf) {
        for(j in 1:nx[p]) {
          fnobj <- function(k) {
            f <- d
            f[[p]][i,j] <- k
            modmat <- model.matrix(object=formula, data=f)
            for (w in (npf+1):(npf+no.pols+no.ints)) {
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
    for (q in (npf+1):(npf+no.pols+no.ints)) {
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


