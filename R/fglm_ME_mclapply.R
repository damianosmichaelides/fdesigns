fglm.me <- function(d, tbounds, formula, n, npf, dx, knotsx, pars, db, knotsb, v, nb,
                    nx, abc, wei, fgeneralised, inter, ptm, tol, dlbound, dubound, progress){
  
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
    
    for (p in 1:npf) {
      for (i in 1:n) {
        for (j in 1:nx[p]) {
          fnobj <- function(k) {
            f <- d
            f[[p]][i,j] <- k
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
  
  output <- list(objval=final.opt, design=d, nits=it, time=ptm)
  output
}

