flm.ce.me <- function(d, t, n, npf, dx, knotsx, pars, db, knotsb, v, nb,
                  tol, dlbound, dubound, nx, flinear, progress, inter, ptm){
  
  tu <- max(t)
  eps <- tol 
  diff <- eps + 1
  
  jmat <- list()
  djmat <- list()
  for (i in 1:npf) {
    if (pars[i] == "power") {
      jmat[[i]] <- jcbfunccpp(db=db[i], degree=dx[i], knots=knotsx[[i]], t=t)
    } else {
      jmat[[i]] <- bs2jcbcpp(t=t, dbeta=db[i], dx=dx[i], knots_beta=knotsb[[i]], knots_x=knotsx[[i]])
    }
  }
  
  djmat <- mlist(d, jmat)
  if (inter == 1)  {
    z <- cbind(1, djmat)
  } else {
    z <- djmat
  }
  curr.eval <- flinear(z=z, v=v) 
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
              q <- cbind(1, mlist(f, jmat))
            } else {
              q <- mlist(f,jmat)
            }            
            result <- flinear(z=q, v=v)
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
    final.opt <- flinear(z=z, v=v)
    diff <- abs(curr.eval - final.opt)
    curr.eval <- final.opt
    
    if (progress == TRUE) {
      cat("Iteration", c(it,"( Current Value = ", curr.eval,")"), "\n")
    }
    it <- it + 1
  }
  
  ptm <- proc.time()[3] - ptm
  
  output <- list(objval=final.opt, design=d, nits=it, time=ptm)
  output
}



