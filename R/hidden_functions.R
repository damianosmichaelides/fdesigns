mlist <- function(l1, l2){
  ll1 <- length(l1)
  # ll2 <- length(l2) - ll1 and ll2 should be the same anyway
  l12 <- list()
  for (i in 1:ll1) {
    l12[[i]] <- l1[[i]] %*% l2[[i]]
  }
  matrix(unlist(l12), nrow=nrow(l1[[1]]))
}

basis <- function(x, degree, i, knots) {
  if(degree == 0){
    B <- ifelse((x >= knots[i]) & (x < knots[i+1]), 1, 0)
  } else {
    if((knots[degree+i] - knots[i]) == 0) {
      alpha1 <- 0
    } else {
      alpha1 <- (x - knots[i])/(knots[degree+i] - knots[i])
    }
    if((knots[i+degree+1] - knots[i+1]) == 0) {
      alpha2 <- 0
    } else {
      alpha2 <- (knots[i+degree+1] - x)/(knots[i+degree+1] - knots[i+1])
    }
    B <- alpha1*basis(x, (degree-1), i, knots) + alpha2*basis(x, (degree-1), (i+1), knots)
  }
  return(B)
}

mybs <- function(x, degree=3, interior.knots=NULL, intercept=FALSE, Boundary.knots = c(0,1)) {
  if(missing(x)) stop("You must provide x")
  Boundary.knots <- sort(Boundary.knots)
  interior.knots.sorted <- NULL
  if(!is.null(interior.knots)) interior.knots.sorted <- sort(interior.knots)
  knots <- c(rep(Boundary.knots[1], (degree+1)), interior.knots.sorted, rep(Boundary.knots[2], (degree+1)))
  K <- length(interior.knots) + degree + 1
  B.mat <- matrix(0,length(x),K)
  for(j in 1:K) B.mat[,j] <- basis(x, degree, j, knots)
  if(any(x == Boundary.knots[2])) B.mat[x == Boundary.knots[2], K] <- 1 # why
  if(intercept == FALSE) {
    return(B.mat[,-1])
  } else {
    return(B.mat)
  }
}



# Vpower <- function(db, tu){
#   nb <- db + 1
#   Vpow <- matrix(0, nrow=nb, ncol=nb)
#   for (b in 1:nb){
#     for (p in 1:nb){
#       Vpow[b,p] <- (tu^(b+(p-1))) / (b+(p-1))
#     }
#   }
#   Vpow
# }

# Vpower <- function(db, tu){
#   nb <- db + 1
#   Vpow <- matrix(0, nrow=nb, ncol=nb)
#   for (b in 0:(nb-1)){
#     for (p in 0:(nb-1)){
#       Vpow[b+1,p+1] <- ifelse( (b*p*(b-1)*(p-1))==0, 0, b*p*(b-1)*(p-1) * (tu^(b+p-3)) / (b+p-3) )
#     }
#   }
#   Vpow
# }

# bspline.secder.prod <- function(t, pos, deg, knots){
#   "/" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y))
#   
#   if(deg<2){
#     result <- 0
#   } else{
#     
#     dt2.deg <- deg - 2
#     
#     exknot <- c(rep(t[1],deg+1), knots, rep(max(t),deg+1))
#     
#     pos1 <- pos[1]
#     pos2 <- pos[2]
#     
#     A <- exknot[pos1+deg] - exknot[pos1]
#     A1 <- exknot[pos1+deg-1] - exknot[pos1]
#     A2 <- exknot[pos1+deg] - exknot[pos1+1]
#     B <- exknot[pos1+deg+1] - exknot[pos1+1]
#     B1 <- exknot[pos1+deg] - exknot[pos1+1]
#     B2 <- exknot[pos1+deg+1] - exknot[pos1+2]
#     
#     C <- exknot[pos2+deg] - exknot[pos2]
#     C1 <- exknot[pos2+deg-1] - exknot[pos2]
#     C2 <- exknot[pos2+deg] - exknot[pos2+1]
#     D <- exknot[pos2+deg+1] - exknot[pos2+1]
#     D1 <- exknot[pos2+deg] - exknot[pos2+1]
#     D2 <- exknot[pos2+deg+1] - exknot[pos2+2]
#     
#     AA1 <- A*A1
#     AA2 <- A*A2
#     BB1 <- B*B1
#     BB2 <- B*B2
#     CC1 <- C*C1
#     CC2 <- C*C2
#     DD1 <- D*D1
#     DD2 <- D*D2
#     
#     alld <- rep(dt2.deg,2)
#     lexk <- list(exknot,exknot)
#     
#     int1 <- bsallintcpp(t=t, j=0, nbs=2, allpos=c(pos1,pos2), alldeg=alld, exknots=lexk)
#     int2 <- bsallintcpp(t=t, j=0, nbs=2, allpos=c(pos1,pos2+1), alldeg=alld, exknots=lexk)
#     int3 <- bsallintcpp(t=t, j=0, nbs=2, allpos=c(pos1,pos2+2), alldeg=alld, exknots=lexk)
#     int4 <- bsallintcpp(t=t, j=0, nbs=2, allpos=c(pos1+1,pos2), alldeg=alld, exknots=lexk)
#     int5 <- bsallintcpp(t=t, j=0, nbs=2, allpos=c(pos1+1,pos2+1), alldeg=alld, exknots=lexk)
#     int6 <- bsallintcpp(t=t, j=0, nbs=2, allpos=c(pos1+1,pos2+2), alldeg=alld, exknots=lexk)
#     int7 <- bsallintcpp(t=t, j=0, nbs=2, allpos=c(pos1+2,pos2), alldeg=alld, exknots=lexk)
#     int8 <- bsallintcpp(t=t, j=0, nbs=2, allpos=c(pos1+2,pos2+1), alldeg=alld, exknots=lexk)
#     int9 <- bsallintcpp(t=t, j=0, nbs=2, allpos=c(pos1+2,pos2+2), alldeg=alld, exknots=lexk)
#     
#     dconst <- (deg * (deg-1))^2
#     
#     result <- dconst * (
#       int1 / (AA1*CC1)    -   (DD1+CC2) * int2 / (AA1*CC2*DD1)    +     int3 / (AA1*DD2)  -
#         (BB1+AA2) * int4 / (AA2*BB1*CC1)    +     (BB1+AA2) * (DD1+CC2) * int5 / (AA2*BB1*CC2*DD1)    -
#           (BB1+AA2) * int6 / (AA2*BB1*DD2)    +     int7 / (BB2*CC1)    -     
#             (DD1+CC2) * int8 / (BB2*CC2*DD1)    +     int9 / (BB2*DD2)
#     )
#     
#     result
#   }
# }
# 
# Vbspline <- function(db, tu, knots){
#   nb <- db + length(knots) + 1
#   Vpow <- matrix(0, nrow=nb, ncol=nb)
#   for (b in 1:nb){
#     for (p in 1:nb){
#       Vpow[b,p] <- bspline.secder.prod(t=c(0,tu), pos=c(b,p), deg=db, knots=knots)
#     }
#   }
#   Vpow
# }

# system.time(Vbspline(db=3, tu=1, knots=c(0.2,0.4,0.6,0.8)))
# system.time(Vbsplinecpp(db=3, tu=1, knots=c(0.2,0.4,0.6,0.8)))
# 
# Vbspline(db=5, tu=2, knots=c(0.24,0.43,0.62,0.84)) - Vbsplinecpp(db=5, tu=2, knots=c(0.24,0.43,0.62,0.84))
# 
# system.time(bsplinesecderprodcpp(t=c(0,1), degree=5, pos=c(5,3), knots=c(0.2,0.4,0.6,0.8)))
# system.time(bspline.secder.prod(t=c(0,1), deg=5, pos=c(5,3), knots=c(0.2,0.4,0.6,0.8)))


# bs2jcbRcpp <- function(t, dbeta, dx, knots_beta, knots_x){
#   bs2jcbcpp(t=t, dbeta=dbeta, dx=dx, knots_beta=as.numeric(knots_beta), 
#             knots_x=as.numeric(knots_x))
# }
# 
# 
# jcbfuncRcpp <- function(db, degree, knots, t){
#   jcbfunccpp(db=db, degree=degree, knots=as.numeric(knots), t=t)
# }
# 
# 
# bsalljcbpowerRcpp <- function(t, j, nbs, aposx, alldegx, allknotsx, nx, nb){
#         bsalljcbpowercpp(t=t, j=j, nbs=nbs, aposx=aposx, alldegx=alldegx, 
#                   allknotsx=as.numeric(allknotsx), nx=nx, nb=nb)
# }
# 
# bsalljcbRcpp <- function(t, j, nbs, aposx, alldeg, allknotsxb, nx, nb){
#         bsalljcbcpp(t=t, j=j, nbs=nbs, aposx=aposx, alldeg=alldeg, 
#                   allknotsxb=as.numeric(allknotsxb), nx=nx, nb=nb)
# }







