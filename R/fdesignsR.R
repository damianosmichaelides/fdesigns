linearcpp.aopt <- function(z, v){
	.Call( "linearcppaopt", z, v, PACKAGE = "fdesigns" )
}

linearcpp.dopt <- function(z, v){
  .Call( "linearcppdopt", z, v, PACKAGE = "fdesigns" )
}

logisticcpp.dopt <- function(z, v, a, w){
  .Call( "logisticcppdopt", z, v, a, w, PACKAGE = "fdesigns" )
}

logisticcpp.aopt <- function(z, v, a, w){
  .Call( "logisticcppaopt", z, v, a, w, PACKAGE = "fdesigns" )
}

poissoncpp.dopt <- function(z, v, a, w){
  .Call( "poissoncppdopt", z, v, a, w, PACKAGE = "fdesigns" )
}

poissoncpp.aopt <- function(z, v, a, w){
  .Call( "poissoncppdopt", z, v, a, w, PACKAGE = "fdesigns" )
}

# bs2intcpp <- function(t, pbeta, px, dbeta, dx, ex_lambdabeta, ex_lambdax){
#   .Call( "bs2intcpp", t, pbeta, px, dbeta, dx, ex_lambdabeta, ex_lambdax, PACKAGE = "funcmods" )
# }
# 
# bs2jcbcpp <- function(t, dbeta, dx, knots_beta, knots_x){
#   .Call( "bs2jcbcpp", t, dbeta, dx, knots_beta, knots_x, PACKAGE = "funcmods" )
# }
# 
# bsintcpp <- function(j, p, degree, exl, t){
#   .Call( "bsintcpp", j, p, degree, exl, t, PACKAGE = "funcmods" )
# }
# 
# jcbfunccpp <- function(db, degree, knots, t){
#   .Call( "jcbfunccpp", db, degree, knots, t, PACKAGE = "funcmods" )
# }
# 
# Vpowercpp <- function(db, tu){
#   .Call( "Vpowercpp", db, tu, PACKAGE = "funcmods" )
# }
# 
# bsallintcpp <- function(t, j, nbs, allpos, alldeg, exknots){
#   .Call( "bsallintcpp", t, j, nbs, allpos, alldeg, exknots, PACKAGE = "funcmods" )
# }
# 
# bsplinesecderprodcpp <- function(t, degree, pos, knots){
#   .Call( "bsplinesecderprodcpp", t, degree, pos, knots, PACKAGE = "funcmods" )
# }
# 
# Vbsplinecpp <- function(db, tu, knots){
#   .Call( "Vbsplinecpp", db, tu, knots, PACKAGE = "funcmods" )
# }
# 
# bsalljcbcpp <- function(t, j, nbs, aposx, alldeg, allknotsxb, nx, nb){
#   .Call( "bsalljcbcpp", t, j, nbs, aposx, alldeg, allknotsxb, nx, nb, PACKAGE = "funcmods" )
# }
# 
# bsalljcbpowercpp <- function(t, j, nbs, aposx, alldegx, allknotsx, nx, nb){
#   .Call( "bsalljcbpowercpp", t, j, nbs, aposx, alldegx, allknotsx, nx, nb, PACKAGE = "funcmods" )
# }
