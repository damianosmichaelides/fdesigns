#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export]]
double divide (double x, double y) {
  if(y==0){
    return 0;
  } else {
    return x/y ; }}




// [[Rcpp::export]]
double bs2intcpp(int j, NumericVector t, int pbeta, int px, int dbeta, int dx, NumericVector ex_lambdabeta, NumericVector ex_lambdax) {
  
  if(dbeta==0 & dx==0) {
    NumericVector lowerlimits {ex_lambdabeta(pbeta-1), ex_lambdax(px-1)};
    NumericVector upperlimits {ex_lambdabeta(pbeta), ex_lambdax(px)};
    double llim = max(lowerlimits);
    double ulim = min(upperlimits);
    if(ulim <= llim){
      return 0;
    } else{
      return divide( pow(ulim,(j+1)) - pow(llim,(j+1)), (j+1) );
    }}
  
  if(dbeta==0 & dx>0)
    return divide( bs2intcpp(j+1, t, pbeta, px, dbeta, dx-1, ex_lambdabeta, ex_lambdax), (ex_lambdax(px+dx-1)-ex_lambdax(px-1)) ) - 
      (ex_lambdax(px-1)) * divide( bs2intcpp(j, t, pbeta, px, dbeta, dx-1, ex_lambdabeta, ex_lambdax), (ex_lambdax(px+dx-1)-ex_lambdax(px-1)) ) + 
      (ex_lambdax(px+dx)) * divide( bs2intcpp(j, t, pbeta, px+1, dbeta, dx-1, ex_lambdabeta, ex_lambdax), (ex_lambdax(px+dx)-ex_lambdax(px)) ) - 
      divide( bs2intcpp(j+1, t, pbeta, px+1, dbeta, dx-1, ex_lambdabeta, ex_lambdax), (ex_lambdax(px+dx)-ex_lambdax(px)) ) ;
  
  if(dbeta>0 & dx==0){
    return divide( bs2intcpp(j+1, t, pbeta, px, dbeta-1, dx,  ex_lambdabeta, ex_lambdax), (ex_lambdabeta(pbeta+dbeta-1)-ex_lambdabeta(pbeta-1)) ) - 
      (ex_lambdabeta(pbeta-1)) * divide( bs2intcpp(j, t, pbeta, px, dbeta-1, dx, ex_lambdabeta, ex_lambdax), (ex_lambdabeta(pbeta+dbeta-1)-ex_lambdabeta(pbeta-1)) ) + 
      (ex_lambdabeta(pbeta+dbeta)) * divide( bs2intcpp(j, t, pbeta+1, px, dbeta-1, dx, ex_lambdabeta, ex_lambdax), (ex_lambdabeta(pbeta+dbeta)-ex_lambdabeta(pbeta)) ) - 
      divide( bs2intcpp(j+1, t, pbeta+1, px, dbeta-1, dx, ex_lambdabeta, ex_lambdax), (ex_lambdabeta(pbeta+dbeta)-ex_lambdabeta(pbeta)) ) ;
  } else{
    double A = (ex_lambdax(px+dx-1) - ex_lambdax(px-1)) * (ex_lambdabeta(pbeta+dbeta-1) - ex_lambdabeta(pbeta-1));
    double B = (ex_lambdax(px+dx-1) - ex_lambdax(px-1)) * (ex_lambdabeta(pbeta+dbeta) - ex_lambdabeta(pbeta));
    double C = (ex_lambdax(px+dx) - ex_lambdax(px)) * (ex_lambdabeta(pbeta+dbeta-1) - ex_lambdabeta(pbeta-1));
    double D = (ex_lambdax(px+dx) - ex_lambdax(px)) * (ex_lambdabeta(pbeta+dbeta) - ex_lambdabeta(pbeta));
    
    return 
      
      divide(1,A) * bs2intcpp(j+2, t, pbeta, px, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) - 
        divide(ex_lambdabeta(pbeta-1),A) * bs2intcpp(j+1, t, pbeta, px, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) + 
        divide(ex_lambdabeta(pbeta+dbeta),B) * bs2intcpp(j+1, t, pbeta+1, px, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) -
        divide(1,B) * bs2intcpp(j+2, t, pbeta+1, px, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) -
        
        divide(ex_lambdax(px-1),A) * bs2intcpp(j+1, t, pbeta, px, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) + 
        divide(ex_lambdax(px-1) * ex_lambdabeta(pbeta-1),A) * bs2intcpp(j, t, pbeta, px, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) - 
        divide(ex_lambdax(px-1) * ex_lambdabeta(pbeta+dbeta),B) * bs2intcpp(j, t, pbeta+1, px, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) + 
        divide(ex_lambdax(px-1),B) * bs2intcpp(j+1, t, pbeta+1, px, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) +
        
        divide(ex_lambdax(px+dx),C) * bs2intcpp(j+1, t, pbeta, px+1, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) - 
        divide(ex_lambdax(px+dx) * ex_lambdabeta(pbeta-1),C) * bs2intcpp(j, t, pbeta, px+1, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) +
        divide(ex_lambdax(px+dx) * ex_lambdabeta(pbeta+dbeta),D) * bs2intcpp(j, t, pbeta+1, px+1, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) - 
        divide(ex_lambdax(px+dx),D) * bs2intcpp(j+1, t, pbeta+1, px+1, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) - 
        
        divide(1,C) * bs2intcpp(j+2, t, pbeta, px+1, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) +
        divide(ex_lambdabeta(pbeta-1),C) * bs2intcpp(j+1, t, pbeta, px+1, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) - 
        divide(ex_lambdabeta(pbeta+dbeta),D) * bs2intcpp(j+1, t, pbeta+1, px+1, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) + 
        divide(1,D) * bs2intcpp(j+2, t, pbeta+1, px+1, dbeta-1, dx-1, ex_lambdabeta, ex_lambdax) ;
  } 
}





// [[Rcpp::export]]
NumericMatrix bs2jcbcpp(NumericVector t, int dbeta, int dx, NumericVector knots_beta, NumericVector knots_x) {
  double tl = t.length();
  double klen_x = knots_x.length();
  double klen_beta = knots_beta.length();
  int nx = klen_x + dx + 1;
  int nb = klen_beta + dbeta + 1;
  
  NumericVector ext1_x (dx+1,t(0));
  NumericVector ext2_x (dx+1,t(tl-1));
  NumericVector exl_x (2*dx + 2 + klen_x);
  int m;
  int k;
  for(m=0; m<(dx+1); ++m) {
    exl_x(m) = ext1_x(m);
    exl_x(dx+1+klen_x+m) = ext2_x(m);}
  for(k=0; k<klen_x; ++k) {
    exl_x(k+dx+1) = knots_x(k);}
  
  NumericVector ext1_beta (dbeta+1,t(0));
  NumericVector ext2_beta (dbeta+1,t(tl-1));
  NumericVector exl_beta (2*dbeta + 2 + klen_beta);
  int q;
  int r;
  for(q=0; q<(dbeta+1); ++q) {
    exl_beta(q) = ext1_beta(q);
    exl_beta(dbeta+1+klen_beta+q) = ext2_beta(q);}
  for(r=0; r<klen_beta; ++r) {
    exl_beta(r+dbeta+1) = knots_beta(r);}
  
  NumericMatrix jcb (nx,nb);
  int i;
  int l;
  for(i=0; i<nx; ++i) {
    for(l=0; l<nb; ++l) {
      jcb(i,l) = bs2intcpp(0,t,l+1,i+1,dbeta,dx,exl_beta,exl_x);}}
  return jcb;
}





// [[Rcpp::export]]
double bsintcpp(int j,int p, int degree, NumericVector exl, NumericVector t ) {
  if (degree == 0) {
    return divide(pow(exl(p),j+1) - pow(exl(p-1),j+1), j+1);
  } else {
    return divide(bsintcpp(j+1,p,degree-1,exl,t),(exl(p+degree-1)-exl(p-1))) - 
      exl(p-1) * divide(bsintcpp(j,p,degree-1,exl,t),(exl(p+degree-1)-exl(p-1))) + 
      exl(p+degree) * divide(bsintcpp(j,p+1,degree-1,exl,t),(exl(p+degree)-exl(p))) - 
      divide(bsintcpp(j+1,p+1,degree-1,exl,t),(exl(p+degree)-exl(p))) ;}}




// [[Rcpp::export]]
NumericMatrix jcbfunccpp(int db, int degree, NumericVector knots, NumericVector t ) {
  int nb = db + 1;
  double tt = t.length();
  double klen = knots.length();
  int nx = degree + klen + 1;
  
  NumericVector ext1 (degree+1,t(0));
  NumericVector ext2 (degree+1,t(tt-1));
  NumericVector exl (2*degree + 2 + klen);
  
  int m;
  int k;
  
  for(m=0; m<(degree+1); ++m) {
    exl(m) = ext1(m);
    exl(degree+1+klen+m) = ext2(m);}
  
  for(k=0; k<klen; ++k) {
    exl(k+degree+1) = knots(k);}
  
  NumericMatrix jcb (nx,nb);
  int i;
  int j;
  for(i=0; i<nx; ++i) {
    for(j=0; j<nb; ++j) {
      jcb(i,j) = bsintcpp(j,i+1,degree,exl,t);}}
  return jcb;}




// [[Rcpp::export]]
NumericMatrix Vpowercpp(int db, int tu) {
  int nb = db + 1;
  
  NumericMatrix Vpow (nb,nb);
  
  int b;
  int p;
  
  for(b=0; b<nb; ++b) {
    for(p=0; p<nb; ++p) {
      if((b*p*(b-1)*(p-1))==0){
        Vpow(b,p) = 0;
      } else{
        Vpow(b,p) = (b*p*(b-1)*(p-1)) * pow(tu,(b+p-3)) / (b+p-3);
      }
      
    }
  }
  return Vpow;
}

  

  
// [[Rcpp::export]]
  double bsallintcpp(NumericVector t, double j, double nbs, NumericVector allpos, NumericVector alldeg, List exknots) {
    
    bool b1 = is_true(any(alldeg>0));
    
    if(b1 == 1){
      int i = which_max(alldeg);
      int deg = alldeg(i);
      int pos = allpos(i);
      NumericVector allposup = clone(allpos);
      NumericVector allupdateddeg = clone(alldeg);
      allposup(i) = allposup(i) + 1;
      int posup = pos + 1;
      allupdateddeg(i) = alldeg(i) - 1;
      NumericVector exk = exknots[i];
      
      return  
        divide( bsallintcpp(t, j+1, nbs, allpos, allupdateddeg, exknots), (exk(deg+pos-1) - exk(pos-1)) ) - 
          exk(pos-1) * divide( bsallintcpp(t, j, nbs, allpos, allupdateddeg, exknots), (exk(deg+pos-1) - exk(pos-1)) ) +
          exk(deg+posup-1) * divide( bsallintcpp(t, j, nbs, allposup, allupdateddeg, exknots), (exk(deg+posup-1) - exk(posup-1)) ) -
          divide( bsallintcpp(t, j+1, nbs, allposup, allupdateddeg, exknots), (exk(deg+posup-1) - exk(posup-1)) );
    } 
    
    if(b1 == 0){
      NumericVector llimv (nbs);
      NumericVector ulimv (nbs);
      int k;
      int ak;
      for (k=0; k<(nbs); ++k) {
        NumericVector exkk = exknots[k];
        ak = allpos(k);
        llimv(k) = exkk(ak-1);
        ulimv(k) = exkk(ak);
      }
      double llim = max(llimv);
      double ulim = min(ulimv);
      if(ulim <= llim){
        return 0;
      } else{
        return divide( pow(ulim,(j+1)) - pow(llim,(j+1)), (j+1) );
      }
    } else{return 1000000000;}
  }



  
// [[Rcpp::export]]
  double bsplinesecderprodcpp(NumericVector t, int degree, NumericVector pos, NumericVector knots) {
    
    double tu = max(t);
    double klen = knots.length();

    NumericVector ext1 (degree+1,t(0));
    NumericVector ext2 (degree+1,tu);
    NumericVector exl (2*degree + 2 + klen);
    
    int m;
    int k;
    
    for(m=0; m<(degree+1); ++m) {
      exl(m) = ext1(m);
      exl(degree+1+klen+m) = ext2(m);}
    
    for(k=0; k<klen; ++k) {
      exl(k+degree+1) = knots(k);}
    
    double dt2deg;
    
    if(degree<2){
      return(0);
    } else{
      dt2deg = degree - 2;
      
      double pos1 = pos(0);
      double pos2 = pos(1);
      
      double A = exl(pos1+degree-1) - exl(pos1-1);
      double A1 = exl(pos1+degree-2) - exl(pos1-1);
      double A2 = exl(pos1+degree-1) - exl(pos1);
      double B = exl(pos1+degree) - exl(pos1);
      double B1 = exl(pos1+degree-1) - exl(pos1);
      double B2 = exl(pos1+degree) - exl(pos1+1);
      
      double C = exl(pos2+degree-1) - exl(pos2-1);
      double C1 = exl(pos2+degree-2) - exl(pos2-1);
      double C2 = exl(pos2+degree-1) - exl(pos2);
      double D = exl(pos2+degree) - exl(pos2);
      double D1 = exl(pos2+degree-1) - exl(pos2);
      double D2 = exl(pos2+degree) - exl(pos2+1);
      
      double AA1 = A*A1;
      double AA2 = A*A2;
      double BB1 = B*B1;
      double BB2 = B*B2;
      double CC1 = C*C1;
      double CC2 = C*C2;
      double DD1 = D*D1;
      double DD2 = D*D2;
      
      NumericVector alld = {dt2deg,dt2deg};
      List lexk(2);
      
      int l;
      for(l=0; l<2; ++l) {
        lexk(l) = exl;}
      
      NumericVector apos1 = {pos1,pos2};
      NumericVector apos2 = {pos1,pos2+1};
      NumericVector apos3 = {pos1,pos2+2};
      NumericVector apos4 = {pos1+1,pos2};
      NumericVector apos5 = {pos1+1,pos2+1};
      NumericVector apos6 = {pos1+1,pos2+2};
      NumericVector apos7 = {pos1+2,pos2};
      NumericVector apos8 = {pos1+2,pos2+1};
      NumericVector apos9 = {pos1+2,pos2+2};
      
        
        double int1 = bsallintcpp(t, 0, 2, apos1, alld, lexk);
        double int2 = bsallintcpp(t, 0, 2, apos2, alld, lexk);
        double int3 = bsallintcpp(t, 0, 2, apos3, alld, lexk);
        double int4 = bsallintcpp(t, 0, 2, apos4, alld, lexk);
        double int5 = bsallintcpp(t, 0, 2, apos5, alld, lexk);
        double int6 = bsallintcpp(t, 0, 2, apos6, alld, lexk);
        double int7 = bsallintcpp(t, 0, 2, apos7, alld, lexk);
        double int8 = bsallintcpp(t, 0, 2, apos8, alld, lexk);
        double int9 = bsallintcpp(t, 0, 2, apos9, alld, lexk);
        
        double dconst = pow((degree * (degree-1)), 2);
          
          return
          
          dconst * (
          divide(int1,(AA1*CC1))    -   divide((DD1+CC2) * int2, (AA1*CC2*DD1))    +     divide(int3,  (AA1*DD2))  -
            divide((BB1+AA2) * int4, (AA2*BB1*CC1))    +     divide((BB1+AA2) * (DD1+CC2) * int5, (AA2*BB1*CC2*DD1))  -
            divide((BB1+AA2) * int6, (AA2*BB1*DD2))    +     divide(int7, (BB2*CC1))    -     
            divide((DD1+CC2) * int8, (BB2*CC2*DD1))    +     divide(int9, (BB2*DD2)) );
      }
  }





// [[Rcpp::export]]
NumericMatrix Vbsplinecpp(double db, double tu, NumericVector knots) {
  
  double klen = knots.length();
  double nb = db + 1 + klen;
  
  NumericMatrix Vpow (nb,nb);
  NumericVector tt = {0,tu};
  
  double b;
  double p;
  
  NumericVector pos2(2);
  
  for(b=0; b<nb; ++b) {
    for(p=0; p<nb; ++p) {
      pos2 = {b+1,p+1};
      Vpow(b,p) = bsplinesecderprodcpp(tt, db, pos2, knots);
    }
  }
  return Vpow;
}





// [[Rcpp::export]]
NumericMatrix bsalljcbcpp(NumericVector t, double j, double nbs, NumericMatrix aposx,
                          NumericVector alldeg, List allknotsxb, NumericVector nx, double nb) {
  
  int k;
  int l;
  int ll;
  int p;
  double nxprod = 1;
  double nxl = nx.length();
  NumericVector apvec(nbs-1);
  NumericVector apvl(nbs);
  
  for(k=0; k<(nxl); ++k) {
    nxprod = nxprod * nx(k);}
  
  double tl = t.length();
  double kl = allknotsxb.length();
  List extknotsxb(kl);
  
  int c;
  for(c=0; c<kl; ++c) {
    NumericVector knotsxb = allknotsxb[c];
    double kxl = knotsxb.length();
    int dxb = alldeg(c);
    NumericVector ext1_x (dxb+1,t(0));
    NumericVector ext2_x (dxb+1,t(tl-1));
    NumericVector exl_x (2*dxb + 2 + kxl);
    int m;
    int k;
    for(m=0; m<(dxb+1); ++m) {
      exl_x(m) = ext1_x(m);
      exl_x(dxb+1+kxl+m) = ext2_x(m);}
    for(k=0; k<kxl; ++k) {
      exl_x(k+dxb+1) = knotsxb(k);}
    extknotsxb[c] = exl_x;
  }
  
  NumericMatrix jcb(nxprod,nb);
  
  for(l=0; l<(nxprod); ++l) {
    apvec = aposx(l,_);
    for(ll=0; ll<(nbs-1); ++ll) {
      apvl(ll) = apvec(ll);}
    for(p=0; p<(nb); ++p) {
      apvl(nbs-1) = p+1;
      jcb(l,p) = bsallintcpp(t,0,nbs,apvl,alldeg,extknotsxb);}}
  
  return jcb;
}





// [[Rcpp::export]]
NumericMatrix bsalljcbpowercpp(NumericVector t, double j, double nbs, NumericMatrix aposx,
                               NumericVector alldegx, List allknotsx, NumericVector nx, double nb) {
  
  int k;
  int l;
  int p;
  double nxprod = 1;
  double nxl = nx.length();
  NumericVector apvec(nbs);
  
  for(k=0; k<(nxl); ++k) {
    nxprod = nxprod * nx(k);}
  
  double tl = t.length();
  double kl = allknotsx.length();
  List extknotsx(kl);
  
  int c;
  for(c=0; c<kl; ++c) {
    NumericVector knotsx = allknotsx[c];
    double kxl = knotsx.length();
    int dx = alldegx(c);
    NumericVector ext1_x (dx+1,t(0));
    NumericVector ext2_x (dx+1,t(tl-1));
    NumericVector exl_x (2*dx + 2 + kxl);
    int m;
    int k;
    for(m=0; m<(dx+1); ++m) {
      exl_x(m) = ext1_x(m);
      exl_x(dx+1+kxl+m) = ext2_x(m);}
    for(k=0; k<kxl; ++k) {
      exl_x(k+dx+1) = knotsx(k);}
    extknotsx[c] = exl_x;
  }
  
  NumericMatrix jcb(nxprod,nb);
  
  for(l=0; l<(nxprod); ++l) {
    apvec = aposx(l,_);
    for(p=0; p<(nb); ++p) {
      jcb(l,p) = bsallintcpp(t,p,nbs,apvec,alldegx,extknotsx);}}
  
  return jcb;
} 










