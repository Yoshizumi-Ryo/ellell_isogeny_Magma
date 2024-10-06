



//use this.
function ell_to_torsion_basis_2(E,N)
  //assert(IsSquare(#E));
  _,sqrt_orderE:=IsSquare(#E);
  //assert(sqrt_orderE^2 eq #E);
  //assert(IsDivisibleBy(sqrt_orderE,N));
  _,r:=IsDivisibleBy(sqrt_orderE,N);
  //assert(#Generators(E) eq 2); 
  basis1:=Generators(E)[1];
  basis2:=Generators(E)[2];
  //assert(Order(basis1) eq sqrt_orderE);
  //assert(Order(basis2) eq sqrt_orderE);
  P:=r*basis1;
  Q:=r*basis2;
  //assert(Order(P) eq N);
  //assert(Order(Q) eq N);
  return P,Q;
end function;





//from y^2=x(x-1)(x-lm) to level 2 theta null point.
function LegendreEll_to_lv2tnp(lm)
  _<x>:=PolynomialRing(Parent(lm));
  sqrt_lm:=RootsInSplittingField(x^2-lm)[1][1];
  sqrt_lmm1:=RootsInSplittingField(x^2-(lm-1))[1][1];
  thnp0_sq:=sqrt_lm;
  thnp1_sq:=sqrt_lmm1;
  thnp2_sq:=1;
  thnp3_sq:=0;
  lv2tnp:=AssociativeArray();
  lv2tnp0:=thnp0_sq+thnp2_sq;
  lv2tnp1:=thnp1_sq;
  lv2tnp:=[lv2tnp0,lv2tnp1];
  return lv2tnp,sqrt_lm,sqrt_lmm1;
end function;


//level 2 theta null point of elliptic product.
function product_theta(lv2tc_1,lv2tc_2)
  lv2tc:=AssociativeArray();
  lv2tc[[0,0]]:=lv2tc_1[1]*lv2tc_2[1];
  lv2tc[[0,1]]:=lv2tc_1[1]*lv2tc_2[2];
  lv2tc[[1,0]]:=lv2tc_1[2]*lv2tc_2[1];
  lv2tc[[1,1]]:=lv2tc_1[2]*lv2tc_2[2];
  return lv2tc;
end function;




//the inverse function.
function lv2tnp_to_LegendreEll(lv2tnp)
  lm:=((lv2tnp[1]^2+lv2tnp[2]^2)/(lv2tnp[1]^2-lv2tnp[2]^2))^2;
  r:=lv2tnp[1]/lv2tnp[2];
  sqrt_lm:=((lm-1)*r^2-lm-1)/2;
  sqrt_lmm1:=(sqrt_lm+1)/r;
  //assert(sqrt_lm^2 eq lm);
  //assert(sqrt_lmm1^2 eq lm-1);
  return lm,sqrt_lm,sqrt_lmm1;
end function;


//for P in E, we give level 2 theta coordinate.
function Legendre_to_lv2tc(lv2tnp,sqrt_lm,sqrt_lmm1,P)
  if P eq Parent(P)!0 then
    return lv2tnp;
  end if;
  u:=P[1];
  lm:=lv2tnp_to_LegendreEll(lv2tnp);
  //assert(sqrt_lm^2 eq lm);
  //assert(sqrt_lmm1^2 eq lm-1);
  thc0_sq:=sqrt_lm*(u-1);
  thc1_sq:=sqrt_lmm1*u;
  thc2_sq:=sqrt_lm*thc0_sq-sqrt_lmm1*thc1_sq;
  thc3_sq:=sqrt_lm*thc1_sq-sqrt_lmm1*thc0_sq;
  return [1,(thc1_sq+thc3_sq)/(thc0_sq+thc2_sq)];
end function;



function prod_lv2tc(E1,E2,S1,S2,lv2tnp_1,sqrt_lm_1,sqrt_lmm1_1,lv2tnp_2,sqrt_lm_2,sqrt_lmm1_2)
  //assert(S1 in E1);
  //assert(S2 in E2);
  lv2_S1:=Legendre_to_lv2tc(lv2tnp_1,sqrt_lm_1,sqrt_lmm1_1,S1);
  lv2_S2:=Legendre_to_lv2tc(lv2tnp_2,sqrt_lm_2,sqrt_lmm1_2,S2);
  lv2_S:=product_theta(lv2_S1,lv2_S2);
  return lv2_S;
end function;

