

function A_theta_transform(lv2tc,zeta_4)
    lv2tc[[1,0]]*:=zeta_4;
    lv2tc[[0,1]]*:=zeta_4;
    return lv2tc;
end function;



function lv2tc_to_frac(lv2tc)
    frac_lv2tc:=[lv2tc[[0,0]],lv2tc[[1,0]],lv2tc[[0,1]],lv2tc[[1,1]],1];
    return frac_lv2tc;
end function;




function Count_prepare(p,l)
    N_B:=2;
    //assert(IsPrime(p));
    //assert((p mod 4) eq 3);
    //assert(IsOdd(l));
    //assert(IsPrime(l));
    //assert(IsDivisibleBy(p+1,l));
    //assert(IsDivisibleBy(p+1,N_B));
    K:=GF(p^2);
    //the 8th primitive root of 1.====================
    _<x>:=PolynomialRing(GF(p^2));
    zeta_8:=RootsInSplittingField(x^4+1)[1][1];
    //=======================================
    //public construction.=================
    //E_0: y^2=x^3-x=x(x-1)(x+1).
    _<x>:=PolynomialRing(GF(p^2));
    lmd:=K!(-1);
    E_0:=EllipticCurve(x*(x-1)*(x-lmd));
    //assert(IsSupersingular(E_0));
    //assert(#E_0 eq (p+1)^2);
    //take one basis of E_0[N_A]=(Z/N_A Z)^2.
    P_A,Q_A:=ell_to_torsion_basis_2(E_0,l);
    //========================================
    lv2tnp_lmd,sqrt_lm,sqrt_lmm1:=LegendreEll_to_lv2tnp(lmd);
    lv2_0:=product_theta(lv2tnp_lmd,lv2tnp_lmd);
    //e_1=(P_A,0). e_2=(0,P_A).
    lv2_e1 :=prod_lv2tc(E_0,E_0,P_A,E_0!0,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    lv2_e2 :=prod_lv2tc(E_0,E_0,E_0!0,P_A,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    lv2_e12:=prod_lv2tc(E_0,E_0,P_A,P_A  ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    //theta transformation.
    zeta_4:=zeta_8^2;
    //assert(zeta_4^2 eq -1);
    tt_lv2_0  :=A_theta_transform(lv2_0  ,zeta_4);
    tt_lv2_e1 :=A_theta_transform(lv2_e1 ,zeta_4);
    tt_lv2_e2 :=A_theta_transform(lv2_e2 ,zeta_4);
    tt_lv2_e12:=A_theta_transform(lv2_e12,zeta_4);
    //change to fractions.
    tt_lv2_0  :=lv2tc_to_frac(tt_lv2_0);
    tt_lv2_e1 :=lv2tc_to_frac(tt_lv2_e1);
    tt_lv2_e2 :=lv2tc_to_frac(tt_lv2_e2);
    tt_lv2_e12:=lv2tc_to_frac(tt_lv2_e12);
    return tt_lv2_0,tt_lv2_e1,tt_lv2_e2,tt_lv2_e12;
end function;





procedure Time_for_isogeny(num_samples,kind)
    p:=1536977231315969305425444639443786099;
    ell_list:=[3,5,7,11,13,17,19];
    //p_list:= [2^4*3^58*5*29-1, 2^4*3^62*7-1, 2^4*3^55*11-1, 2^4*3^57*13-1, 2^4*3^54*17-1, 2^4*3^61*19-1];
    for i in [1..#ell_list] do
        l:=ell_list[i];
        //p:=p_list[i];
        times:=[];
        for j in [1..num_samples] do
            tc_0,tc_e1,tc_e2,tc_e12:=Count_prepare(p,l);
            if kind eq "CodOne" then
                time_start:= Cputime();
                codomain_tc_0:=CodOne(tc_0,[tc_e1,tc_e2,tc_e12],l);
                time_fin:=Cputime(time_start);
            end if;
            if kind eq "CodSq" then
                time_start:= Cputime();
                codomain_tc_0:=CodSq(tc_0,[tc_e1,tc_e2,tc_e12],l);
                time_fin:=Cputime(time_start);
            end if;
            Append(~times,time_fin);
            //print("Sample:",j,". Time(sec):",time_fin);
        end for;
        average_time:=(&+times)/num_samples;
        log2p:=Floor(Log(p) / Log(2)) + 1;
        //print(kind);
        printf "ell= %o\n", l;
        //printf "log_2(p)= %o\n", log2p;
        printf "Average time(sec): %o\n", average_time;
        print ("");
    end for;
end procedure;



