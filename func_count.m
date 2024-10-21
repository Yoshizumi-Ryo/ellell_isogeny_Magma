

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
    P_A,Q_A:=ell_to_torsion_basis_2(E_0,3*l);
    x1:=Random(E_0);
    x2:=Random(E_0);
    //========================================
    lv2tnp_lmd,sqrt_lm,sqrt_lmm1:=LegendreEll_to_lv2tnp(lmd);
    lv2_0:=product_theta(lv2tnp_lmd,lv2tnp_lmd);
    //e_1=(P_A,0). e_2=(0,P_A).
    lv2_e1  :=prod_lv2tc(E_0,E_0,P_A    ,E_0!0  ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    lv2_e2  :=prod_lv2tc(E_0,E_0,E_0!0  ,P_A    ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    lv2_e12 :=prod_lv2tc(E_0,E_0,P_A    ,P_A    ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);

    lv2_x   :=prod_lv2tc(E_0,E_0,x1    ,x2   ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    lv2_xpe1   :=prod_lv2tc(E_0,E_0,x1+P_A    ,x2   ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    lv2_xpe2   :=prod_lv2tc(E_0,E_0,x1    ,x2+P_A   ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    lv2_xpe12   :=prod_lv2tc(E_0,E_0,x1+P_A    ,x2+P_A   ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);


    //theta transformation.
    zeta_4:=zeta_8^2;
    //assert(zeta_4^2 eq -1);
    tt_lv2_0   :=A_theta_transform(lv2_0   ,zeta_4);
    tt_lv2_e1  :=A_theta_transform(lv2_e1  ,zeta_4);
    tt_lv2_e2  :=A_theta_transform(lv2_e2  ,zeta_4);
    tt_lv2_e12 :=A_theta_transform(lv2_e12 ,zeta_4);
    tt_lv2_x   :=A_theta_transform(lv2_x   ,zeta_4);
    tt_lv2_xpe1:=A_theta_transform(lv2_xpe1,zeta_4);
    tt_lv2_xpe2:=A_theta_transform(lv2_xpe2,zeta_4);
    tt_lv2_xpe12:=A_theta_transform(lv2_xpe12,zeta_4);

    //change to fractions.----------
    tt_lv2_0   :=lv2tc_to_frac(tt_lv2_0   );
    tt_lv2_e1  :=lv2tc_to_frac(tt_lv2_e1  );
    tt_lv2_e2  :=lv2tc_to_frac(tt_lv2_e2  );
    tt_lv2_e12 :=lv2tc_to_frac(tt_lv2_e12 );
    tt_lv2_x   :=lv2tc_to_frac(tt_lv2_x );
    tt_lv2_xpe1 :=lv2tc_to_frac(tt_lv2_xpe1 );
    tt_lv2_xpe2 :=lv2tc_to_frac(tt_lv2_xpe2 );
    tt_lv2_xpe12 :=lv2tc_to_frac(tt_lv2_xpe12 );

    //prepare.--------
    lv2_le1  :=Mult(tt_lv2_0,tt_lv2_e1,l);
    lv2_le2  :=Mult(tt_lv2_0,tt_lv2_e2,l);
    lv2_le12 :=Mult(tt_lv2_0,tt_lv2_e12,l);

    lv2_e1ple1 :=Mult(tt_lv2_0,tt_lv2_e1,l+1);
    lv2_e1ple2 :=Mxpy(tt_lv2_0,l,tt_lv2_e2,tt_lv2_e1,tt_lv2_e12);
    lv2_e1ple12:=Mxpy(tt_lv2_0,l,lv2_le12,tt_lv2_e1,Mxpy(tt_lv2_0,2,tt_lv2_e1,tt_lv2_e2,tt_lv2_e12));

    lv2_e2ple1 :=Mxpy(tt_lv2_0,l,tt_lv2_e1,tt_lv2_e2,tt_lv2_e12);
    lv2_e2ple2 :=Mult(tt_lv2_0,tt_lv2_e2,l+1);
    lv2_e2ple12:=Mxpy(tt_lv2_0,l,lv2_le12,tt_lv2_e2,Mxpy(tt_lv2_0,2,tt_lv2_e2,tt_lv2_e1,tt_lv2_e12));
    
    lv2_e12ple1 :=Mxpy(tt_lv2_0,l+1,tt_lv2_e1,tt_lv2_e2,tt_lv2_e12);
    lv2_e12ple2 :=Mxpy(tt_lv2_0,l+1,tt_lv2_e2,tt_lv2_e1,tt_lv2_e12);
    lv2_e12ple12 :=Mult(tt_lv2_0,tt_lv2_e12,l+1);

    lv2_xple1 :=Mxpy(tt_lv2_0,l,tt_lv2_e1 ,tt_lv2_x,tt_lv2_xpe1);
    lv2_xple2 :=Mxpy(tt_lv2_0,l,tt_lv2_e2 ,tt_lv2_x,tt_lv2_xpe2);
    lv2_xple12 :=Mxpy(tt_lv2_0,l,tt_lv2_e12 ,tt_lv2_x,tt_lv2_xpe12);

    lv2_xpe1ple1 :=Mxpy(tt_lv2_0,l+1,tt_lv2_e1 ,tt_lv2_x,tt_lv2_xpe1);
    lv2_xpe1ple2 :=Mxpy(tt_lv2_0,l,tt_lv2_e2 ,tt_lv2_xpe1,tt_lv2_xpe12);
    lv2_xpe1ple12 :=Mxpy(tt_lv2_0,l,tt_lv2_e12 ,tt_lv2_xpe1,Diff_Add(tt_lv2_0,tt_lv2_xpe12,tt_lv2_e1,lv2_xple2));

    lv2_xpe2ple1 :=Mxpy(tt_lv2_0,l,tt_lv2_e2 ,tt_lv2_xpe1,tt_lv2_xpe12);
    lv2_xpe2ple2 :=Mxpy(tt_lv2_0,l+1,tt_lv2_e2 ,tt_lv2_x,tt_lv2_xpe2);
    lv2_xpe2ple12 :=Mxpy(tt_lv2_0,l,tt_lv2_e12 ,tt_lv2_xpe2,Diff_Add(tt_lv2_0,tt_lv2_xpe12,tt_lv2_e2,lv2_xple1));
    //---------
    codomain_tc_0,lmd_data_1212:=CodSq(tt_lv2_0,[lv2_le1,lv2_le2,lv2_le12],3);
    lmd_data_pow12:=Product_power_lambda([lv2_le1,lv2_le2,lv2_le12],3,lmd_data_1212);

    img_e1:=EvalOne(tt_lv2_0,[lv2_le1,lv2_le2,lv2_le12],[lv2_e1ple1,lv2_e1ple2,lv2_e1ple12],lmd_data_pow12,3,K);

    img_e2:=EvalOne(tt_lv2_0,[lv2_le1,lv2_le2,lv2_le12],[lv2_e2ple1,lv2_e2ple2,lv2_e2ple12],lmd_data_pow12,3,K);

    img_e12:=EvalOne(tt_lv2_0,[lv2_le1,lv2_le2,lv2_le12],[lv2_e12ple1,lv2_e12ple2,lv2_e12ple12],lmd_data_pow12,3,K);

    img_x  :=EvalOne(tt_lv2_0,[lv2_le1,lv2_le2,lv2_le12],[lv2_xple1,lv2_xple2,lv2_xple12],lmd_data_pow12,3,K);

    img_xpe1:=EvalOne(tt_lv2_0,[lv2_le1,lv2_le2,lv2_le12],[lv2_xpe1ple1,lv2_xpe1ple2,lv2_xpe1ple12],lmd_data_pow12,3,K);

    img_xpe2:=EvalOne(tt_lv2_0,[lv2_le1,lv2_le2,lv2_le12],[lv2_xpe2ple1,lv2_xpe2ple2,lv2_xpe2ple12],lmd_data_pow12,3,K);

    return codomain_tc_0,img_e1,img_e2,img_e12,img_x,img_xpe1,img_xpe2;
end function;





procedure Time_for_isogeny_1(l,num_samples,method)
    p:=102282731615980594196068022554337214440490321465300216431359430531731842529319;
    K:=GF(p^2);
    //p_list:= [2^4*3^58*5*29-1, 2^4*3^62*7-1, 2^4*3^55*11-1, 2^4*3^57*13-1, 2^4*3^54*17-1, 2^4*3^61*19-1];
        //p:=p_list[i];
    times:=[];
    for j in [1..num_samples] do
        tc_0,tc_e1,tc_e2,tc_e12,_,_,_:=Count_prepare(p,l);
        //codsq.
        time_start:= Cputime();
        if method eq "CodSq" then
            codomain_tc_0:=CodSq(tc_0,[tc_e1,tc_e2,tc_e12],l)[1];
        elif method eq "CodOne" then
            codomain_tc_0,lmd_data_1212:=CodOne(tc_0,[tc_e1,tc_e2,tc_e12],l);
        end if;
        time_fin:=Cputime(time_start);
        Append(~times,time_fin);
    end for;
    average_time:=(&+times)/num_samples;
    log2p:=Floor(Log(p) / Log(2)) + 1;
    //print(kind);
    printf "log_2(p)= %o\n", log2p;
    printf "ell= %o\n", l;
    printf "Samples: %o\n", num_samples;
    printf "method: %o\n", method;
    printf "Average time(sec): %o\n", average_time;
end procedure;




/*
function Count_prepare_2(p,l)
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
    x_1:=Random(E_0);
    x_2:=Random(E_0);
    //========================================
    lv2tnp_lmd,sqrt_lm,sqrt_lmm1:=LegendreEll_to_lv2tnp(lmd);
    lv2_0:=product_theta(lv2tnp_lmd,lv2tnp_lmd);
    //e_1=(P_A,0). e_2=(0,P_A).
    lv2_e1 :=prod_lv2tc(E_0,E_0,P_A     ,E_0!0   ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    lv2_e2 :=prod_lv2tc(E_0,E_0,E_0!0   ,P_A     ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    lv2_e12:=prod_lv2tc(E_0,E_0,P_A     ,P_A     ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    lv2_xs   :=[];
    lv2_xpe1s:=[];
    lv2_xpe2s:=[];
    for i in [1..12] do
        lv2_xs   [i]:=prod_lv2tc(E_0,E_0,pts[i]    ,pts[i+12]    ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
        lv2_xpe1s[i]:=prod_lv2tc(E_0,E_0,pts[i]+P_A,pts[i+12]    ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
        lv2_xpe2s[i]:=prod_lv2tc(E_0,E_0,pts[i]    ,pts[i+12]+P_A,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    end for;
    zeta_4:=zeta_8^2;
    //assert(zeta_4^2 eq -1);
    tt_lv2_0  :=A_theta_transform(lv2_0    ,zeta_4);
    tt_lv2_e1 :=A_theta_transform(lv2_e1   ,zeta_4);
    tt_lv2_e2 :=A_theta_transform(lv2_e2   ,zeta_4);
    tt_lv2_e12:=A_theta_transform(lv2_e12  ,zeta_4);
    tt_lv2_x1s   :=[];
    tt_lv2_x1pe1s:=[];
    tt_lv2_x1pe2s:=[];
    for i in [1..12] do
        tt_lv2_x1s   [i]:=A_theta_transform(lv2_xs   [i],zeta_4);
        tt_lv2_x1pe1s[i]:=A_theta_transform(lv2_xpe1s[i],zeta_4);
        tt_lv2_x1pe2s[i]:=A_theta_transform(lv2_xpe2s[i],zeta_4);
    end for;
    //change to fractions.
    tt_lv2_0    :=lv2tc_to_frac(tt_lv2_0    );
    tt_lv2_e1   :=lv2tc_to_frac(tt_lv2_e1   );
    tt_lv2_e2   :=lv2tc_to_frac(tt_lv2_e2   );
    tt_lv2_e12  :=lv2tc_to_frac(tt_lv2_e12  );
    f_lv2_x1s   :=[];
    f_lv2_x1pe1s:=[];
    f_lv2_x1pe2s:=[];
    for i in [1..12] do
        f_lv2_x1s   [i]:=lv2tc_to_frac(tt_lv2_x1s   [i]);
        f_lv2_x1pe1s[i]:=lv2tc_to_frac(tt_lv2_x1pe1s[i]);
        f_lv2_x1pe2s[i]:=lv2tc_to_frac(tt_lv2_x1pe2s[i]);
    end for;
    return tt_lv2_0,tt_lv2_e1,tt_lv2_e2,tt_lv2_e12,f_lv2_x1s,f_lv2_x1pe1s,f_lv2_x1pe2s;
end function;
*/



procedure Time_for_isogeny_2(l,num_samples,method)
    //assert (num_points ge 1);
    p:=102282731615980594196068022554337214440490321465300216431359430531731842529319;
    K:=GF(p^2);
    //p_list:= [2^4*3^58*5*29-1, 2^4*3^62*7-1, 2^4*3^55*11-1, 2^4*3^57*13-1, 2^4*3^54*17-1, 2^4*3^61*19-1];
    //p:=p_list[i];
    times:=[];
    for k in [1..12] do
        times[k]:=[];
    end for;
    for j in [1..num_samples] do
        tc_0,tc_e1,tc_e2,tc_e12,tc_x,tc_xpe1,tc_xpe2:=Count_prepare(p,l);
        //codone.
        time_start:= Cputime();
        if method eq "CodOne" then
            codomain_tc_0,lmd_data_1212:=CodOne(tc_0,[tc_e1,tc_e2,tc_e12],l);
        elif method eq "CodSq" then
            codomain_tc_0,lmd_data_1212:=CodSq(tc_0,[tc_e1,tc_e2,tc_e12],l);
        end if;
        lmd_data_pow12:=Product_power_lambda([tc_e1,tc_e2,tc_e12],l,lmd_data_1212);
        for k in [1..12] do
            img_tc:=EvalOne(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data_pow12,l,K);
            time_end:=Cputime(time_start);
            Append(~times[k],time_end);
        end for;
    end for;
    average_time:=[];
    for k in [1..12] do
        average_time[k]:=(&+times[k])/num_samples;
    end for;
    log2p:=Floor(Log(p) / Log(2)) + 1;
    //print(kind);
    printf "log_2(p)= %o\n", log2p;
    printf "ell= %o\n", l;
    printf "Samples: %o\n", num_samples;
    printf "Codomain: %o\n", method; 
    for k in [3,6,9,12] do
        printf "%o points, Average time(sec): %o\n", k,average_time[k];
    end for;  
end procedure;


