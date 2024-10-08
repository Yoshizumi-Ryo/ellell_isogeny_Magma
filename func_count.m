

function A_theta_transform(lv2tc,zeta_4)
    lv2tc[[1,0]]*:=zeta_4;
    lv2tc[[0,1]]*:=zeta_4;
    return lv2tc;
end function;



function lv2tc_to_frac(lv2tc)
    frac_lv2tc:=[lv2tc[[0,0]],lv2tc[[1,0]],lv2tc[[0,1]],lv2tc[[1,1]],1];
    return frac_lv2tc;
end function;




function Count_prepare_1(p,l)
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
    lv2_e1  :=prod_lv2tc(E_0,E_0,P_A    ,E_0!0  ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    lv2_e2  :=prod_lv2tc(E_0,E_0,E_0!0  ,P_A    ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    lv2_e12 :=prod_lv2tc(E_0,E_0,P_A    ,P_A    ,lv2tnp_lmd,sqrt_lm,sqrt_lmm1,lv2tnp_lmd,sqrt_lm,sqrt_lmm1);
    //theta transformation.
    zeta_4:=zeta_8^2;
    //assert(zeta_4^2 eq -1);
    tt_lv2_0   :=A_theta_transform(lv2_0   ,zeta_4);
    tt_lv2_e1  :=A_theta_transform(lv2_e1  ,zeta_4);
    tt_lv2_e2  :=A_theta_transform(lv2_e2  ,zeta_4);
    tt_lv2_e12 :=A_theta_transform(lv2_e12 ,zeta_4);
    //change to fractions.
    tt_lv2_0   :=lv2tc_to_frac(tt_lv2_0   );
    tt_lv2_e1  :=lv2tc_to_frac(tt_lv2_e1  );
    tt_lv2_e2  :=lv2tc_to_frac(tt_lv2_e2  );
    tt_lv2_e12 :=lv2tc_to_frac(tt_lv2_e12 );
    return tt_lv2_0,tt_lv2_e1,tt_lv2_e2,tt_lv2_e12;
end function;





procedure Time_for_isogeny_1(num_samples)
    p:=102282731615980594196068022554337214440490321465300216431359430531731842529319;
    ell_list:=[5,7,11,13];
    K:=GF(p^2);
    //p_list:= [2^4*3^58*5*29-1, 2^4*3^62*7-1, 2^4*3^55*11-1, 2^4*3^57*13-1, 2^4*3^54*17-1, 2^4*3^61*19-1];
    for i in [1..#ell_list] do
        l:=ell_list[i];
        //p:=p_list[i];
        times_codone:=[];
        times_codsq :=[];
        for j in [1..num_samples] do
            tc_0,tc_e1,tc_e2,tc_e12:=Count_prepare_1(p,l);
            //codsq.
            time_start:= Cputime();
            codomain_tc_0:=CodSq(tc_0,[tc_e1,tc_e2,tc_e12],l)[1];
            time_fin:=Cputime(time_start);
            Append(~times_codsq,time_fin);
            //codone.
            time_start:= Cputime();
            codomain_tc_0,lmd_data_1212:=CodOne(tc_0,[tc_e1,tc_e2,tc_e12],l);
            time_fin:=Cputime(time_start);
            Append(~times_codone,time_fin);
        end for;
        average_time_codone:=(&+times_codone)/num_samples;
        average_time_codsq :=(&+times_codsq) /num_samples;
        log2p:=Floor(Log(p) / Log(2)) + 1;
        //print(kind);
        printf "log_2(p)= %o\n", log2p;
        printf "ell= %o\n", l;
        printf "Samples: %o\n", num_samples;
        printf "CodSq : Average time(sec): %o\n", average_time_codsq;
        printf "CodOne: Average time(sec): %o\n", average_time_codone;
        print ("");
    end for;
end procedure;





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
    x_1a:=Random(E_0);
    x_1b:=Random(E_0);
    x_2a:=Random(E_0);
    x_2b:=Random(E_0);
    x_3a:=Random(E_0);
    x_3b:=Random(E_0);
    pts:=[];
    for i in [1..24] do
        pts[i]:=Random(E_0);
    end for;
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




procedure Time_for_isogeny_2(num_samples)
    //assert (num_points ge 1);
    p:=102282731615980594196068022554337214440490321465300216431359430531731842529319;
    ell_list:=[5,7,11,13];
    K:=GF(p^2);
    //p_list:= [2^4*3^58*5*29-1, 2^4*3^62*7-1, 2^4*3^55*11-1, 2^4*3^57*13-1, 2^4*3^54*17-1, 2^4*3^61*19-1];
    methods:=["CodOne","CodOne","CodOne","CodSq"];
    for i in [1..#ell_list] do
        method:=methods[i];
        l:=ell_list[i];
        //p:=p_list[i];
        times:=[];
        for k in [1..12] do
            times[k]:=[];
        end for;
        for j in [1..num_samples] do
            tc_0,tc_e1,tc_e2,tc_e12,tc_xs,tc_xpe1s,tc_xpe2s:=Count_prepare_2(p,l);
            //codone.
            time_start:= Cputime();
            if method eq "CodOne" then
                codomain_tc_0,lmd_data_1212:=CodOne(tc_0,[tc_e1,tc_e2,tc_e12],l);
            elif method eq "CodSq" then
                codomain_tc_0:=CodSq(tc_0,[tc_e1,tc_e2,tc_e12],l)[1];
            end if;
            lmd_data_pow12:=Product_power_lambda([tc_e1,tc_e2,tc_e12],l,lmd_data_1212);
            for k in [1..12] do
                img_tc:=EvalOne(tc_0,[tc_e1,tc_e2,tc_e12],[tc_xs[k],tc_xpe1s[k],tc_xpe2s[k]],lmd_data_pow12,l,K);
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
        print "Codomain:", method; 
        for k in [3,6,9,12] do
            printf "%o points, Average time(sec): %o\n", k,average_time[k];
        end for;
        print ("");
    end for;
end procedure;


