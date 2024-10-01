



function SumOfSquare(l)
    assert IsPrime(l);
    assert l ne 2;
    if l eq 3 then
        return <1, 1, 1>;
    end if;
    if l mod 4 eq 1 then
        for a in [0..l] do
            is_square, b := IsSquare(l - a^2);
            if is_square then
                return <a, b>;
            end if;
        end for;
    end if;
    if l mod 8 eq 3 then
        for a in [1..l] do
            is_square, b := IsSquare(l - 2 * a^2);
            if is_square then
                return <a, a, b>;
            end if;
        end for;
    end if;
    if l mod 8 eq 7 and l mod 3 eq 1 then
        for a in [1..l] do
            is_square, b := IsSquare(l - 3 * a^2);
            if is_square then
                return <a, a, a, b>;
            end if;
        end for;
    end if;
    for a in [0..l] do
        for b in [a..l] do
            for c in [b..l] do
                is_square, d := IsSquare(l - a^2 - b^2 - c^2);
                if is_square then
                    if a eq 0 then
                        assert b ne 0;
                        return <b, c, d>;
                    else
                        return <a, b, c, d>;
                    end if;
                end if;
            end for;
        end for;
    end for;
end function;




function HalfCoeffWithout0(l)
    ld := (l - 1) div 2;
    coeff_set := { <1, 0>, <0, 1>, <1, 1> };
    for k2 in [2..ld] do
        Include(~coeff_set, <0, k2>);
    end for;
    for k2 in [2..l-2] do
        Include(~coeff_set, <1, k2>);
    end for;
    for k1 in [2..ld] do
        Include(~coeff_set, <k1, 0>);
    end for;
    for k1 in [2..l-1] do
        Include(~coeff_set, <k1, 1>);
    end for;
    for k1 in [2..l-2] do
        Include(~coeff_set, <k1, 2>);
    end for;
    for k1 in [2..ld] do
        for k2 in [3..l - k1 - 1] do
            Include(~coeff_set, <k1, k2>);
        end for;
    end for;
    for k1 in [ld + 1..l - 1] do
        for k2 in [3..l - k1] do
            Include(~coeff_set, <k1, k2>);
        end for;
    end for;
    assert #coeff_set eq (l^2 - 1) div 2;
    return coeff_set;
end function;





function HalfLinCom(tc_0, tc_e1, tc_e2, tc_e12, l)
    assert IsOdd(l);
    assert IsPrime(l);
    ld := (l - 1) div 2;
    lincom := AssociativeArray();
    lincom[<0,0>] := tc_0;
    lincom[<1,0>] := tc_e1;
    lincom[<0,1>] := tc_e2;
    lincom[<1,1>] := tc_e12;
    for k2 in [2..ld] do
        lincom[<0, k2>] := Diff_Add(tc_0,lincom[<0, k2-1>], tc_e2, lincom[<0, k2-2>]);
    end for;
    for k2 in [2..l-2] do
        lincom[<1, k2>] := Diff_Add(tc_0,lincom[<1, k2-1>], tc_e2, lincom[<1, k2-2>]);
    end for;
    for k1 in [2..ld] do
        lincom[<k1, 0>] := Diff_Add(tc_0,lincom[<k1-1, 0>], tc_e1, lincom[<k1-2, 0>]);
    end for;
    for k1 in [2..l-1] do
        lincom[<k1, 1>] := Diff_Add(tc_0,lincom[<k1-1, 1>], tc_e1, lincom[<k1-2, 1>]);
    end for;
    for k1 in [2..l-2] do
        lincom[<k1, 2>] := Diff_Add(tc_0,lincom[<k1-1, 2>], tc_e1, lincom[<k1-2, 2>]);
    end for;
    for k1 in [2..ld] do
        for k2 in [3..l-k1-1] do
            lincom[<k1, k2>] := Diff_Add(tc_0,lincom[<k1, k2-1>], tc_e2, lincom[<k1, k2-2>]);
        end for;
    end for;
    for k1 in [ld+1..l-1] do
        for k2 in [3..l-k1] do
            lincom[<k1, k2>] := Diff_Add(tc_0,lincom[<k1, k2-1>], tc_e2, lincom[<k1, k2-2>]);
        end for;
    end for;
    Remove(~lincom, <0, 0>);
    assert #lincom eq (l^2 - 1) div 2;
    assert HalfCoeffWithout0(l) eq Keys(lincom);
    return lincom;
end function;



function Remain_Half_coeff_without0(l)
    remained_set := { <k1, k2> : k1 in [0..l-1], k2 in [0..l-1] };
    remained_set := remained_set diff HalfCoeffWithout0(l);
    remained_set := remained_set diff { <0, 0> };
    assert 2 * #remained_set + 1 eq l^2;
    return remained_set;
end function;





function LmdLpow(tc_0, tc_e, tc_lde, tc_ldm1e)
    tc_ldp1e := Diff_Add(tc_0, tc_lde, tc_e, tc_ldm1e);
    numer_1 := tc_lde[1];  
    denom_1 := tc_ldp1e[5];   
    numer_2 := tc_ldp1e[1];  
    denom_2 := tc_lde[5];     
    lmd_lpow_value := [numer_1 * denom_1, numer_2 * denom_2];
    return lmd_lpow_value;
end function;




// Construct k_1 * e_1 + k_2 * e_2 for 0 <= k_1, k_2 < l
function XpLinCom(tc_0, basis, pt_list,l,K)
    tc_e1   := basis[1];
    tc_e2   := basis[2];
    tc_e12  := basis[3];
    tc_x    := pt_list[1];
    tc_xpe1 := pt_list[2];
    tc_xpe2 := pt_list[3];
    // Compute tc_xpe12
    tc_xpe12 := Extended_Addition(tc_0, tc_x, tc_e1, tc_e2, tc_xpe1, tc_e12, tc_xpe2,K);
    // Initialize dictionary of theta coordinates
    xplincom := AssociativeArray();
    xplincom[<0, 0>] := tc_x;
    xplincom[<1, 0>] := tc_xpe1;
    xplincom[<0, 1>] := tc_xpe2;
    xplincom[<1, 1>] := tc_xpe12;
    // Compute values for (0, k2) and (1, k2)
    for k2 in [2..l-1] do
        xplincom[<0, k2>] := Diff_Add(tc_0,xplincom[<0, k2-1>], tc_e2, xplincom[<0, k2-2>]);
        xplincom[<1, k2>] := Diff_Add(tc_0,xplincom[<1, k2-1>], tc_e2, xplincom[<1, k2-2>]);
    end for;
    // Compute values for (k1, k2)
    for k2 in [0..l-1] do
        for k1 in [2..l-1] do
            assert not IsDefined(xplincom, <k1, k2>);
            xplincom[<k1, k2>] := Diff_Add(tc_0,xplincom[<k1-1, k2>], tc_e1, xplincom[<k1-2, k2>]);
        end for;
    end for;
    // Check that all values are computed
    assert #Keys(xplincom) eq l^2;   
    return xplincom;
end function;



//=============================
//=============================
//=============================



function CodomainCommon(tc_0, basis,l)
    tc_e1  := basis[1];
    tc_e2  := basis[2];
    tc_e12 := basis[3];
    h_lincom := HalfLinCom(tc_0, tc_e1, tc_e2, tc_e12, l);
    assert IsPrime(l);
    assert #h_lincom eq (l^2 - 1) div 2;
    ld := (l - 1) div 2;
    h_lincom[<0, 0>] := tc_0;
    lmd1_lpow  := LmdLpow(tc_0, tc_e1 , h_lincom[<ld,  0>], h_lincom[<ld-1, 0   >]);
    lmd2_lpow  := LmdLpow(tc_0, tc_e2 , h_lincom[<0 , ld>], h_lincom[<0   , ld-1>]);
    lmd12_lpow := LmdLpow(tc_0, tc_e12, h_lincom[<ld, ld>], h_lincom[<ld-1, ld-1>]);
    lmd_data_1212:=[lmd1_lpow,lmd2_lpow,lmd12_lpow];
    Remove(~h_lincom, <0, 0>);
    lmd_div_lpow := [
        lmd12_lpow[1] * lmd1_lpow[2] * lmd2_lpow[2],
        lmd12_lpow[2] * lmd1_lpow[1] * lmd2_lpow[1]
    ];
    lmds_with_common_denom, den := Common_denom_frac_3([lmd1_lpow, lmd2_lpow, lmd_div_lpow]);
    return l, h_lincom, lmds_with_common_denom[1], lmds_with_common_denom[2], lmds_with_common_denom[3], den,lmd_data_1212;
end function;










function CodOne(tc_0, basis,l)
    l, h_lincom, lmd1_lpow, lmd2_lpow, lmd_div_lpow, den,lmd_data_1212 := CodomainCommon(tc_0, basis,l);
    den_pow := Multpower_straight(den, 3 * (l - 1)^2);
    lmd1_lpow_pow :=  Multpower_sq(lmd1_lpow[1], l);
    lmd2_lpow_pow :=  Multpower_sq(lmd2_lpow[1], l);
    lmd_div_lpow_pow := Multpower_straight(lmd_div_lpow[1], (l - 1)^2);
    h_lincom_lpow := AssociativeArray();
    for key in Keys(h_lincom) do
        h_lincom_lpow[key] := [h_lincom[key][i]^l : i in [1..5]];
    end for;
    assert #h_lincom_lpow eq (l^2 - 1) div 2;
    for key in Keys(h_lincom) do
        k1 := key[1];
        k2 := key[2];
        coeff := [
            lmd1_lpow_pow[k1^2] * lmd2_lpow_pow[k2^2] * lmd_div_lpow_pow[k1 * k2],
            den_pow[k1^2 + k2^2 + k1 * k2]
        ];
        h_lincom_lpow[key] := Mult_frac(h_lincom_lpow[key],coeff);
    end for;
    assert #h_lincom_lpow eq (l^2 - 1) div 2;
    tc_00 := [tc_0[i]^l : i in [1..5]];
    h_lincom_lpow[<0, 0>] := tc_00;
    h_lincom_lpow, _ := Dict_common_denom(h_lincom_lpow);
    assert #h_lincom_lpow eq (l^2 + 1) div 2;
    tc_f0 := h_lincom_lpow[<0, 0>];
    assert #HalfCoeffWithout0(l) eq (l^2 - 1) div 2;
    for key in HalfCoeffWithout0(l) do
        for i in [1..4] do
            tc_f0[i] +:= 2 * h_lincom_lpow[key][i];
        end for;
    end for;
    return tc_f0,lmd_data_1212;
end function;





function CodSq(tc_0, basis,l)
    l, lincom, lmd1_lpow, lmd2_lpow, lmd_div_lpow, den,lmd_data_1212 := CodomainCommon(tc_0, basis,l);
    a_u := SumOfSquare(l);
    r := #a_u;    
    ss1 := [(a_u[u]*(l-1)) mod l : u in [1..r]];
    tt1 := [(a_u[u]*(l-1) - ss1[u]) div l : u in [1..r]];
    max_exp := (l-1)^2 + l * (&+[tt1[u]^2 : u in [1..r]]) - 2 * (l-1) * (&+[a_u[u]*tt1[u] : u in [1..r]]);
    assert max_exp le r*(l-1);
    lmd1_lpow_pow    := Multpower_straight(lmd1_lpow[1], max_exp);
    lmd2_lpow_pow    := Multpower_straight(lmd2_lpow[1], max_exp);
    lmd_div_lpow_pow := Multpower_straight(lmd_div_lpow[1], max_exp);
    den_pow          := Multpower_straight(den, max_exp * 3); 
    assert #lincom eq (l^2 - 1) div 2;
    lincom[<0, 0>] := tc_0;
    assert #lincom eq (l^2 + 1) div 2;
    for k1k2 in Remain_Half_coeff_without0(l) do
        k1:=k1k2[1];
        k2:=k1k2[2];
        a1 := l - 2 * k1;
        a2 := l - 2 * k2;
        adiv := l - k1 - k2;
        aden := a1 + a2 + adiv;
        assert aden eq 3 * adiv;
        if k1 eq 0 then
            assert a2 le 0;
            lincom[<0, k2>] := Mult_frac(lincom[<0, l - k2>] , [den_pow[-a2], lmd2_lpow_pow[-a2]]);
        elif k2 eq 0 then
            assert a1 le 0;
            lincom[<k1, 0>] := Mult_frac(lincom[<l - k1, 0>] , [den_pow[-a1], lmd1_lpow_pow[-a1]]);
        elif (a1 ge 0) and (a2 ge 0) then
            lincom[<k1, k2>] := Mult_frac(lincom[<l - k1, l - k2>] , [lmd1_lpow_pow[a1] * lmd2_lpow_pow[a2] * den_pow[-aden], lmd_div_lpow_pow[-adiv]]);
        elif (a1 ge 0) and (a2 lt 0) then
            lincom[<k1, k2>] := Mult_frac(lincom[<l - k1, l - k2>] , [lmd1_lpow_pow[a1] * den_pow[-aden], lmd2_lpow_pow[-a2] * lmd_div_lpow_pow[-adiv]]);
        elif (a1 lt 0) and (a2 ge 0) then
            lincom[<k1, k2>] := Mult_frac(lincom[<l - k1, l - k2>] , [lmd2_lpow_pow[a2] * den_pow[-aden], lmd1_lpow_pow[-a1] * lmd_div_lpow_pow[-adiv]]);
        else
            assert a1 lt 0 and a2 lt 0;
            lincom[<k1, k2>] := Mult_frac(lincom[<l - k1, l - k2>] , [den_pow[-aden], lmd1_lpow_pow[-a1] * lmd2_lpow_pow[-a2] * lmd_div_lpow_pow[-adiv]]);
        end if;
    end for;
    assert #lincom eq l^2;
    pre_tc_f0 :=[];
    denom:=tc_0[5]^r;
    for i in [1..4] do
        pre_tc_f0[i] := [&*[tc_0[(a_u[u] mod 2) * i+1] : u in [1..r]],denom];
    end for;
    for k1k2 in HalfCoeffWithout0(l) do
        k1:=k1k2[1];
        k2:=k1k2[2];
        s1 := [(a_u[u] * k1) mod l : u in [1..r]];
        t1 := [(a_u[u] * k1 - s1[u]) div l : u in [1..r]];
        s2 := [(a_u[u] * k2) mod l : u in [1..r]];
        t2 := [(a_u[u] * k2 - s2[u]) div l : u in [1..r]];
        h_1 := k1^2 + l * (&+[t1[u]^2 : u in [1..r]]) - 2*k1 * (&+[a_u[u]*t1[u] : u in [1..r]]);
        h_2 := k2^2 + l * (&+[t2[u]^2 : u in [1..r]]) - 2*k2 * (&+[a_u[u]*t2[u] : u in [1..r]]);
        h_div := k1*k2 + l * (&+[t1[u]*t2[u] : u in [1..r]]) - k1 * (&+[a_u[u]*t2[u] : u in [1..r]]) - k2 * (&+[a_u[u]*t1[u] : u in [1..r]]);
        h_den := h_1 + h_2 + h_div;
        coeff := [lmd1_lpow_pow[h_1] * lmd2_lpow_pow[h_2] * lmd_div_lpow_pow[h_div], den_pow[h_den]];
        for i in [1..4] do
            plus_term_num := coeff[1] * (&* [lincom[<s1[u], s2[u]>][(a_u[u] mod 2) * i+1] : u in [1..r]]);
            plus_term_den := coeff[2] * (&* [lincom[<s1[u], s2[u]>][5]                  : u in [1..r]]);
            pre_tc_f0[i] := Frac_add(pre_tc_f0[i], [2 * plus_term_num, plus_term_den]);
        end for;
    end for;
    proj_tc_f0 := Projective_Theta(pre_tc_f0);
    tc_f0 := [proj_tc_f0[1],proj_tc_f0[2],proj_tc_f0[3],proj_tc_f0[4],1];
    return tc_f0,lmd_data_1212;
end function;



//====================================
//====================================
//====================================

// Product_power_lambda function in Magma
function Product_power_lambda(basis,l,lmd_data_1212)
    assert #basis eq 3;
    tc_e1  := basis[1];
    tc_e2  := basis[2];
    tc_e12 := basis[3];
    // Extract lambda power values
    lmd1_lpow  := lmd_data_1212[1];
    lmd2_lpow  := lmd_data_1212[2];
    lmd12_lpow := lmd_data_1212[3];
    // Compute lambda division power
    lmd_div_lpow := [lmd12_lpow[1] * lmd1_lpow[2] * lmd2_lpow[2], lmd12_lpow[2] * lmd1_lpow[1] * lmd2_lpow[1]];
    // Common denominator fraction
    lmd_seq, den := Common_denom_frac_3([lmd1_lpow, lmd2_lpow, lmd_div_lpow]);
    lmd1_lpow   :=lmd_seq[1];
    lmd2_lpow   :=lmd_seq[2];
    lmd_div_lpow:=lmd_seq[3];
    // Compute den_pow and lambda powers
    den_pow          := Multpower_straight(den,            3 * (l - 1)^2);
    lmd1_lpow_pow    := Multpower_straight(lmd1_lpow[1],   l^2);
    lmd2_lpow_pow    := Multpower_straight(lmd2_lpow[1],   l^2);
    lmd_div_lpow_pow := Multpower_straight(lmd_div_lpow[1],(l - 1)^2);
    // Initialize product dictionary
    lmd_pow_product := AssociativeArray();
    // Compute lambda power products
    for k1 in [0..l-1] do
        for k2 in [0..l-1] do
            lmd_pow_product[<k1, k2>] := [lmd1_lpow_pow[k1^2] * lmd2_lpow_pow[k2^2] * lmd_div_lpow_pow[k1*k2],den_pow[k1^2+k2^2+k1*k2]];
        end for;
    end for;
    // Compute lambda square for l
    lm_1_lsq := [lmd1_lpow_pow[l], den_pow[l]];
    lm_2_lsq := [lmd2_lpow_pow[l], den_pow[l]];
    // Ensure the length of the product dictionary is correct
    assert #Keys(lmd_pow_product) eq l^2;
    lmd_data_pow12:=<lmd_pow_product, lm_1_lsq, lm_2_lsq>;
    return lmd_data_pow12;
end function;






// Take_lmd_power function in Magma
function Take_lmd_power(tc_x, xple1, xple2, lm_1_lsq, lm_2_lsq, l)
    // Compute ml1_lpow and ml2_lpow
    ml1_lpow := [tc_x[1] * lm_1_lsq[2] * xple1[5], tc_x[5] * lm_1_lsq[1] * xple1[1]];
    ml2_lpow := [tc_x[1] * lm_2_lsq[2] * xple2[5], tc_x[5] * lm_2_lsq[1] * xple2[1]];
    // Common denominator fraction
    mu_seq, m_den := Common_denom_frac_2([ml1_lpow, ml2_lpow]);
    ml1_lpow:=mu_seq[1];
    ml2_lpow:=mu_seq[2];
    // Compute powers
    ml1_lpow_pow   := Multpower_straight(ml1_lpow[1], (l - 1));
    ml2_lpow_pow   := Multpower_straight(ml2_lpow[1], (l - 1));
    muden_lpow_pow := Multpower_straight(m_den, 2 * (l - 1));
    return ml1_lpow_pow, ml2_lpow_pow, muden_lpow_pow;
end function;







// EvalOne function in Magma
function EvalOne(tc_0, basis, pt_list,lmd_data_pow12,l,K)
    lmd_pow_product :=lmd_data_pow12[1];
    lm_1_lsq        :=lmd_data_pow12[2];
    lm_2_lsq        :=lmd_data_pow12[3];
    // Assertions to check the structure of inputs
    assert Type(lmd_pow_product) eq Assoc;
    assert Type(lm_1_lsq) eq SeqEnum;
    assert Type(lm_2_lsq) eq SeqEnum;
    assert #basis eq 3;
    assert #pt_list eq 3;
    assert #Keys(lmd_pow_product) eq l^2;
    // Compute linear combination of x and basis
    xplincom := XpLinCom(tc_0, basis, pt_list,l,K);
    assert #Keys(xplincom) eq l^2;
    tc_e1 := basis[1];
    tc_e2 := basis[2];
    tc_x  := pt_list[1];
    // Compute xple1 and xple2 using Diff_Add
    xple1 := Diff_Add(tc_0,xplincom[<l-1,0>], tc_e1, xplincom[<l-2,0>]);
    xple2 := Diff_Add(tc_0,xplincom[<0,l-1>], tc_e2, xplincom[<0,l-2>]);
    // Take lambda power
    ml1_lpow_pow, ml2_lpow_pow, muden_lpow_pow := Take_lmd_power(tc_x, xple1, xple2, lm_1_lsq, lm_2_lsq, l);
    assert #Keys(xplincom) eq l^2;
    // Normalize xplincom using common denominator
    xplincom, den := Dict_common_denom(xplincom);
    den_lpow := den^l;
    // Compute xplincom_lpow
    xplincom_lpow := AssociativeArray();
    for key in Keys(xplincom) do
        xplincom_lpow[key] :=[];
        for i in [1..4] do
            xplincom_lpow[key][i] := xplincom[key][i]^l;
        end for;
        xplincom_lpow[key][5]:=den_lpow;
    end for;
    assert #Keys(xplincom_lpow) eq l^2;
    // Modify xplincom_lpow with coefficients
    for key in Keys(xplincom_lpow) do
        k1, k2 := Explode(key);
        coeff := [lmd_pow_product[<k1,k2>][1] * ml1_lpow_pow[k1] * ml2_lpow_pow[k2],
                  lmd_pow_product[<k1,k2>][2] * muden_lpow_pow[k1+k2]];
        xplincom_lpow[key] := Mult_frac(xplincom_lpow[key],coeff);
    end for;
    assert #Keys(xplincom_lpow) eq l^2;
    // Normalize again with common denominator
    xplincom_lpow, den_lpow := Dict_common_denom(xplincom_lpow);
    pre_tc_fx := [K!0, K!0, K!0, K!0];
    // Sum up the numerators for final result
    for k1 in [0..l-1] do
        for k2 in [0..l-1] do
            assert xplincom_lpow[<k1,k2>][5] eq den_lpow;
            for i in [1..4] do
                pre_tc_fx[i] +:= xplincom_lpow[<k1,k2>][i];
            end for;
        end for;
    end for;
    tc_fx := [pre_tc_fx[1],pre_tc_fx[2],pre_tc_fx[3],pre_tc_fx[4],1];
    return tc_fx;
end function;






