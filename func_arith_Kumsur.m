


function Hadamard(tcn)
    Htc:= [
        tcn[1] + tcn[2] + tcn[3] + tcn[4],
        tcn[1] - tcn[2] + tcn[3] - tcn[4],
        tcn[1] + tcn[2] - tcn[3] - tcn[4],
        tcn[1] - tcn[2] - tcn[3] + tcn[4],
        tcn[5]
    ];
    return Htc;
end function;



function Square(tc)
    return [tc[i]^2: i in [1..5]];
end function;


function Multiply(tc1,tc2)
    return [tc1[i]*tc2[i]: i in [1..5]];
end function;





function Divide(tc_1, tc_2)
    fracs:=[[tc_1[i]*tc_2[5],tc_1[5]*tc_2[i]]:i in [1..4]];
    new_fracs:=Common_denom_frac_4(fracs);
    return [new_fracs[1][1],new_fracs[2][1],new_fracs[3][1],new_fracs[4][1],new_fracs[1][2]];
end function;



function Kappa_ii(tc_0,tc_x,tc_y)
    tc_0_Hsq := Hadamard(Square(tc_0));
    tc_x_Hsq := Hadamard(Square(tc_x));
    tc_y_Hsq := Hadamard(Square(tc_y));
    prod_tc  := Multiply(tc_x_Hsq, tc_y_Hsq);
    z0_chi := Divide(prod_tc, tc_0_Hsq);
    kappa_ii := Hadamard(z0_chi);
    kappa_ii[5]:=4*kappa_ii[5];
    return kappa_ii;
end function;



function Kappa_ii_for_double(tc_0,tc_x)
    tc_0_Hsq := Hadamard(Square(tc_0));
    tc_x_Hsq := Hadamard(Square(tc_x));
    prod_tc  := Square(tc_x_Hsq);
    z0_chi := Divide(prod_tc, tc_0_Hsq);
    kappa_ii := Hadamard(z0_chi);
    kappa_ii[5]:=4*kappa_ii[5];
    return kappa_ii;
end function;



function Double(tc_0, tc_x)
    kappa_ii := Kappa_ii_for_double(tc_0,tc_x);
    tc_2x := Divide(kappa_ii, tc_0);
    return tc_2x;
end function;




function Diff_Add(tc_0, tc_x, tc_y, tc_xmy)
    if tc_0 eq tc_xmy then
        return Double(tc_0, tc_x);
    end if;
    kappa_ii := Kappa_ii(tc_0, tc_x, tc_y);
    tc_xpy := Divide(kappa_ii, tc_xmy);
    return tc_xpy;
end function;




// Three way addition function in Magma
function Extended_Addition(tc_0, tc_x, tc_y, tc_z, tc_xpy, tc_ypz, tc_zpx,K)
    prod_0_ypz   := [tc_0[i] * tc_ypz[i]  : i in [1..4]];
    prod_zpx_xpy := [tc_zpx[i] * tc_xpy[i] : i in [1..4]];
    prod_y_z     := [tc_y[i] * tc_z[i]     : i in [1..4]];
    sum_prod_0_ypz   := [K!0 : i in [1..4]];
    sum_prod_zpx_xpy := [K!0 : i in [1..4]];
    sum_prod_y_z     := [K!0 : i in [1..4]];
    for chi in [1..4] do
        chi_t := [K!((-1)^((Floor((chi-1)/2) mod 2) * (Floor((t-1)/2) mod 2) + 
                        (chi-1) mod 2 * (t-1) mod 2)) : t in [1..4]];
        for t in [1..4] do
            sum_prod_0_ypz[chi]   +:= chi_t[t] * prod_0_ypz[t];
            sum_prod_zpx_xpy[chi] +:= chi_t[t] * prod_zpx_xpy[t];
            sum_prod_y_z[chi]     +:= chi_t[t] * prod_y_z[t];
        end for;
    end for;
    E_chi := [[sum_prod_0_ypz[chi] * sum_prod_zpx_xpy[chi], sum_prod_y_z[chi]] : chi in [1..4]];
    E_chi, den1 := Common_denom_frac_4(E_chi);
    pre_tc_xpypz := [];
    for i in [1..4] do
        chi_i := [K!((-1)^((Floor((chi-1)/2) mod 2) * (Floor((i-1)/2) mod 2) + 
                        (chi-1) mod 2 * (i-1) mod 2)) : chi in [1..4]];
        sum_part_i := K!0;
        for chi in [1..4] do
            sum_part_i +:= chi_i[chi] * E_chi[chi][1];
        end for;
        pre_tc_xpypz[i] := [sum_part_i, 4 * tc_x[i]];
    end for;
    pre_tc_xpypz, den_xpypz := Common_denom_frac_4(pre_tc_xpypz);
    den_xpypz *:= tc_xpy[5] * tc_ypz[5] * tc_zpx[5] * den1;
    dx_dy_dz := tc_x[5] * tc_y[5] * tc_z[5] * tc_0[5];
    for i in [1..4] do
        pre_tc_xpypz[i][1] *:= dx_dy_dz;
    end for;
    tc_xpypz := [];
    for i in [1..4] do
        tc_xpypz[i]:=pre_tc_xpypz[i][1];
    end for;
    tc_xpypz[5]:= den_xpypz;
    return tc_xpypz;
end function;


