


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



function Chi_t(chi,t)
    if chi eq 2 then
        if t eq 2 or t eq 4 then
            return -1;
        end if;
    end if;
    if chi eq 3 then
        if t eq 3 or t eq 4 then
            return -1;
        end if;
    end if;
    if chi eq 4 then
        if t eq 2 or t eq 3 then
            return -1;
        end if;
    end if;
    return 1;
end function;




function IpJ(i,j)
    if i eq j then
        return 1;
    end if;
    if ([i,j] eq [1,2]) or ([i,j] eq [2,1]) or  ([i,j] eq [3,4]) or ([i,j] eq [4,3]) then
        return 2;
    end if;
    if ([i,j] eq [1,3]) or ([i,j] eq [3,1]) or  ([i,j] eq [2,4]) or ([i,j] eq [4,2]) then
        return 3;
    end if;
    return 4;
end function;





function U_ichisq(tc,i,chi)
    result:=0;
    for t in [1..4] do
        result+:=Chi_t(chi,t)*tc[IpJ(i,t)]*tc[t];
    end for;
    return result;
end function;




function Z_i_chi(tc_0,tc_x,tc_y,i,chi)
    //assert(#tc_0 eq 4);
    //assert(#tc_x eq 4);
    //assert(#tc_y eq 4);
    tc_0_U := U_ichisq(tc_0,i,chi);
    tc_x_U := U_ichisq(tc_x,i,chi);
    tc_y_U := U_ichisq(tc_y,i,chi);
    return (tc_x_U*tc_y_U)/tc_0_U;
end function;




function Kappa_ij(tc_0,tc_x,tc_y)
    //assert(#tc_0 eq 4);
    //assert(#tc_x eq 4);
    //assert(#tc_y eq 4);
    kappa:=[];
    for i in [1..4] do
        kappa[i]:=[];
        for j in [1..4] do
            xx:=0;
            for chi in [1..4] do
                if Chi_t(chi,i) eq Chi_t(chi,j) then
                    assert(Chi_t(chi,IpJ(i,j)) eq 1);
                    xx+:=(1/4)*Chi_t(chi,i)*Z_i_chi(tc_0,tc_x,tc_y,IpJ(i,j),chi);
                end if;
            end for;
            kappa[i][j]:=xx;
        end for;
    end for;
    return kappa;
end function;





function Take_sqrt(a)
    _<x>:=PolynomialRing(Parent(a));
    sqrt_a:=RootsInSplittingField(x^2-a)[1][1];
    assert(sqrt_a^2 eq a);
    sqrt_a:=Parent(a)!sqrt_a;
    return sqrt_a;
end function;




function Normal_add(tc_0,tc_x,tc_y)
    //assert(#tc_0 eq 4);
    //assert(#tc_x eq 4);
    //assert(#tc_y eq 4);
    kappa:=Kappa_ij(tc_0,tc_x,tc_y);
    tc_xpy:=[];
    tc_xpy[1]:=Parent(tc_0[1])!1;
    D_2:=kappa[1][2]^2-kappa[1][1]*kappa[2][2];
    sq_D2:=Take_sqrt(D_2);
    assert(sq_D2^2 eq D_2);
    assert(Parent(sq_D2) eq Parent(D_2));
    tc_xpy[2]:=((kappa[1][2]+sq_D2)/kappa[1][1]);
    tc_xpy[3]:=(tc_xpy[2]*kappa[3][1]-kappa[3][2])/sq_D2;
    tc_xpy[4]:=(tc_xpy[2]*kappa[4][1]-kappa[4][2])/sq_D2;
    tc_xpy[5]:=1;
    return tc_xpy;
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
    tc_xpy   := Divide(kappa_ii, tc_xmy);
    return tc_xpy;
end function;


function Mult(tc_0,tc_x,k)
    tc_mx:=[];
    tc_mx[1]:=tc_x;
    tc_mx[2]:=Double(tc_0, tc_x);
    for m in [3..k] do
        tc_mx[m]:=Diff_Add(tc_0,tc_mx[m-1],tc_x,tc_mx[m-2]);
    end for;
    return tc_mx[k];
end function;



function Mxpy(tc_0,k,tc_x,tc_y,tc_xpy)
    tc_mxpy:=[];
    tc_mxpy[1]:=tc_xpy;
    tc_mxpy[2]:=Diff_Add(tc_0,tc_mxpy[1],tc_x,tc_0);
    for m in [3..k] do
        tc_mxpy[m]:=Diff_Add(tc_0,tc_mxpy[m-1],tc_x,tc_mxpy[m-2]);
    end for;
    return tc_mxpy[k];
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


