


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




