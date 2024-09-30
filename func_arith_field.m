


//give bianry representation of N with i-digits.
function Binary_Exp(N, i)
    N_list := Intseq(N, 2);
    N_tup := N_list;
    //assert #N_tup le i;
    for m in [1..i - #N_tup] do
        Append(~N_tup, 0);
    end for;
    return N_tup;
end function;





function fatoriztion_seq(Num)
  fact:=Factorization(Num);
  fact_seq:=[];
  c:=1;
  for i in {1..#fact} do
    for j in {1..fact[i][2]} do
      fact_seq[c]:=fact[i][1];
      c+:=1;
    end for;
  end for;
  return fact_seq;
end function;
    




function decomposition_to_seq(N)
  seq:=fatoriztion_seq(N);
  inv_seq:=[];
  for i in {1..#seq} do
    inv_seq[i]:=seq[#seq-i+1];
  end for;
  comp_seq:=[];
  for i in {1..(#seq-1)} do
    comp_seq[i]:=inv_seq[i+1];
  end for;
  comp_seq[#seq]:=inv_seq[1];
  return comp_seq;
end function;




//denoms=[1,2,3,4,5]
function Common_Denom(denoms)
    //assert &and([denoms[i] ne 0 : i in [1..#denoms]]);
    N := #denoms;
    d_list := Intseq(N-1, 2);
    n := #d_list;
    a := AssociativeArray();
    b := AssociativeArray();
    for m in [1..N] do
        a[Binary_Exp(m-1, n)] := denoms[m];
    end for;
    function B_val(t)
        if t eq [] then
            return 0;
        else 
            return &+[t[j] * (2^(j-1)) : j in [1..#t]];
        end if;
    end function;
    for i in [1..n] do
        d_i_list := [d_list[j+1] : j in [i..n-1]];
        //assert #d_i_list eq n - i;
        d_i_val := B_val(d_i_list);
        for e_val in [0..d_i_val] do
            e_list := Binary_Exp(e_val, n-i);
            if e_list eq d_i_list and d_list[i] eq 0 then
                a[e_list] := a[[0] cat e_list];
            else
                a[e_list] := a[[0] cat e_list] * a[[1] cat e_list];
            end if;
        end for;
    end for;
    total_product := a[[]];
    b[[0]] := a[[1]];
    b[[1]] := a[[0]];
    for i in Reverse([0..n-2]) do
        d_i_list := [d_list[j+1] : j in [i..n-1]];
        d_i_val := B_val(d_i_list);
        for e_val in [0..d_i_val] do
            e_list := Binary_Exp(e_val, n-i);
            if e_list eq d_i_list and e_list[1] eq 0 and d_i_list[1] eq 0 then
                r_e_list := e_list[2..#e_list];
                b[e_list] := b[r_e_list];
            else
                r_e_list := e_list[2..#e_list];
                c := (e_list[1] + 1) mod 2;
                c_e_list := [c] cat r_e_list;
                b[e_list] := b[r_e_list] * a[c_e_list];
            end if;
        end for;
    end for;
    b_list := [b[Binary_Exp(i, n)] : i in [0..N-1]];
    //assert total_product ne 0;
    return b_list, total_product;
end function;





function Dict_common_denom(lincom)
    keys := Keys(lincom);
    sorted_keys := SetToSequence(keys);
    k:=#lincom;
    list_val:=[lincom[sorted_keys[i]][5] : i in [1..k]];
    list_num,den:=Common_Denom(list_val);
    dict_num := AssociativeArray();
    for i in [1..k] do
        dict_num[sorted_keys[i]] := list_num[i];
    end for;
    new_lincom:=AssociativeArray();
    for key in keys do
        new_lincom[key]:=[];
        for i in [1..4] do
            new_lincom[key][i]:=dict_num[key]*lincom[key][i]; 
        end for;
        new_lincom[key][5]:=den;
    end for;
    return new_lincom, den;
end function;





function Common_denom_frac_3(fracs)
    d_12 :=fracs[1][2]*fracs[2][2];
    d_123:=d_12*fracs[3][2];
    n_1:=fracs[1][1]*fracs[2][2]*fracs[3][2];
    n_2:=fracs[2][1]*fracs[1][2]*fracs[3][2];
    n_3:=fracs[3][1]*d_12;
    return [[n_1,d_123],[n_2,d_123],[n_3,d_123]],d_123;
end function;



// e.g. fracs=[[2,1],[3,1],[4,1],[5,1]]
function Common_denom_frac_4(fracs)
    d_12:=fracs[1][2]*fracs[2][2];
    d_34:=fracs[3][2]*fracs[4][2];
    d_1234:=d_12*d_34;
    n_1:=fracs[1][1]*fracs[2][2]*d_34;
    n_2:=fracs[2][1]*fracs[1][2]*d_34;
    n_3:=fracs[3][1]*fracs[4][2]*d_12;
    n_4:=fracs[4][1]*fracs[3][2]*d_12;
    return [[n_1,d_1234],[n_2,d_1234],[n_3,d_1234],[n_4,d_1234]],d_1234;
end function;



function Multpower_straight(value, max_exp)
    pow_list := [Parent(value)!1, value];
    value_pow := value;
    for exp in [2..max_exp] do
        value_pow *:= value;
        Append(~pow_list, value_pow);
    end for;
    //assert #pow_list eq max_exp + 1;
    pow_ass:=AssociativeArray();
    for i in [1..#pow_list] do
        pow_ass[i-1]:=pow_list[i];
    end for;
    return pow_ass;
end function;



function Multpower_sq(value, l)
    max_index_index := 2 * #Intseq(l-1, 2) - 1;
    value_2_pow := [value];
    for index_index in [1..max_index_index] do
        a := value_2_pow[index_index];
        Append(~value_2_pow, a^2);
    end for;
    result_dict := AssociativeArray();
    result_dict[0^2] := 1;
    result_dict[1^2] := value;
    for k in [2..l-1] do
        bin_exp_ksq := Intseq(k^2, 2);
        pow_value := 1;
        for i in [1..#bin_exp_ksq] do
            if bin_exp_ksq[i] eq 1 then
                pow_value *:= value_2_pow[i];
            end if;
        end for;
        result_dict[k^2] := pow_value;
    end for;
    return result_dict;
end function;


function Mult_frac(tc, frac)
    new_tc:=[tc[1]*frac[1],tc[2]*frac[1],tc[3]*frac[1],tc[4]*frac[1],tc[5]*frac[2]];
    return new_tc;
end function;





function Frac_add(frac1, frac2)
    //assert #frac1 eq 2;
    //assert #frac2 eq 2;
    //assert frac1[2] ne 0;
    //assert frac2[2] ne 0;
    return [frac1[1] * frac2[2] + frac1[2] * frac2[1], frac1[2] * frac2[2]];
end function;



function Projective_Theta(frac_coord)
    //assert #frac_coord eq 4;
    n1, d1 := Explode(frac_coord[1]);
    n2, d2 := Explode(frac_coord[2]);
    n3, d3 := Explode(frac_coord[3]);
    n4, d4 := Explode(frac_coord[4]);
    d12 := d1 * d2;
    d34 := d3 * d4;
    return [n1 * d2 * d34, d1 * n2 * d34, d12 * n3 * d4, d12 * d3 * n4];
end function;


