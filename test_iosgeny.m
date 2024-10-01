

p:=4696592703174116824165605645274228079;

l:=5;
K:=GF(p^2);

tc_0,tc_e1,tc_e2,tc_x,tc_xpe1,tc_xpe2:=Count_prepare(p,l);

img_tc_0_one,lmd_data_1212:=CodOne(tc_0,[tc_e1,tc_e2,tc_e12],l);

img_tc_0_sq,lmd_data_1212 :=CodSq (tc_0,[tc_e1,tc_e2,tc_e12],l);
 

lmd_data_pow12:=Product_power_lambda([tc_e1,tc_e2,tc_e12],l,lmd_data_1212);


img_tc_x:=EvalOne(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data_pow12,l,K);


Time_for_isogeny(5);



//================================================


for t in [2^56..2^100] do
    p:=-1+(2^2)*3*5*7*11*13*17*19*23*29*31*37*41*43*47*53*t;
    if IsPrime(p) then
        if Log(2,p) ge 120 then
            Log(2,p);
            print(p);
            break;
        end if;
    end if;
end for;



for l in [2..55] do
    if IsPrime(l) then
        IsDivisibleBy(p+1,l);
    end if;
end for;




