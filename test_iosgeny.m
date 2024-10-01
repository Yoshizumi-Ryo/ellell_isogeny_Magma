

p:=1536977231315969305425444639443786099;
l:=5;
K:=GF(p^2);

tc_0,tc_e1,tc_e2,tc_x,tc_xpe1,tc_xpe2:=Count_prepare(p,l);

img_tc_0_one,lmd_data_1212:=CodOne(tc_0,[tc_e1,tc_e2,tc_e12],l);

img_tc_0_sq,lmd_data_1212 :=CodSq (tc_0,[tc_e1,tc_e2,tc_e12],l);
 

lmd_data_pow12:=Product_power_lambda([tc_e1,tc_e2,tc_e12],l,lmd_data_1212);


img_tc_x:=EvalOne(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data_pow12,l,K);


Time_for_isogeny(10);
