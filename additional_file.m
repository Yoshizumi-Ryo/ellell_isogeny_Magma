

load "functions.m";
load "isoNN.m";



procedure Time_SFalg(N,num_samples,method,num_points)
    assert (N ge 5);
    assert (method eq 2 or method eq 3);
    assert (num_points ge 0);
    //this p is 256 bit.
    p:=102282731615980594196068022554337214440490321465300216431359430531731842529319;
    Fp:=GF(p);
    Fp2<i>:=ExtensionField<Fp,x|x^2+1>;
    _<x>:=PolynomialRing(Fp2);
    C:=HyperellipticCurve(x^6+1);
    C:=random_walk(C);
    K,rosen:=KummerFromCurve(C);
    C:=HyperellipticCurve(x*(x-1)*(x-rosen[1])*(x-rosen[2])*(x-rosen[3]));
    J:=Jacobian(C);
    poly<X,Y,Z,T>:=PolynomialRing(Fp2,4);
    P := [X,Y,Z,T];
    times := [];
    errors := 0; // counting number of errors we run into 
    for i in [1..num_samples] do
        repeat
            repeat
                R:=((p+1) div N)*Random(J);
            until N*R eq J!0 and R ne J!0;
            repeat
                S:=((p+1) div N)*Random(J);
            until N*S eq J!0 and S notin [J!0,R,-R];
        until WeilPairing(R,S,N) eq 1;
        R:=JtoK(R,rosen,K);
        S:=JtoK(S,rosen,K);
        points:=[];
        for nn in [1..num_points] do
            point:=RandomKummerPoint(K);
            Append(~points,point);
        end for;
        assert (#(points) eq num_points);
        time_start:=Cputime();
        phi, _ := GetIsogeny(P, R, S, K, N, method : timing := true);
        // Checking if GetImage will work (i.e., not landing on a product of ECs)
        image_thetas := Evaluate(phi, K[2]);
        for nn in [1..num_points] do
            image_P := Evaluate(phi, points[nn]);
        end for;
        time_end:=Cputime(time_start);
        Append(~times,time_end);
        "Sample:", i, time_end, "(sec)";
        //"Average:", &+times/(i-errors), ".\n";
    end for;
    av_time := &+times/(num_samples);
    log2p:=Floor(Log(p) / Log(2)) + 1;
    print("");
    printf "log_2(p)= %o\n", log2p;
    printf "ell= %o\n", N;
    print "method:", method; 
    printf "Average time(sec): %o\n", av_time;
end procedure;




