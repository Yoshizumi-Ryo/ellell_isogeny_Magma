# Measuring the time of $(\ell,\ell)$-isogenies in Magma.


Here algorithms are proposed in a paper: "Efficient theta-based algorithms for computing $(\ell,\ell)$-isogenies on Kummer surfaces for arbitrary odd $\ell$"". The ePrint is [here](https://eprint.iacr.org/2024/1519). 
These are written in [Magma](http://magma.maths.usyd.edu.au/magma/).


We can measure time to implement two algorithms $\mathtt{CodOne}$, $\mathtt{CodSq}$.
In addition, we can mesure the total time to compute the theta-null point of codomain and theta coordiantes of imeges of $n$ points.

Here, $p$ is 102282731615980594196068022554337214440490321465300216431359430531731842529319 whose bit length is 256, 
and we compute $(\ell,\ell)$-isogeny between Kummer surfaces over $\mathbb{F}_{p^2}$ for $\ell=5,7,11,13$.


By implementing algorithms $s$ times for each, we provide the average implementation times of the algorithms. Here, $1\le n\le 5$. 

## usage

First, we load ```main.m```  in Magma as follows:
```
load "main.m";
```

### Time to compute a theta-null point of codomain.

Write as follows:
```
Time_for_isogeny_1({the number of samples s});
```
For example, if $s=100$,  
```
> Time_for_isogeny_1(100);
log_2(p)= 256
ell= 5
Samples: 100
CodSq : Average time(sec): 0.003
CodOne: Average time(sec): 0.002

log_2(p)= 256
ell= 7
Samples: 100
CodSq : Average time(sec): 0.005
CodOne: Average time(sec): 0.004

log_2(p)= 256
ell= 11
Samples: 100
CodSq : Average time(sec): 0.012
CodOne: Average time(sec): 0.011

log_2(p)= 256
ell= 13
Samples: 100
CodSq : Average time(sec): 0.015
CodOne: Average time(sec): 0.016
```
### Time to compute a theta-null point of codomain and theta coordinates of images of n points.

Write as follows:
```
Time_for_isogeny_2({the number of samples s});
```

For example, if $s=100$, 
```
> Time_for_isogeny_2(100); 
log_2(p)= 256
ell= 5
Samples: 100
Codomain: CodOne
3 points, Average time(sec): 0.016
6 points, Average time(sec): 0.032
9 points, Average time(sec): 0.046
12 points, Average time(sec): 0.061

log_2(p)= 256
ell= 7
Samples: 100
Codomain: CodOne
3 points, Average time(sec): 0.034
6 points, Average time(sec): 0.064
9 points, Average time(sec): 0.092
12 points, Average time(sec): 0.122

log_2(p)= 256
ell= 11
Samples: 100
Codomain: CodOne
3 points, Average time(sec): 0.085
6 points, Average time(sec): 0.156
9 points, Average time(sec): 0.229
12 points, Average time(sec): 0.301

log_2(p)= 256
ell= 13
Samples: 100
Codomain: CodSq
3 points, Average time(sec): 0.118
6 points, Average time(sec): 0.219
9 points, Average time(sec): 0.321
12 points, Average time(sec): 0.422
```



# Comparison with  [NN_isogenies](https://github.com/mariascrs/NN_isogenies)

Here, we compre with the algorithm of the paper [Isogenies on Kummer Surfaces](https://arxiv.org/abs/2409.14819).

We can implement for the same $p$ as above and for $\ell=5,7,11,13$.

Under [NN_isogenies](https://github.com/mariascrs/NN_isogenies)  in Magma, load  ```additional_file.m``` in this folder: 
```
load "additional_file.m";
```

You choose the degree $\ell$, the number of samples $s$, the method to compute isogeny from GE(2) and sqrt(3).
For Table 1, write as follows: 

```
Time_SFalg_1({degree ell},{the number of samples s},{method 2(GE) or 3(sqrt)});
```
For example, if $\ell=5, s=50$, method 2(GE), then
```
> Time_SFalg_1(5,50,2);
log_2(p)= 256
ell= 5
Samples: 50
method: 2
Average time(sec): 0.006
```

For Table 2, write as follows: 
```
Time_SFalg_2({degree ell},{the number of samples s},{method 2(GE) or 3(sqrt)});
```

For example, if $\ell=7, s=30$, and method 2(GE), then

```
> Time_SFalg_2(7,30,2);
log_2(p)= 256
ell= 7
Samples: 30
method: 2
2 points, Average time(sec): 0.030
3 points, Average time(sec): 0.031
4 points, Average time(sec): 0.031
5 points, Average time(sec): 0.033

```
