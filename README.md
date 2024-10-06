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
For example, if $s=150$,  
```
> Time_for_isogeny_1(150);
log_2(p)= 256
ell= 5
Samples: 150
CodSq : Average time(sec): 0.003
CodOne: Average time(sec): 0.002

log_2(p)= 256
ell= 7
Samples: 150
CodSq : Average time(sec): 0.006
CodOne: Average time(sec): 0.004

log_2(p)= 256
ell= 11
Samples: 150
CodSq : Average time(sec): 0.012
CodOne: Average time(sec): 0.011

log_2(p)= 256
ell= 13
Samples: 150
CodSq : Average time(sec): 0.015
CodOne: Average time(sec): 0.016
```
### Time to compute a theta-null point of codomain and theta coordinates of images of n points.

Write as follows:
```
Time_for_isogeny_2({the number of samples s});
```

For example, if $s=50$, 
```
> Time_for_isogeny_2(150); 
log_2(p)= 256
ell= 5
Samples: 150
Codomain: CodOne
1 points, Average time(sec): 0.007
2 points, Average time(sec): 0.011
3 points, Average time(sec): 0.016
4 points, Average time(sec): 0.021
5 points, Average time(sec): 0.025

log_2(p)= 256
ell= 7
Samples: 150
Codomain: CodOne
1 points, Average time(sec): 0.014
2 points, Average time(sec): 0.024
3 points, Average time(sec): 0.033
4 points, Average time(sec): 0.043
5 points, Average time(sec): 0.052

log_2(p)= 256
ell= 11
Samples: 150
Codomain: CodOne
1 points, Average time(sec): 0.036
2 points, Average time(sec): 0.060
3 points, Average time(sec): 0.085
4 points, Average time(sec): 0.108
5 points, Average time(sec): 0.134

log_2(p)= 256
ell= 13
Samples: 150
Codomain: CodSq
1 points, Average time(sec): 0.050
2 points, Average time(sec): 0.084
3 points, Average time(sec): 0.119
4 points, Average time(sec): 0.153
5 points, Average time(sec): 0.188
```



# Comparison with  [NN_isogenies](https://github.com/mariascrs/NN_isogenies)

Here, we compre with the algorithm of the paper [Isogenies on Kummer Surfaces](https://arxiv.org/abs/2409.14819).

We can implement for the same $p$ as above and for $\ell=5,7,11,13$.

Under [NN_isogenies](https://github.com/mariascrs/NN_isogenies)  in Magma, load  ```additional_file.m``` in this folder: 
```
load "additional_file.m";
```

You choose the degree $\ell$, the number of samples $s$, the method to compute isogeny from GF and sqrt, and the number of sending points $n\ge 0$. 
In Table 1, we implement $n=0$, and in Table 2, we implements $n=1,\dots,5$.
Then, write as follows: 

```
Time_SFalg({degree ell},{the number of samples s},{method 2(GE) or 3(sqrt)},{the number of points n});
```

For example, if $\ell=5, s=10$, method GE, $n=3$; then,  

```
>Time_SFalg(5,10,2,3);
Sample: 1 0.010 (sec)
Sample: 2 0.000 (sec)
Sample: 3 0.010 (sec)
Sample: 4 0.010 (sec)
Sample: 5 0.000 (sec)
Sample: 6 0.010 (sec)
Sample: 7 0.010 (sec)
Sample: 8 0.010 (sec)
Sample: 9 0.000 (sec)
Sample: 10 0.010 (sec)

log_2(p)= 256
ell= 5
method: 2
Average time(sec): 0.007
```





