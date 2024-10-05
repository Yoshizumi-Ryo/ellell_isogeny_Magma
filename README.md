# Measuring the time of $(\ell,\ell)$-isogenies in Magma.

We can measure time to implement two algorithms $\mathtt{CodOne}$, $\mathtt{CodSq}$.
In addition, we can mesure the total time to compute the theta-null point of codomain and theta coordiantes of imeges of $n$ points.

Here, $p$ is 102282731615980594196068022554337214440490321465300216431359430531731842529319 whose bit length is 256, 
and we compute $(\ell,\ell)$-isogeny between Kummer surfaces over $\mathbb{F}_{p^2}$ for $\ell=5,7,11,13$.


By implementing algorithms $s$ times for each, we provide the average implementation times of the algorithms. Here, $s$ is an inputed positive integer. 

## usage

First, we load ```main.m```  in Magma as follows:
```
load "main.m";
```

### Time to compute a theta-null point of codomain.

Write as follows:
```
Time_for_isogeny_1({sample numbers s});
```
For example, if $s=10$,  
```
> Time_for_isogeny(10);
ell= 5
CodOne: Average time(sec): 0.002
CodSq : Average time(sec): 0.002

ell= 7
CodOne: Average time(sec): 0.003
CodSq : Average time(sec): 0.008

ell= 11
CodOne: Average time(sec): 0.011
CodSq : Average time(sec): 0.012

ell= 13
CodOne: Average time(sec): 0.013
CodSq : Average time(sec): 0.013
```
### Time to compute a theta-null point of codomain and theta coordinates of images of n points.

Time_for_isogeny_2(10,2);


Write as follows:
```
Time_for_isogeny_2({sample numbers s},{the number of points n});
```

For example, if $s=10$ and $n=2$, 
```
>Time_for_isogeny_2(10,2);
log_2(p)= 256
ell= 5
Codomain CodOne
Average time(sec): 0.011

log_2(p)= 256
ell= 7
Codomain CodSq
Average time(sec): 0.029

log_2(p)= 256
ell= 11
Codomain CodSq
Average time(sec): 0.065

log_2(p)= 256
ell= 13
Codomain CodSq
Average time(sec): 0.085
```



# Comparison with  [NN_isogenies](https://github.com/mariascrs/NN_isogenies)

Here, we compre with the algorithm of the paper [Isogenies on Kummer Surfaces](https://arxiv.org/abs/2409.14819).

We can implement for the same $p$ as above and for $\ell=5,7,11,13$.

Under [NN_isogenies](https://github.com/mariascrs/NN_isogenies), write in Magma as follows:
```
load "benchmarks.m";
```
Then, define the following $p$ which is the same as our implementasion:
```
p:=4696592703174116824165605645274228079;
```
For computing codomain, i.e., for measuring time of $\mathtt{GetImage}$, write
```
benchmark_getimage(100, 50, 10 : Ns := [3,5,7,11,13,17,19], primes := [p,p,p,p,p,p,p]);
```
For computing isogeny,  i.e., for measuring time of $\mathtt{GetIsogeny}$, write
```
benchmark_table1(100, 50, 10 : Ns := [5,7,11,13] , primes := [p,p,p,p]);
```
For computing evaluation.,  i.e., for measuring time of $\mathtt{Evaluation}$, 
we use 
```
load "runner.m";
```
You need to rewrite a prime number $p$ and degree $N$ in ```"runner.m"```.




