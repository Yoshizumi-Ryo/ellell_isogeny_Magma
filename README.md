# Measuring the time of $(\ell,\ell)$-isogenies in Magma.

We can measure time to implement two algorithms $\mathtt{CodOne}$, $\mathtt{CodSq}$, and  $\mathtt{EvalOne}$ for $\ell=3,5,7,11,13,17,19$ over $\mathbb{F}_{p^2}$. 
Here, $p$ is 4696592703174116824165605645274228079 whose bit length is 122.


By implementing algorithms $n$ times for each, we provide the average implementation times of the algorithms. Here, $n$ is an inputed positive integer. 


## usage

First, we load ```main.m```  in Magma as follows:
```
load "main.m";
```
Then, write as follows:
```
Time_for_isogeny({the above n});
```
For example, if $n=30$,  
```
> Time_for_isogeny(30);
ell= 3
CodOne: Average time(sec): 0.001
CodSq : Average time(sec): 0.002
EvaOne: Average time(sec): 0.002

ell= 5
CodOne: Average time(sec): 0.001
CodSq : Average time(sec): 0.003
EvaOne: Average time(sec): 0.005

...

ell= 19
CodOne: Average time(sec): 0.037
CodSq : Average time(sec): 0.037
EvaOne: Average time(sec): 0.081

```


# Comparison with GetIsogeny in [NN_isogenies](https://github.com/mariascrs/NN_isogenies)

Here, we compre with the algorithm $\mathtt{GetIsogeny}$ of the paper [Isogenies on Kummer Surfaces](https://arxiv.org/abs/2409.14819).

We can implement for the same $p$ as above and for $\ell=5,7,11,13$.

Under [NN_isogenies](https://github.com/mariascrs/NN_isogenies), write in Magma as follows:
```
load "benchmarks.m";
```
Then, define the following $p$ which is the same as our implementasion:
```
p:=4696592703174116824165605645274228079;
```
For computing codomain, 
```
benchmark_getimage(100, 50, 10 : Ns := [3,5,7,11,13,17,19], primes := [p,p,p,p,p,p,p]);
```
For computing evaluation, 
```
benchmark_table1(100, 50, 10 : Ns := [5,7,11,13] , primes := [p,p,p,p]);
```




