# Measuring the time of $(\ell,\ell)$-isogenies in Magma.

We can measure time to implement two algorithms $\mathtt{CodOne}$, $\mathtt{CodSq}$, and  $\mathtt{EvaOne}$ for $\ell=3,5,7,11,13,17,19$ over $\mathbb{F}_{p^2}$. 
The bit length of $p$ is 121.

By implementing an algorithm $n$ times, we provide the average implementation time of the algorithm. Here, $n$ is an inputed positive integer. 


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
>load "benchmarks.m";
Loading "benchmarks.m"
Loading "functions.m"
Loading "isoNN.m"
Loading "scalings.m"
Loading "partition.m"
>p:=1536977231315969305425444639443786099;
>primes_logp100 := [p,p,p,p];
>benchmark_table1(100, 50, 10 : Ns := [5, 7, 11, 13] , primes := primes_logp100);
Using pre-generated primes.

Getting results for Method 2.
Working with N = 5
Sample: 1
Sample: 2
Sample: 3

...

[ 1536977231315969305425444639443786099, 1536977231315969305425444639443786099, 
1536977231315969305425444639443786099, 1536977231315969305425444639443786099 ]
[ 121, 121, 121, 121 ]
[
    [ 0.003, 0.001, 0.000, 0.006 ],
    [ 0.016, 0.004, 0.005, 0.030 ],
    [ 0.535, 0.011, 0.090, 0.706 ],
    [ 2.859, 0.019, 0.287, 3.334 ]
]
[
    [ 0.017, 0.002, 0.008, 0.036 ],
    [ 0.542, 0.010, 0.019, 0.644 ],
    [ 2.875, 0.018, 0.025, 3.088 ]
]

```




