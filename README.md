# Measuring the time of $(\ell,\ell)$-isogenies in Magma.

We can measure time to implement two algorithms $\mathtt{CodOne}$ and $\mathtt{CodSq}$ for $\ell=3,5,7,11,13,17,19$ over $\mathbb{F}_{p^2}$. 
The bit length of $p$ is 121.

By implementing an algorithm $n$ times, we provide the average implementation time of the algorithm. Here, $n$ is an inputed positive integer. 


## usage

First, we load ```main.m```  in Magma as follows:
```
load "main.m";
```
Then, write as follows:
```
Time_for_isogeny({the above n},{"CodOne" or "CodSq"});
```
For example, if $n=30$ and for $\mathtt{CodOne}$, 
```
>Time_for_isogeny(30,"CodOne");
ell= 3
Average time(sec): 0.001

ell= 5
Average time(sec): 0.002

...

ell= 19
Average time(sec): 0.036
```
