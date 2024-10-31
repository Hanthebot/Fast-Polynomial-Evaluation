# UNIVARIATE FFT EVALUATION

## Overview
1. Compute the reduction $\bar f$ of $f$ modulo $X^p_j − X_j$ for $j = 0, \dots, m − 1$.
2. Use a fast Fourier transform to compute $\bar f(\alpha )=f(\alpha)$ for all $\alpha \in F_p^m$
3. Look up and return $f(\alpha_i)$ for all $i=0,\dots,N-1$

- Initially, we set $m=1$ (univariate) and only perform step 2 (and partially step 3) to dig into the idea.

## Interpretation / implementation of discrete FFT
1. For a field $\mathbb{F}_p$, find prime number $p'$ of form $2^n+1$ s.t. $p\leq p'$
   - as we consider 32-bit (unsigned) integer for now, let $p'=PMOD=2^{16}+1$
     - generally: it's form of $2^{2^n}$
   - it must be $2^n+1$ as it should have $2^n$ elements in field.
     - which can be iterated by RoU: 3, 952 or 59355 for example
2. Here, univariate polynomial $P(x)=\sum_{i=0}^{d\leq 2^{31}}a_ix^i$ can be expressed as
   - $P(x)=<a_0,a_1,\dots,a_d>$ (coefficient representation)
   - $P(x)=P_e(x^2)+x\cdot P_o(x^2)$
   - $P(-x)=P_e(x^2)-x\cdot P_o(x^2)$
     - where $P_e(x)=<a_0,a_2,\dots,a_{2i},\dots, a_d>$ assume $2|d$, just for convenience
     - and $P_o(x)=<a_1,a_3,\dots,a_{2i-1},\dots, a_{d-1}>$
     - $x+(-x)=PMOD$

## Challenges
1. ### Univariate evaluation
   - divide and conquer: how to implement the base case?
2. ### Multivariate evaluation
   - coefficient: received in following format
     - $a_{i_1,i_2,\dots,i_k}$: coefficient of $\prod_{j=1} X_j^k$
     - each $i_j<d$
   - store $a$ in $k$-dimension array (size: $d^k$)
     - where $A[i_1][i_2]\dots[i_k] = a_{i_1,i_2,\dots,i_k}$

### Remarks
- N/A

### Reference material
1. (Video by Reducible)[https://www.youtube.com/watch?v=h7apO7q16V0]
2. [Stack Overflow](https://mathoverflow.net/questions/115560/primitive-kth-root-of-unity-in-a-finite-field-mathbbf-p)