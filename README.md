# MULTIVARIATE FFT EVALUATION

## Overview
1. Compute the reduction $\bar f$ of $f$ modulo $X^p_j − X_j$ for $j = 0, \dots, m − 1$.
2. Use a fast Fourier transform to compute $\bar f(\alpha )=f(\alpha)$ for all $\alpha \in F_p^m$
3. Look up and return $f(\alpha_i)$ for all $i=0,\dots,N-1$

- Initially, we set $m=1$ (univariate) and only perform step 2 (and partially step 3) to dig into the idea.

## Interpretation / implementation of discrete FFT
1. For a field $\mathbb{F}_p$, find prime number $p$ of form $2^n+1$ s.t. $p$
   - while it must be prime, it may not need to be [Fermat prime](https://en.wikipedia.org/wiki/Fermat_number)
     - to be studied, but it might only mean extra base cases in recursion
   - the field contains $2^n$ elements in field, excluding 0.
     - which can be iterated by RoU: 3, 952 or 59355 for example
2. Here, univariate polynomial $P(x)=\sum_{i=0}^{d\leq 2^{31}}a_ix^i$ can be expressed as
   - $P(x)=<a_0,a_1,\dots,a_d>$ (coefficient representation)
   - $P(x)=P_e(x^2)+x\cdot P_o(x^2)$
   - $P(-x)=P_e(x^2)-x\cdot P_o(x^2)$
     - where $P_e(x)=<a_0,a_2,\dots,a_{2i},\dots, a_d>$ assume $2|d$, just for convenience
     - and $P_o(x)=<a_1,a_3,\dots,a_{2i-1},\dots, a_{d-1}>$
     - $x+(-x)=p$

## File structure
```
.
├── nd_vector              # Numpy's ndarray-like wrapper over memory
├── galois                 # modified version of Saied's GaloisCPP
├── fft_finite.h           # header file
├── fft_finite.cpp         # includes main()
├── fft.cpp                # FFT implementation
├── util.cpp               # initialization / IO functions
├── sample_generator.py    # generates sample input for 1D FFT
├── fft_multivar.h         # header file
├── fft_multivar.cpp       # includes main()
├── nd_fft.cpp             # FFT implementation
├── util_multivar.cpp      # initialization / IO functions
├── sample_generator_nd.py # generates sample input for ND FFT
└── ...
```

## To-Do's
1. ### Multivariate evaluation
   - coefficient: received in following format
     - $a_{i_1,i_2,\dots,i_k}$: coefficient of $\prod_{j=1} X_j^k$
     - each $i_j<d$
   - store $a$ in $n$-dimension array (size: $d^m$)
     - where $A[i_1][i_2]\dots[i_k] = a_{i_1,i_2,\dots,i_k}$
   - after operation: store $A[i_1][i_2]\dots[i_k]=f(X_1=w^{i_1},X_2=w^{i_2},\dots,X_k=w^{i_k})$
2. ### Implementing $n$-d FFT
   - implementation idea: run FFT over $n$-d array, considering each $(n-1)$-d array as a single element
   - consider implementation of: temp storage, operation overload (addition, assignment, etc.)
     - define a dedicated addition, subtraction, swap operation
       - without requiring too much extra storage & work, etc.
       - multiplication by $F$ (same type as nd vector element), member wise addition, etc. to be implemented
     - for temporary storage: assign two $(n-1)$ degree array at wrapper, pass it by reference

### Remarks
- [01/11/2024] Performs 1 million multiplication $\approx$ 0.7 seconds
  - Lenovo Yoga 6, AMD Ryzen 7 5700U
- [21/11/2024] Performs 1 million multiplication $\approx$ 0.52 seconds
  - Lenovo Yoga 6, AMD Ryzen 7 5700U
  - changed C++ version to `C++20`, for the use of `std::span`

### Reference material (univariate)
1. [Video by Reducible](https://www.youtube.com/watch?v=h7apO7q16V0)
2. [Stack Overflow](https://mathoverflow.net/questions/115560/primitive-kth-root-of-unity-in-a-finite-field-mathbbf-p)

### Reference material (multivariate)
