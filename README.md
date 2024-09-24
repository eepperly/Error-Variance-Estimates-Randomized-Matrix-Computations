# Efficient error and variance estimation for randomized matrix computations

This repository contains code for the paper [Efficient error and variance estimation for randomized matrix computations](https://doi.org/10.1137/23M1558537) by [Ethan N. Epperly](https://www.ethanepperly.com) and [Joel A. Tropp](https://tropp.caltech.edu).

## Leave-one-out error estimation

Leave-one-out error estimation is a method for estimating the ([Frobenius norm](https://mathworld.wolfram.com/FrobeniusNorm.html)) error of a randomized low-rank approximation, such as those produced by the randomized SVD or randomized Nyström approximation.
To call use the leave-one-out error estimator with the randomized SVD, write

```matlab
[U, S, V, est] = randsvd(A, s[, q]);
```

The output of this procedure is an approximate [truncated singular value decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition#Truncated_SVD) of rank `s`, $A \approx USV^*$, and an estimate `est` of the Frobenius norm error $\\|A - USV^\*\\|_{\\rm F}$. The (optional) parameter `q` is a nonnegative integer setting the number of steps of subspace iteration; `q` is set to zero by default.

For a [(symmetric) positive semidefinite matrix](https://en.wikipedia.org/wiki/Definite_matrix) $A$, one can also use the Nyström approximation:

```matlab
[V,D,est] = nystrom(A, s[, q]);
```

The output of this procedure is an approximate truncated eigenvalue decomposition of rank `s`, $A \approx VDV^*$, and an estimate `est` of the Frobenius norm error $\\|A - VDV^\*\\|_{\\rm F}$. The (optional) parameter `q` again sets the number of steps of subspace iteration, zero by default.

**Warning:** For matrices with rapid singular value decay, these procedures can become unstable when combined with subspace iteration.
We are unaware of a way of stabilizing these procedures.

## Jackknife variance estimation

Often, the low-rank approximation $\hat{A} = USV^\*$ or $\hat{A} = VDV^\*$ is not the most important output of a randomized low-rank approximation.
Rather, we are interested in derived quantities such as eigenvectors (or, closely related, spectral projectors $V(:,1:r)V(:,1:r)^*$).
For these purposes, one can use the jackknife variance estimate, which provides an estimate of how sensitive these derived quantities are to the randomness used in the algorithm.
Variance estimates are not the same as error estimates, but a variance estimate can be a useful proxy for an error estimate in a pinch.

### Randomized Nyström approximation and positive semidefinite matrices

For sake of generality, the jackknife variance estimate is designed to work with a general quantity of interest $Q(\hat{A})$ derived from the low-rank aopproximation $\hat{A}$.
For technical reasons, we require the function $Q(\cdot)$ to be defined for matrices of all sizes and to be unitarily covariant: for a matrix $U$ with orthonormal columns, $Q(UMU^\*) = UQ(M)U^\*$.
We represent $Q(\cdot)$ in MATLAB by a function `transform` such where $Q(M) = $`transform(U,D)` with `[U,D] = eig(M)`.
The jackknife variance estimate takes the form

```matlab
[V,D,est] = nystrom(A, s, transform[, q]);
```

The output `est` is an estimate of the _standard deviation_ of $Q(\hat{A})$, defined as

$$
\\mathrm{SD} = \left( \mathbb{E} \left[\\|Q(\hat{A}) - \mathbb{E}[Q(\hat{A})]\\|_{\\rm F}^2\right] \right)^{1/2},
$$

where $\hat{A}$ is a rank-$s$ Nyström approximation.

### Randomized SVD and rectangular matrices

The interface for jackknife variance estimation is similar for the randomized SVD.
Here, we require the quantity of interest $Q(\cdot)$ be unitarily covariant on both the left and right $Q(UMV^\*) = UQ(M)V\^*$, and $Q$ is represented in MATLAB as a function $Q(M)=$`transform(U,S,V)` where `[U,S,V] = svd(M)`.
The jackknife variance estimate takes the form

```matlab
[U,S,V,est] = randsvd(A, s, transform[, q]);
```

Again, we emphasize that `est` is an estimate of the _standard deviation_ of $Q(\hat{A})$, as defined above.

## Spectral clustering

As an application of matrix jackknife variance estimation, we can use it to diagnose the quality of Nyström-accelerated spectral clustering. 
To do so, use command

```matlab
[clusters, jack, coords, V, D] = spectral_clustering(pts, kernel, num_clusters, s[, q]);
```

The `clusters` output of this procedure is a assignment of each row of the array `pts` into a cluster.
The number of clusters is set by `num_clusters`,
The `kernel` function is a positive semidefinite kernel function measuring similarity; if a numeric value is provided for `kernel`, this will be set to the _bandwidth_ for a Gaussian kernel.
This procedure uses Nyström approximation to accelerate the spectral clustering process, with `s` setting the rank of the Nyström approximation and `q` setting the number of subspace iteration steps (set to _three_ by default).
This procedure also outputs the coordinates (`coords`) outputted by spectral clustering, as well as `V` and `D` defining the Nyström approximation as above.

This procedure outputs an estimate `jack` of the standard deviation of the coordinates used for spectral clustering; see [our paper](https://doi.org/10.1137/23M1558537) for a precise definition of how this quantity is defined.
A value of `jack` around one indicates that the coordinates are highly sensitive to the randomness use in the procedure, where smaller values (e.g., `jack`$\le 0.1$ or lower) indicate that the coordinates are somewhat stable under the randomness used in the algorithm.
As we demonstrate in [the paper](https://doi.org/10.1137/23M1558537), there is a meaningful correlation between the `jack` parameter decreasing and the clustering procedure identifying the "correct" clusters.

**Warning:** There are multiple different nonequivalent spectral clustering algorithms.
See [our paper](https://doi.org/10.1137/23M1558537) for details on the procedure we use.
