# pqEDMDm

pqEDMD Matlab: Repository for the development of the pqEDMD algorithm in Matlab.
Original version and citation:

[Garcia-Tenorio, C.; Vande Wouwer, A. A Matlab Toolbox for Extended Dynamic Mode Decomposition Based on Orthogonal Polynomials and p-q Quasi-Norm Order Reduction. Mathematics 2022, 10, 3859. https://doi.org/10.3390/math10203859](https://www.mdpi.com/2227-7390/10/20/3859)

## Installation

1. Clone this repository
```bash
git clone https://github.com/garten-cam/pqEDMDm.git
``` 

2. Add the directory to your path
```matlab
>> addpath('~/pqEDMDm')
```

## Minimal Working Example

Consider a Duffing oscillator, where the differential equation of the system is,

$$
\begin{align}
	\dot{x}_1(t)&=x_2(t)\\
	\dot{x}_2(t)&=-\delta{}x_2(t)-x_1(t)(\beta+\alpha{}x_1^{2}(t)). \\
y(t) &=x(t)
\end{align}
$$

After integrating the system from six arbitrary initial conditions, and selecting two of those trajectories as the training set, and the remaining for testing, the result of plotting the states against each other is:

![Sample trajectories](examples/figures/tr_ts.png)

From the sample trajectories, apply the algorithm with a `legendreObservable` and a `pqDecomposition`, that implements a traditional
least squares regression. The result of the approximation is:

![Results](examples/figures/approx.png)

## Second Example

Consider the same Duffing oscillator with a forcing signal $u=\cos(\omega t)$ on the second state, i.e.,

$$
\begin{align}
	\dot{x}_1(t)&=x_2(t)\\
	\dot{x}_2(t)&=-\delta{}x_2(t)-x_1(t)(\beta+\alpha{}x_1^{2}(t)) + \cos(\omega t). \\
y(t) &=x(t)
\end{align}
$$

After selecting two trajectories for the estimating the system, four to test the result, and performing the approximation with a `legendreObservable` and a `svdDecomposition` gives the following results:

![Forced Duffing](examples/figures/forced_duff.png)

## TODO

The current architecture for a decomposition is not versatile for extending to another type of decomposition different from LSQ.

- [x] Change the architecture of the code.
- [x] The observable and type of decomposition should enter the pqEDMD as an argument and become an attribute (property)
- [x] Remove the huge p matrix functionality.
- [ ] Finish the sidDecomposition.
- [ ] Complete the documentation.
