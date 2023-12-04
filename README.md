# pqEDMDm

pqEDMD Matlab: Repository for the development of the pqEDMD algorithm in Matlab.
Original version and citation:

[Garcia-Tenorio, C.; Vande Wouwer, A. A Matlab Toolbox for Extended Dynamic Mode Decomposition Based on Orthogonal Polynomials and p-q Quasi-Norm Order Reduction. Mathematics 2022, 10, 3859. https://doi.org/10.3390/math10203859](https://www.mdpi.com/2227-7390/10/20/3859)

## Installation

Add all the files to your Matlab path

## Minimal Working Example

Consider a Duffing oscillator, where the vector field of the system is

$$
 \begin{align}
	\dot{x}_1&=x_2\\
	\dot{x}_2&=-\delta{}x_2-x_1(\beta+\alpha{}x_1^{2}).
\end{align}
$$

## TODO

The current architecture for a decomposition is not versatile for extending to another type of decomposition different from LSQ.

- [ ] Change the architecture of the code.
      -- [ ] The observable and type of decomposition should enter the code as an argument
