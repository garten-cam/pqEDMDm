# pqEDMDm

pqEDMD Matlab: Repository for the development of the pqEDMD algorithm in Matlab.
Original version and citation:

[Garcia-Tenorio, C.; Vande Wouwer, A. A Matlab Toolbox for Extended Dynamic Mode Decomposition Based on Orthogonal Polynomials and p-q Quasi-Norm Order Reduction. Mathematics 2022, 10, 3859. https://doi.org/10.3390/math10203859](https://www.mdpi.com/2227-7390/10/20/3859)

The pqEDMD is an evolution of the extended dynamic mode decomposition (EDMD) algorithm. The premise is to use a p-q-quasi norm reduction method to select the orders of the polynomials. From the orders, the products of univariate polynomials create the set of _observables_ that _extend_ the measurements of the system, and provide the possibility of having an operator on this set of functions. And finally, the action of the operator on the set of functions is related to the evolution of the states in a nonlinear system.

Consider an arbitrary nonlinear system $(\mathcal{M};\mathcal{U};T(x);k)$,

$$
\begin{align}
 x(k+1) &= T(x(k)) + Bu\\
   y(k) &= C x(k) + Du(k)
\end{align}
$$

where the state $x\in\mathcal{M}\subseteq\mathbb{R}^{n}$ is the nonempty compact state space, with forcing signal $u\in\mathcal{U}\subseteq\mathbb{R}^{m}$ which is a nonempty compact input space, $y\in\mathcal{Y}\subseteq\mathbb{R}^{l}$ is the output space, $k\in\mathbb{Z}_{0}^{+}$ is the discrete time, $B\in\mathbb{R}^{n\times m}$, and the nonlinear transition operator is $T\colon{}\mathcal{M}\rightarrow{}\mathcal{M}$

## Installation

Add all the files to your Matlab path

## Minimal Working Example

Consider a Duffing oscillator, where the differential equation of the system is,

$$
\begin{align}
	\dot{x}_1(t)&=x_2(t)\\
	\dot{x}_2(t)&=-\delta{}x_2(t)-x_1(t)(\beta+\alpha{}x_1^{2}(t)). \\
y(t) &=x(t)
\end{align}
$$

Notice that the output is the state and that it is an autonomous system. The objective of the pqEDMD is to get the approximation of the system in discrete-time as a linear operator of a function space.

$$
\Psi(x(k+1)) = A \Psi(x(k))
$$

$$
\Psi(x) = \begin{bmatrix}\psi_1(x) & \psi_2(x) & \cdots & \psi_l(x)\end{bmatrix}^{T}
$$

To generate the data for the approximation of the Duffing equation (in this case, where the system has two asymptotically stable points), integrate the system from some random initial conditions.

```MATLAB
% Two asymptotically stable points response
% parameters
tas.alpha = -1;
tas.beta = 1;
tas.delta = 0.5;
% The pqEDMD class accepts a structire array where the only necessary field
% in the state variables. It is not a tensor, because not all the
% trajectories are of the same lenght.

% preallocate the structure of tas orbits
tas_o = repmat(struct('sv', zeros(n_points, 2), ...
                            't', zeros(n_points, 1)), num_ics,1);
% I am saving the 't' time array only for plotting purposes. The algorithm
% does not mind if that field is in there

% Loop for all initial conditions
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
    [tas_o(orb).t, tas_o(orb).sv] = ode23(@(t,x)DuffEqODE(t,x,tas),...
        0:tfin/n_points:tfin, ...
        ics(orb,:), ...
        odeSettings);
end
```

## TODO

The current architecture for a decomposition is not versatile for extending to another type of decomposition different from LSQ.

- [ ] Change the architecture of the code.

- [ ] The observable and type of decomposition should enter the pqEDMD as an argument and become an attribute (property)
