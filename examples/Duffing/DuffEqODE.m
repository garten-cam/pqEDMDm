function dx = DuffEqODE(~,x,p)
% DuffEqODE implements the Duffing oscillator according to the array of
%parameters p.
%
% \dot(x)_1 = x_2;
% \got(x)_2 = -\delta x_2 - \alpha x_1 - \beta x_1^3
dx1 = x(2);
dx2 = - p.delta*x(2) - p.alpha*x(1) - p.beta*x(1)^3;
dx = [dx1;dx2];
end
