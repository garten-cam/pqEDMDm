function [Dx] = DuffEqODE(t,X,P,n,u)
%DuffEqODE is the default ODE for the Koopman operator
%algorithm. It serves an an exmaple for the implementation of
%other ODEs.
Dx1 = X(2);
Dx2 = -P.delta*X(2) - P.alpha*X(1) - P.beta*X(1)^3;
Dx = [Dx1;Dx2] + n(t,X) + u(t,X);
end