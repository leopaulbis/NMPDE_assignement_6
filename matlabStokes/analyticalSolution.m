function [u,p] = analyticalSolution(X)
global Re;
x = X(:,1); y = X(:,2);

%u=[2*y;x]; p=x-1;

u = [y.*sin(x);-0.5*cos(x).*y.^2]; p = sin(x);

