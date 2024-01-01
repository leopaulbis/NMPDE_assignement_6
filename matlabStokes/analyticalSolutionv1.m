function [u] = analyticalSolutionv1(X)

[u,p] = analyticalSolution(X);
u = u(1:length(u)/2);
