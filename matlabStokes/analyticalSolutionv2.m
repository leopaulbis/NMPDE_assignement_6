function [u] = analyticalSolutionv2(X)

[u,p] = analyticalSolution(X);
u = u(length(u)/2+1:end);