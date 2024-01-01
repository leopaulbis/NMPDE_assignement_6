function N=evaluateNodalBasisQuawithoutDerivatives(XI,nodesCoord,degree)
% N=evaluateNodalBasisQuawithoutDerivatives(XI,nodesCoord,degree)
% Evaluates at XI the nodal basis of polynomials for the given degree
% with nodal coodinates at nodesCoord

% Vandermonde matrix
V=orthogonalPolynomialsQua(degree,nodesCoord);

%Orthogonal basis at XIs and change of polynomial basis
P=orthogonalPolynomialsQua(degree,XI);
N=P/V; 


