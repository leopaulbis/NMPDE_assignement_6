function val = boundaryValueFunctionCavity(X)

x = X(:,1); y = X(:,2);

tol=1.e-10;
nodesone = find(abs(y-1) < tol );
val = zeros(length(x)+length(y),1);
val(nodesone) = 1;

