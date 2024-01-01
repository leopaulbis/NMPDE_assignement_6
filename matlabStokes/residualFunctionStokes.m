function reducedResidual=residualFunctionStokes(x,A,b,prescribedValues,prescribedDOF,unknowns)

nDOF=size(A,1);
z = zeros(nDOF,1);
z(unknowns)=x;
z(prescribedDOF)=prescribedValues;

residual = A*z-b;

reducedResidual = residual(unknowns);


