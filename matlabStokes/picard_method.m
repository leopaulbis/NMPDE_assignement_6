function [xk]=picard_method(C,b,x0,tol,max_it)
xk=x0; %initial guess 
xk_1=C\b;
i=0;
while (norm(xk_1-xk)<tol*norm(xk) & i<max_it)
    i=i+1;
    xk=xk1;
    xk1=C\b;
end 