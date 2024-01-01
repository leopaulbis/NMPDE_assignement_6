function [u,p,i,error]=solve_NS(viscosity,X,T,Xp,Tp,referenceElement,u0)
%__Definition of boundary conditions
    boundaryValue=@boundaryValueFunctionCavity; sourceFunction=@sourceCavity;
    nOfNodes = size(X,1); nOfNodesp = size(Xp,1);
    x = X(:,1); y = X(:,2); tol=1.e-10;
    nodesCCD = find(abs(x)<tol|abs(x-1)<tol|abs(y)<tol|abs(y-1)<tol); 
    XnodesCCD = X(nodesCCD,:);  %coordinates of the nodes on the boundary
    coefficientsCCD = [nodesCCD; nodesCCD+nOfNodes];
    uCCD = boundaryValue(XnodesCCD); %boundary value %returns a (2*nOfNodes x 1) vector
    
    %__System computation
    [K,G,f]=computeSystemStokes(X,T,Xp,Tp,referenceElement,sourceFunction,viscosity);
    C=computeNSconvectionMatrix(u0,X,T,referenceElement);

    nK=size(K,1); A = spalloc(nK+nOfNodesp,nK+nOfNodesp,nnz(K)+2*nnz(G));
    A(1:nK,1:nK)=K+C; A(1:nK,nK+1:end)=G;
    A(nK+1:end,1:nK)=G';
    b = [f; zeros(nOfNodesp,1)];
 
    pend=0; 
    prescribedValues =[uCCD;pend]; prescribedDOF = [coefficientsCCD;size(A,1)];
    
    %__Imposition of Dirichlet boundary conditions (system reduction)
    unknowns= setdiff(1:2*size(X,1)+nOfNodesp,prescribedDOF); 
    
    b = b(unknowns)-A(unknowns,prescribedDOF)*prescribedValues;
    A=A(unknowns,unknowns);
    %__System solution
    
    x0=A\b;

    aux = zeros(2*size(X,1)+nOfNodesp,1);
    aux(unknowns) = x0;
    aux(prescribedDOF) = prescribedValues;
    u1= aux(1:2*nOfNodes); p1 = aux(2*nOfNodes+1:end);
    sol1=[u1;p1];

    max_it=100;
    i=0;
    error=norm(u1-u0)/norm(u0);

    while (norm(u1-u0)>10^-3*norm(u0) && i<max_it)
        % disp(norm(u1-u0));
        %disp(i);
        
        sol0=sol1;
        u0=sol0(1:2*nOfNodes);
        i=i+1;
        
        C=computeNSconvectionMatrix(u0,X,T,referenceElement);
    
        nK=size(K,1); A = spalloc(nK+nOfNodesp,nK+nOfNodesp,nnz(K)+2*nnz(G));
        A(1:nK,1:nK)=K+C; A(1:nK,nK+1:end)=G;
        A(nK+1:end,1:nK)=G';
        b = [f; zeros(nOfNodesp,1)];

        pend=0; 
        prescribedValues =[uCCD;pend]; prescribedDOF = [coefficientsCCD;size(A,1)];
    
        %__Imposition of Dirichlet boundary conditions (system reduction)
        unknowns= setdiff(1:2*size(X,1)+nOfNodesp,prescribedDOF); 
    
        b = b(unknowns)-A(unknowns,prescribedDOF)*prescribedValues;
        A=A(unknowns,unknowns);

        x1=A\b;
        %x1=fsolve(A*x0'-b);

        aux = zeros(2*size(X,1)+nOfNodesp,1);
        aux(unknowns) = x1;
        aux(prescribedDOF) = prescribedValues;
        u1 = aux(1:2*nOfNodes); p1 = aux(2*nOfNodes+1:end);

        sol1=[u1;p1];
        error=norm(u1-u0)/norm(u0);
    end

    sol=sol1;
    u=sol(1:2*nOfNodes);
    p=sol(2*nOfNodes:end);
   