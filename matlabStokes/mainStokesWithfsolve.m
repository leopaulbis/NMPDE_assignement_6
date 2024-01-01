clear all, clc, close all
setpath

%Cavity flow example definition
viscosity=1000; boundaryValue=@boundaryValueFunctionCavity; sourceFunction=@sourceCavity; analytical=0;
nOfElements1d=20;
 degree=2; %Q2Q1 Taylor-Hood element
    [X,T,Xp,Tp]=CreateMeshAdaptedCavityQ2Q1(40);
    % [X,T] = createUniformRectangleMesh(0,degree,0,1,0,1,2.^meshNum,2.^meshNum);
    % [Xp,Tp] = createUniformRectangleMesh(0,degree-1,0,1,0,1,2.^meshNum,2.^meshNum);
    figure(1), aux=plotMesh(X,T); axis equal, title('Mesh for velocity')
    hold on, plot(X(:,1),X(:,2),'o'), hold off
    figure(2), aux=plotMesh(Xp,Tp); axis equal , title('Mesh for pressure')
    hold on, plot(Xp(:,1),Xp(:,2),'o'), hold off
    nOfNodes = size(X,1); nOfNodesp = size(Xp,1);
    x = X(:,1); y = X(:,2); tol=1.e-10;
    nodesCCD = find(abs(x)<tol|abs(x-1)<tol|abs(y)<tol|abs(y-1)<tol); %Nodes on the boundary for fluids speed
    figure(1), hold on, plot(x(nodesCCD),y(nodesCCD),'b*','LineWidth',4); hold off
    XnodesCCD = X(nodesCCD,:);  %coordinates of the nodes on the boundary
    coefficientsCCD = [nodesCCD; nodesCCD+nOfNodes];
    uCCD = boundaryValue(XnodesCCD); %boundary value %returns a (2*nOfNodes x 1) vector
    referenceElement = createReferenceElementStokesQua(degree);
    
    %__System computation
    [K,G,f]=computeSystemStokes(X,T,Xp,Tp,referenceElement,sourceFunction,viscosity);
    %u0=[ones(nOfNodes,1);zeros(nOfNodes,1)];
    %u0=zeros(2*size(X,1),1);%initial speed 
    u0=zeros(nOfNodes);
    C=computeNSconvectionMatrix(u0,X,T,referenceElement);
    % figure(8);
    % spy(C);

    nK=size(K,1); A = spalloc(nK+nOfNodesp,nK+nOfNodesp,nnz(K)+2*nnz(G));
    A(1:nK,1:nK)=K+C; A(1:nK,nK+1:end)=G;
    A(nK+1:end,1:nK)=G';
    b = [f; zeros(nOfNodesp,1)];

    %zeros(length(uCCD,1));
    %__Pressure at last node = 0 or analytical (only for pure Dirichlet problems) 
    pend=0; 
    prescribedValues =[uCCD;pend]; prescribedDOF = [coefficientsCCD;size(A,1)];
    
    %__Imposition of Dirichlet boundary conditions (system reduction)
    unknowns= setdiff(1:2*size(X,1)+nOfNodesp,prescribedDOF); %actual degrees of freedom (not boundary nodes)
    
    b = b(unknowns)-A(unknowns,prescribedDOF)*prescribedValues;
    A=A(unknowns,unknowns);
    %__System solution
    
    %sol=A\b;
    
    x0=A\b;
    %x0=fsolve(A*u0'-b,u0);

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
        disp(error);
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
        unknowns= setdiff(1:2*size(X,1)+nOfNodesp,prescribedDOF); %actual degrees of freedom (not boundary nodes)
    
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

    sol=sol1;
    u=sol(1:2*nOfNodes);
    p=sol(2*nOfNodes:end);
    %disp(p);
    ux=u(1:nOfNodes); uy=u(nOfNodes+1:end);
%Postprocess
%% Plot of solution
ux=u(1:nOfNodes); uy=u(nOfNodes+1:2*nOfNodes);
figure(2),clf
quiver(X(:,1),X(:,2),ux,uy);
title('FEM velocity'); axis([0,1,0,1]); axis square
figure(3),clf
plotContinuosSolutionPressure(Xp,Tp,p,referenceElement)
colorbar, axis equal
title('FEM pressure')

%Streamlines
Tboundary=connectivityMatrixBoundary(T,referenceElement);
phi=computeStreamFunction(ux,uy,X,T,Tboundary,referenceElement);
figure(6), contourPlot(phi,X,T), title('Stream lines')
hold on, plot([0,1,1,0,0],[0,0,1,1,0],'k-'), hold off


