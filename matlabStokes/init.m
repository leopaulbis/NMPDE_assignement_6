% by Sonia Fernandez-Mendez LaCaN 2016
%
% 2D Stokes problem in [0,1]x[0,1] with Dirichlet or homogeneous Neumann conditions
% Triangular FEs
clear all, clc, close all
setpath

%example definition
analytical=0;
if analytical==1 %example with ANALYTICAL solution (for convergence test)
    viscosity=2; boundaryValue=@boundaryValueFunctionAnalytical; sourceFunction=@sourceAnalytical; analytical=1;
    L2errorS=[]; L2errorsP=[]; meshsizes = [];
    numMallas=1:4;
else %CAVITY flow example
    viscosity=1; boundaryValue=@boundaryValueFunctionCavity; sourceFunction=@sourceCavity; analytical=0;
    numMallas = 4;
end

for meshNum= numMallas
    degree=2; %Q2Q1 Taylor-Hood element
    [X,T] = createUniformRectangleMesh(0,degree,0,1,0,1,2.^meshNum,2.^meshNum);
    [Xp,Tp] = createUniformRectangleMesh(0,degree-1,0,1,0,1,2.^meshNum,2.^meshNum);
    figure(1), aux=plotMesh(X,T); axis equal, title('Mesh for velocity')
    hold on, plot(X(:,1),X(:,2),'o'), hold off
    figure(2), aux=plotMesh(Xp,Tp); axis equal , title('Mesh for pressure')
    hold on, plot(Xp(:,1),Xp(:,2),'o'), hold off
    nOfNodes = size(X,1); nOfNodesp = size(Xp,1);
    
    %Reference element
    referenceElement = createReferenceElementStokesQua(degree);
    
    %__Definition of boundary conditions
    x = X(:,1); y = X(:,2); tol=1.e-10;
    nodesCCD = find(abs(x)<tol|abs(x-1)<tol|abs(y)<tol|abs(y-1)<tol); %Nodes on the boundary for fluids speed
    figure(1), hold on, plot(x(nodesCCD),y(nodesCCD),'b*','LineWidth',4); hold off
    XnodesCCD = X(nodesCCD,:);  %coordinates of the nodes on the boundary
    coefficientsCCD = [nodesCCD; nodesCCD+nOfNodes];
    uCCD = boundaryValue(XnodesCCD); %boundary value %returns a (2*nOfNodes x 1) vector
    
    %__System computation
    [K,G,f]=computeSystemStokes(X,T,Xp,Tp,referenceElement,sourceFunction,viscosity);
    nK=size(K,1); A = spalloc(nK+nOfNodesp,nK+nOfNodesp,nnz(K)+2*nnz(G));
    A(1:nK,1:nK)=K; A(1:nK,nK+1:end)=G;
    A(nK+1:end,1:nK)=G';
    b = [f; zeros(nOfNodesp,1)];
    
    %__Pressure at last node = 0 or analytical (only for pure Dirichlet problems)
    if analytical==1, [kk,pend]=analyticalSolution(Xp(end,:)); else pend=0; end
    prescribedValues =[uCCD;pend]; prescribedDOF = [coefficientsCCD;size(A,1)];
    
    %__Imposition of Dirichlet boundary conditions (system reduction)
    unknowns= setdiff(1:2*size(X,1)+nOfNodesp,prescribedDOF); %actual degrees of freedom (not boundary nodes)
    
    b = b(unknowns)-A(unknowns,prescribedDOF)*prescribedValues;
    A=A(unknowns,unknowns);
    %__System solution
    perm = symrcm(A); sol=zeros(size(A,1),1); sol(perm)=A(perm,perm)\b(perm); %With reordering
    aux = zeros(2*size(X,1)+nOfNodesp,1);
    aux(unknowns) = sol;
    aux(prescribedDOF) = prescribedValues;
    u = aux(1:2*nOfNodes); p = aux(2*nOfNodes+1:end);
    ux=u(1:nOfNodes); uy=u(nOfNodes+1:end);
    
    %Postprocess
    %% Plot of solution and convergence plots
    figure(3),clf
    quiver(X(:,1),X(:,2),ux,uy);
    title('FEM velocity'); axis equal
    figure(4),clf
    plotContinuosSolutionPressure(Xp,Tp,p,referenceElement)
    colorbar, axis equal
    title('FEM pressure')
    pause(0.1)
    if ( analytical )
        meshsize = sqrt(1/size(T,1));
        L2errorp=sqrt(mean((p-analyticalSolutionp(Xp)).^2));
        L2errorsP = [L2errorsP L2errorp];
        L2error=sqrt(mean((ux-analyticalSolutionv1(X)).^2)+mean((uy-analyticalSolutionv2(X)).^2));
        L2errorS = [L2errorS , L2error];
        meshsizes = [meshsizes, meshsize];
        fprintf('  L2 error = %0.2e\n',L2error)
        if (length(numMallas)>1 & meshNum==numMallas(end))
            figure(5), plot(log10(meshsizes),log10(L2errorS),'-o',log10(meshsizes),log10(L2errorsP),'-o');
            xlabel('log_{10} h'); ylabel('log_{10} L2 error')
            legend('u','p')
            plotSlopes(log10(meshsizes),log10(L2errorS))
            plotSlopes(log10(meshsizes),log10(L2errorsP))
            fittingcoef = polyfit(log10(meshsizes(end-2:end)),log10(L2errorS(end-2:end)),1);
            fprintf('\n slope (u)= %0.1f\n',fittingcoef(1))
            fittingcoef = polyfit(log10(meshsizes(end-2:end)),log10(L2errorsP(end-2:end)),1);
            fprintf('\n slope (p)= %0.1f\n',fittingcoef(1))
        end
    else
        %Streamlines
        Tboundary=connectivityMatrixBoundary(T,referenceElement);
        phi=computeStreamFunction(ux,uy,X,T,Tboundary,referenceElement);
        figure(6), contourPlot(phi,X,T), title('Stream lines')
        hold on, plot([0,1,1,0,0],[0,0,1,1,0],'k-'), hold off
    end
end

disp(p);