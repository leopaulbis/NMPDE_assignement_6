% by Sonia Fernandez-Mendez LaCaN 2016
%
% 2D Stokes problem in [0,1]x[0,1] with Dirichlet or homogeneous Neumann conditions
% Triangular FEs
clear all, clc, close all
setpath

%viscosity=1; boundaryValue=@boundaryValueFunctionCavity; sourceFunction=@sourceCavity;
numMallas = 4;

for meshNum= numMallas 
    degree=2; %Q2Q1 Taylor-Hood element
    [X,T,Xp,Tp]=CreateMeshAdaptedCavityQ2Q1(40);
    figure(1), aux=plotMesh(X,T); axis equal, title('Mesh for velocity')
    hold on, plot(X(:,1),X(:,2),'o'), hold off
    figure(2), aux=plotMesh(Xp,Tp); axis equal , title('Mesh for pressure')
    hold on, plot(Xp(:,1),Xp(:,2),'o'), hold off
    nOfNodes = size(X,1); nOfNodesp = size(Xp,1);
    
    %Reference element
    referenceElement = createReferenceElementStokesQua(degree);
    viscosity=linspace(0.1,0.02,10);
    u0=zeros(2*nOfNodes,1);
    
    for i=1:length(viscosity)
        disp(viscosity(i));
        [u,p,j,error]=solve_NS(viscosity(i),X,T,Xp,Tp,referenceElement,u0);
        disp(j);
        disp(error);
        u0=u;

        ux=u(1:nOfNodes);
        uy=u(nOfNodes+1:2*nOfNodes);
        Tboundary=connectivityMatrixBoundary(T,referenceElement);
        phi=computeStreamFunction(ux,uy,X,T,Tboundary,referenceElement);
        figure(6+i), contourPlot(phi,X,T), title('Stream lines')
        hold on, plot([0,1,1,0,0],[0,0,1,1,0],'k-'), hold off
        
        %Save the plot
        % file_name="plot/plot";
        % file_name=strcat(file_name,"_",num2str(i),".png");
        % saveas(gcf,file_name);
    end 

    ux=u(1:nOfNodes);
    uy=u(nOfNodes+1:2*nOfNodes);
    %Postprocess
    %% Plot of the solution
    figure(3),clf
    quiver(X(:,1),X(:,2),ux,uy);
    title('FEM velocity'); axis equal
    pause(0.1)
 
        %Streamlines
        Tboundary=connectivityMatrixBoundary(T,referenceElement);
        phi=computeStreamFunction(ux,uy,X,T,Tboundary,referenceElement);
        figure(6), contourPlot(phi,X,T), title('Stream lines')
        hold on, plot([0,1,1,0,0],[0,0,1,1,0],'k-'), hold off
end

